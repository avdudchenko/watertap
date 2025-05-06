from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.reaktoro_utils import (
    ViablePrecipitantsBase,
    ViableReagentsBase,
    ReaktoroOptionsContainer,
)
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import (
    assert_optimal_termination,
)

from pyomo.environ import (
    Var,
    Constraint,
    value,
    Param,
    Expression,
    Reals,
    units as pyunits,
)
from pyomo.common.config import ConfigValue
from watertap.unit_models.pressure_changer import Pump
from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock

import idaes.core.util.scaling as iscale

from watertap.unit_models.stoichiometric_reactor import (
    StoichiometricReactor,
)
from reaktoro_pse.reaktoro_block import ReaktoroBlock
from collections import OrderedDict

__author__ = "Alexander Dudchenko"


class ViableReagents(ViableReagentsBase):
    def __init__(self):
        self.register_reagent(
            "Na2CO3",
            105.99 * pyunits.g / pyunits.mol,
            {"Na_+": 2, "HCO3_-": 1},
            min_dose=0.1,
            max_dose=3000,
            purity=1,
            cost=0.19,
        )
        self.register_reagent(
            "CaO",
            56.0774 * pyunits.g / pyunits.mol,
            {"Ca_2+": 1, "H2O": 1},
            min_dose=0.1,
            max_dose=3000,
            purity=1,
            cost=0.155,
        )


class ViablePrecipitants(ViablePrecipitantsBase):
    def __init__(self):
        self.register_solid(
            "Calcite",
            100.09 * pyunits.g / pyunits.mol,
            {"Ca_2+": 1, "HCO3_-": 1},
            "Ca_2+",
        )
        self.register_solid(
            "Gypsum",
            172.17 * pyunits.g / pyunits.mol,
            {"Ca_2+": 1, "SO4_2-": 1},
            "Ca_2+",
        )
        self.register_solid(
            "Brucite",
            58.3197 * pyunits.g / pyunits.mol,
            {"Mg_2+": 1, "H2O": 2},
            "Mg_2+",
        )


@declare_process_block_class("PrecipitationUnit")
class PrecipitationUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "selected_reagents",
        ConfigValue(
            default=["CaO", "Na2CO3"],
            description="List of reagents to add to reactor",
            doc="""
            This selects reagents from ViableReagents class and adds them to reactor. 
            Default available reagents are CaO, and Na2CO3. 
            To add additional reagents, pass initialized ViableReagent
            object with registered new reagents into viable_reagents config option.
            """,
        ),
    )

    CONFIG.declare(
        "viable_reagents",
        ConfigValue(
            default=None,
            description="ViableReagents class that defines all possible reagents",
            doc="""
                Should be a ViableReagents class that contains all of the selected reagents
            """,
        ),
    )
    CONFIG.declare(
        "selected_precipitants",
        ConfigValue(
            default=["Calcite"],
            description="List of precipitants that should form in the reactor",
            doc="""
            This selects precipitants to form in reactor from Viableprecipitants class and adds them to reactor. 
            Default available precipitants are Calcite, Gypsum, and Brucite.
            To add additional precipitants, pass initialized Viableprecipitants
            object with registered new precipitants into selected_precipitants config option.
            """,
        ),
    )
    CONFIG.declare(
        "viable_precipitants",
        ConfigValue(
            default=None,
            description="Viableprecipitants class that defines all possible precipitants that can form",
            doc="""
                Should be a Viableprecipitants class that defines information for all formed precipitants
            """,
        ),
    )
    CONFIG.declare(
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="To use Reaktoro-PSE for estimating amount of solids formed",
            doc="""
            If True, builds a reaktoro block and uses it to calculate how much solids form based on feed
            composition adn amount of added reagent. 
            """,
        ),
    )
    CONFIG.declare(
        "reaktoro_options",
        ConfigValue(
            default=None,
            description="User options for configuring Reaktoro-PSE provided as a dict",
            doc="""
            User can provide additional reaktoro options, or override defaults provided by ReaktoroOptionsContainer class
            """,
        ),
    )

    def build(self):
        super().build()
        if self.config.viable_reagents is None:
            self.config.viable_reagents = ViableReagents()
        if self.config.viable_precipitants is None:
            self.config.viable_precipitants = ViablePrecipitants()

        self.selected_precipitants = {
            key: self.config.viable_precipitants[key]
            for key in self.config.selected_precipitants
        }

        self.selected_reagents = {
            key: self.config.viable_reagents[key]
            for key in self.config.selected_reagents
        }

        self.precipitation_reactor = StoichiometricReactor(
            property_package=self.config.default_property_package,
            reagent=self.selected_reagents,
            precipitate=self.selected_precipitants,
        )
        self.precipitation_reactor.pH = Var(
            ["inlet", "outlet"],
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(1, 12),
        )
        if self.config.default_costing_package is not None:
            self.precipitation_reactor.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package
            )

            self.precipitation_reactor.reagent_cost = Param(
                list(self.selected_reagents.keys()),
                units=self.config.default_costing_package.base_currency / pyunits.kg,
                mutable=True,
                domain=Reals,
            )
            for reagent, reagent_config in self.selected_reagents.items():
                self.precipitation_reactor.reagent_cost[reagent].set_value(
                    reagent_config["cost"]
                )

                self.config.default_costing_package.register_flow_type(
                    f"{self.name}_reagent_{reagent}".replace(".", "_"),
                    self.precipitation_reactor.reagent_cost[reagent],
                )
                # we adjust the mass by purity, as the flow_mass_reagent
                self.config.default_costing_package.cost_flow(
                    self.precipitation_reactor.flow_mass_reagent[reagent],
                    f"{self.name}_reagent_{reagent}".replace(".", "_"),
                )
        self.add_sludge_mass_estimation()
        if self.config.add_reaktoro_chemistry:
            self.add_reaktoro_chemistry()
        else:
            self.build_equality_ph_constraints()
        self.register_port(
            "inlet",
            self.precipitation_reactor.inlet,
            {"pH": self.precipitation_reactor.pH["inlet"]},
        )
        self.register_port(
            "outlet",
            self.precipitation_reactor.outlet,
            {"pH": self.precipitation_reactor.pH["outlet"]},
        )
        self.register_port(
            "sludge",
            self.precipitation_reactor.waste,
            {"pH": self.precipitation_reactor.pH["outlet"]},
        )

    def build_equality_ph_constraints(self):
        self.eq_outlet_pH = Constraint(
            expr=self.precipitation_reactor.pH["inlet"]
            == self.precipitation_reactor.pH["outlet"]
        )

    def add_sludge_mass_estimation(self):
        sludge_components = []
        for (
            phase,
            ion,
        ), obj in self.precipitation_reactor.separator.waste_state[
            0.0
        ].flow_mass_phase_comp.items():
            sludge_components.append(obj)
        for (
            precipitants,
            obj,
        ) in self.precipitation_reactor.flow_mass_precipitate.items():
            sludge_components.append(obj)

        self.precipitation_reactor.total_sludge_product = Expression(
            expr=sum(sludge_components)
        )

    def add_reaktoro_chemistry(self):
        solvents = self.config.viable_reagents.create_solvent_constraint(
            self.precipitation_reactor, self.precipitation_reactor.flow_mol_reagent
        )
        reagents = {}
        for r in self.selected_reagents:
            reagents[r] = self.precipitation_reactor.flow_mol_reagent[r]

        if solvents is not None:
            for solvent in solvents:
                reagents[solvent] = self.precipitation_reactor.flow_mol_solvent[solvent]

        outputs = {("pH", None): self.precipitation_reactor.pH["outlet"]}
        for phase, obj in self.precipitation_reactor.flow_mol_precipitate.items():
            outputs[("speciesAmount", phase)] = obj

        self.reaktoro_options = ReaktoroOptionsContainer()
        self.reaktoro_options.system_state_option(
            "temperature",
            self.precipitation_reactor.dissolution_reactor.properties_in[0].temperature,
        )
        self.reaktoro_options.system_state_option(
            "pressure",
            self.precipitation_reactor.dissolution_reactor.properties_in[0].pressure,
        )
        self.reaktoro_options.system_state_option(
            "pH", self.precipitation_reactor.pH["inlet"]
        )
        self.reaktoro_options.aqueous_phase_option(
            "composition",
            self.precipitation_reactor.dissolution_reactor.properties_in[
                0
            ].flow_mol_phase_comp,
        )
        self.reaktoro_options["mineral_phase"] = {
            "phase_components": list(self.selected_precipitants.keys())
        }
        self.reaktoro_options["register_new_chemistry_modifiers"] = (
            self.config.viable_reagents.get_reaktoro_chemistry_modifiers()
        )
        self.reaktoro_options["chemistry_modifier"] = reagents
        self.reaktoro_options["outputs"] = outputs
        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)

        self.precipitation_block = ReaktoroBlock(**self.reaktoro_options)

    def set_fixed_operation(self):
        for reagent, options in self.selected_reagents.items():
            self.precipitation_reactor.reagent_dose[reagent].setlb(
                options["min_dose"] / 1000
            )
            self.precipitation_reactor.reagent_dose[reagent].setub(
                options["max_dose"] / 1000
            )
            self.precipitation_reactor.flow_mol_reagent[reagent].setlb(None)

        self.precipitation_reactor.waste_mass_frac_precipitate.fix(0.2)

        for precip in self.selected_precipitants.keys():
            self.precipitation_reactor.flow_mol_precipitate[precip].setlb(None)
            if self.config.add_reaktoro_chemistry == False:
                self.precipitation_reactor.flow_mol_precipitate[precip].fix(1e-5)
        for reagent, _ in self.selected_reagents.items():
            self.precipitation_reactor.reagent_dose[reagent].fix(10 / 1000)

        self.inlet.fix()
        assert degrees_of_freedom(self) == 0
        self.inlet.unfix()

    def scale_before_initialization(self, **kwargs):
        max_dose = []
        for reagent, options in self.selected_reagents.items():
            dose_scale = options["max_dose"]
            max_dose.append(dose_scale)
            # use mol flow, as thats what will be propagated by default via mcas
            mass_flow_scale = dose_scale / value(
                self.precipitation_reactor.inlet.flow_mol_phase_comp[0.0, "Liq", "H2O"]
                * (18.015 / 1000)
            )
            iscale.set_scaling_factor(
                self.precipitation_reactor.flow_mass_reagent[reagent],
                mass_flow_scale,
            )
            iscale.set_scaling_factor(
                self.precipitation_reactor.flow_mol_reagent[reagent],
                mass_flow_scale
                * value(
                    pyunits.convert(
                        self.config.viable_reagents[reagent]["mw"],
                        pyunits.kg / pyunits.mol,
                    )
                ),
            )
            iscale.set_scaling_factor(
                self.precipitation_reactor.reagent_dose[reagent], dose_scale
            )
        precip_scale = max(max_dose) / value(
            self.precipitation_reactor.inlet.flow_mol_phase_comp[0.0, "Liq", "H2O"]
            * (18.015 / 1000)
        )
        for precip in self.selected_precipitants.keys():
            iscale.set_scaling_factor(
                self.precipitation_reactor.flow_mass_precipitate[precip],
                precip_scale,
            )
        iscale.set_scaling_factor(self.precipitation_reactor.pH, 1)
        if self.config.add_reaktoro_chemistry:
            self.config.viable_reagents.scale_solvent_vars_and_constraints(
                self.precipitation_reactor, self.precipitation_reactor.flow_mol_reagent
            )
        else:
            iscale.constraint_scaling_transform(self.eq_outlet_pH, 1)

    def initialize_unit(self, **kwargs):

        for phase, data in self.selected_precipitants.items():
            # assume that only fraction of ions will actually precipitate
            flow = (
                self.precipitation_reactor.inlet.flow_mol_phase_comp[
                    0.0, "Liq", data["primary_ion"]
                ].value
                * 0.0001
            )
            self.precipitation_reactor.flow_mol_precipitate[phase].fix(flow)
        self.precipitation_reactor.initialize()
        if self.config.add_reaktoro_chemistry:
            # get initial mole flows
            self.precipitation_block.initialize()
            self.precipitation_block.display_jacobian_scaling()
            # recalculate state with updated mol flow values
            self.precipitation_reactor.initialize()
            for phase, data in self.selected_precipitants.items():
                self.precipitation_reactor.flow_mol_precipitate[phase].unfix()

    def get_model_state_dict(self):
        def get_ion_comp(stream, pH):
            data_dict = OrderedDict()
            data_dict["Mass flow of H2O"] = stream.flow_mass_phase_comp["Liq", "H2O"]
            for phase, ion in stream.conc_mass_phase_comp:
                if ion != "H2O":
                    data_dict[ion] = stream.conc_mass_phase_comp[phase, ion]
            data_dict["pH"] = pH
            data_dict["Temperature"] = stream.temperature
            data_dict["Pressure"] = stream.pressure
            return data_dict

        self.inlet.fix()
        unit_dofs = degrees_of_freedom(self)
        self.inlet.unfix()
        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "Inlet state": get_ion_comp(
                self.precipitation_reactor.dissolution_reactor.properties_in[0],
                self.precipitation_reactor.pH["inlet"],
            ),
            "Chemical dosing:": self.precipitation_reactor.reagent_dose,
            "Solids formed:": self.precipitation_reactor.flow_mass_precipitate,
            "Treated state": get_ion_comp(
                self.precipitation_reactor.precipitation_reactor.properties_out[0],
                self.precipitation_reactor.pH["outlet"],
            ),
            "Waste state": get_ion_comp(
                self.precipitation_reactor.separator.waste_state[0],
                self.precipitation_reactor.pH["outlet"],
            ),
        }
        if self.config.default_costing_package is not None:
            model_state_dict["Costs"] = {
                "Capital cost": self.precipitation_reactor.costing.capital_cost
            }
            for reagent in self.precipitation_reactor.flow_mass_reagent:
                model_state_dict[f"Costs"][f"Reagent {reagent} cost"] = (
                    self.config.default_costing_package.aggregate_flow_costs[
                        f"{self.name}_reagent_{reagent}".replace(".", "_")
                    ]
                )
        return self.name, model_state_dict

__author__ = "Alexander Dudchenko"
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.reaktoro_utils import (
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
            "HCl",
            36.46 * pyunits.g / pyunits.mol,
            {"Cl_-": 1, "H2O": 1},
            min_dose=0.1,
            max_dose=3000,
            purity=0.38,
            solvent=("H2O", 18.01 * pyunits.g / pyunits.mol),
            cost=0.17,
            density_reagent=1.18 * pyunits.kg / pyunits.liter,
        )
        self.register_reagent(
            "H2SO4",
            98.08 * pyunits.g / pyunits.mol,
            {"SO4_2-": 1, "H2O": 1},
            min_dose=0.1,
            max_dose=3000,
            purity=0.93,
            solvent=("H2O", 18.01 * pyunits.g / pyunits.mol),
            cost=0.12,
            density_reagent=1.8136 * pyunits.kg / pyunits.liter,
        )


@declare_process_block_class("ChemicalAdditionUnit")
class ChemicalAdditionUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "selected_reagents",
        ConfigValue(
            default=["HCl"],
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
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="To use Reaktoro-PSE for estimate pH change",
            doc="""
            If True, builds a reaktoro block and uses it to calculate change in pH due to the addition of given chemical. 
            """,
        ),
    )
    CONFIG.declare(
        "reaktoro_options",
        ConfigValue(
            default=None,
            description="User options for configuring Reaktoro-PSE",
            doc="""
            User can provide additional reaktoro options, or override defaults provided by ReaktoroOptionsContainer class
            """,
        ),
    )

    def build(self):
        super().build()
        if self.config.viable_reagents is None:
            self.config.viable_reagents = ViableReagents()
        self.selected_reagents = {
            key: self.config.viable_reagents[key]
            for key in self.config.selected_reagents
        }

        self.chemical_reactor = StoichiometricReactor(
            property_package=self.config.default_property_package,
            reagent=self.selected_reagents,
        )
        self.chemical_reactor.pH = Var(
            ["inlet", "outlet"],
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(1, 12),
        )
        if self.config.default_costing_package is not None:
            self.chemical_reactor.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package
            )

            self.chemical_reactor.reagent_cost = Param(
                list(self.selected_reagents.keys()),
                units=self.config.default_costing_package.base_currency / pyunits.kg,
                mutable=True,
                domain=Reals,
            )
            for reagent, reagent_config in self.selected_reagents.items():
                self.chemical_reactor.reagent_cost[reagent].set_value(
                    reagent_config["cost"]
                )

                self.config.default_costing_package.register_flow_type(
                    f"{self.name}_reagent_{reagent}".replace(".", "_"),
                    self.chemical_reactor.reagent_cost[reagent],
                )
                # we adjust the mass by purity, as the flow_mass_reagent
                self.config.default_costing_package.cost_flow(
                    self.chemical_reactor.flow_mass_reagent[reagent],
                    f"{self.name}_reagent_{reagent}".replace(".", "_"),
                )
        if self.config.add_reaktoro_chemistry:
            self.add_reaktoro_chemistry()
        else:
            self.chemical_reactor.eq_ph = Constraint(
                expr=self.chemical_reactor.pH["inlet"]
                == self.chemical_reactor.pH["outlet"]
            )
        self.register_port(
            "inlet",
            self.chemical_reactor.inlet,
            {"pH": self.chemical_reactor.pH["inlet"]},
        )
        self.register_port(
            "outlet",
            self.chemical_reactor.outlet,
            {"pH": self.chemical_reactor.pH["outlet"]},
        )

    def add_reaktoro_chemistry(self):
        solvents = self.config.viable_reagents.create_solvent_constraint(
            self.chemical_reactor, self.chemical_reactor.flow_mol_reagent
        )
        reagents = {}
        for r in self.selected_reagents:
            reagents[r] = self.chemical_reactor.flow_mol_reagent[r]
        if solvents is not None:
            for solvent in solvents:
                reagents[solvent] = self.chemical_reactor.flow_mol_solvent[solvent]

        outputs = {("pH", None): self.chemical_reactor.pH["outlet"]}

        self.reaktoro_options = ReaktoroOptionsContainer()
        self.reaktoro_options.system_state_option(
            "temperature",
            self.chemical_reactor.dissolution_reactor.properties_in[0].temperature,
        )
        self.reaktoro_options.system_state_option(
            "pressure",
            self.chemical_reactor.dissolution_reactor.properties_in[0].pressure,
        )
        self.reaktoro_options.system_state_option(
            "pH", self.chemical_reactor.pH["inlet"]
        )
        self.reaktoro_options.aqueous_phase_option(
            "composition",
            self.chemical_reactor.dissolution_reactor.properties_in[
                0
            ].flow_mol_phase_comp,
        )
        self.reaktoro_options["register_new_chemistry_modifiers"] = (
            self.config.viable_reagents.get_reaktoro_chemistry_modifiers()
        )
        self.reaktoro_options["chemistry_modifier"] = reagents
        self.reaktoro_options["outputs"] = outputs

        self.chemistry_block = ReaktoroBlock(**self.reaktoro_options)

    def set_fixed_operation(self):
        for reagent, options in self.selected_reagents.items():
            self.chemical_reactor.reagent_dose[reagent].setlb(
                options["min_dose"] / 1000
            )
            self.chemical_reactor.reagent_dose[reagent].setub(
                options["max_dose"] / 1000
            )
            self.chemical_reactor.flow_mol_reagent[reagent].setlb(None)
        # self.chemical_reactor.pH.fix()

    def scale_before_initialization(self, **kwargs):
        max_dose = []
        for reagent, options in self.selected_reagents.items():
            dose_scale = options["max_dose"]
            max_dose.append(dose_scale)
            # use mol flow, as thats what will be propagated by default via mcas
            mass_flow_scale = dose_scale / value(
                self.chemical_reactor.inlet.flow_mol_phase_comp[0.0, "Liq", "H2O"]
                * (18.015 / 1000)
            )
            iscale.set_scaling_factor(
                self.chemical_reactor.flow_mass_reagent[reagent],
                mass_flow_scale,
            )
            iscale.set_scaling_factor(
                self.chemical_reactor.reagent_dose[reagent], dose_scale
            )
        iscale.set_scaling_factor(self.chemical_reactor.pH, 1)

    def initialize_unit(self, **kwargs):

        for reagent, _ in self.selected_reagents.items():
            self.chemical_reactor.reagent_dose[reagent].fix(10 / 1000)

        self.chemical_reactor.initialize()
        if self.config.add_reaktoro_chemistry:
            # get intial mole flows
            self.chemistry_block.initialize()
            self.chemistry_block.display_jacobian_scaling()
            # recalcualte state with updated mol flow values
            self.chemical_reactor.initialize()

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

        model_state = {
            "Inlet state": get_ion_comp(
                self.chemical_reactor.dissolution_reactor.properties_in[0],
                self.chemical_reactor.pH["inlet"],
            ),
            "Chemical dosing:": self.chemical_reactor.reagent_dose,
            "Treated state": get_ion_comp(
                self.chemical_reactor.dissolution_reactor.properties_out[0],
                self.chemical_reactor.pH["outlet"],
            ),
        }
        if self.config.default_costing_package is not None:
            model_state["Costs"] = {
                "Capital cost": self.chemical_reactor.costing.capital_cost
            }
            for reagent in self.chemical_reactor.flow_mass_reagent:
                model_state[f"Costs"][f"Reagent {reagent} cost"] = (
                    self.config.default_costing_package.aggregate_flow_costs[
                        f"{self.name}_reagent_{reagent}".replace(".", "_"),
                    ]
                )
        return self.name, model_state

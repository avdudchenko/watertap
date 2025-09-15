__author__ = "Alexander Dudchenko"
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.reaktoro_utils import (
    ViableReagents,
    ReaktoroOptionsContainer,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from pyomo.environ import (
    Var,
    Constraint,
    value,
    Param,
    Reals,
    units as pyunits,
)
from pyomo.common.config import ConfigValue
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
    CONFIG.declare(
        "add_alkalinity",
        ConfigValue(
            default=True,
            description="Defines if to add alkalinity to Reaktoro output",
            doc="""
                Defines if to add alkalinity to Reaktoro output, this will use Reaktoro to calculate alkalinity
            """,
        ),
    )
    CONFIG.declare(
        "track_pE",
        ConfigValue(
            default=False,
            description="if pE should be tracked in the model",
            doc="""
                    Providing True will add pE variable to the model and track it
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
        self.chemical_reactor.dissolution_reactor.properties_in[0].conc_mass_phase_comp[
            ...
        ]
        self.chemical_reactor.dissolution_reactor.properties_out[
            0
        ].conc_mass_phase_comp[...]
        self.chemical_reactor.pH = Var(
            ["inlet", "outlet"],
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(0, 13),
        )
        if self.config.track_pE:
            self.chemical_reactor.pE = Var(
                ["inlet", "outlet"],
                initialize=7,
                units=pyunits.dimensionless,
                bounds=(None, None),
            )
        if self.config.default_costing_package is not None:
            self.chemical_reactor.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package,
                **self.config.default_costing_package_kwargs,
            )

            self.chemical_reactor.reagent_cost = Var(
                list(self.selected_reagents.keys()),
                units=self.config.default_costing_package.base_currency / pyunits.kg,
                domain=Reals,
            )
            for reagent, reagent_config in self.selected_reagents.items():
                self.chemical_reactor.reagent_cost[reagent].fix(reagent_config["cost"])

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
            if self.chemical_reactor.find_component("pE") is not None:
                self.chemical_reactor.eq_pE = Constraint(
                    expr=self.chemical_reactor.pE["inlet"]
                    == self.chemical_reactor.pE["outlet"]
                )
        output_vars = {"pH": self.chemical_reactor.pH["outlet"]}
        if self.config.track_pE:
            output_vars["pE"] = self.chemical_reactor.pE["outlet"]
        input_vars = {"pH": self.chemical_reactor.pH["inlet"]}
        if self.config.track_pE:
            input_vars["pE"] = self.chemical_reactor.pE["inlet"]

        self.register_port(
            "inlet",
            self.chemical_reactor.inlet,
            input_vars,
        )
        self.register_port(
            "outlet",
            self.chemical_reactor.outlet,
            output_vars,
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
        if self.config.track_pE:
            outputs[("pE", None)] = self.chemical_reactor.pE["outlet"]
        if self.config.add_alkalinity:
            self.chemical_reactor.alkalinity = Var(
                initialize=1,
                units=pyunits.mg / pyunits.L,
                doc="Alkalinity (mg/L as CaCO3)",
            )
            outputs[("alkalinityAsCaCO3", None)] = self.chemical_reactor.alkalinity
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
        if self.config.track_pE:
            self.reaktoro_options.system_state_option(
                "pE", self.chemical_reactor.pE["inlet"]
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

        # self.reaktoro_options["chemistry_modifier_log10_basis"] = True
        self.reaktoro_options["chemistry_modifier"] = reagents
        self.reaktoro_options["outputs"] = outputs
        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)
        self.chemistry_block = ReaktoroBlock(**self.reaktoro_options)

    def set_fixed_operation(self):
        """fixes operation point for chemical addition unit model"""
        for reagent, _ in self.selected_reagents.items():
            self.chemical_reactor.reagent_dose[reagent].fix(10 / 1000)
        for reagent, options in self.selected_reagents.items():
            self.chemical_reactor.reagent_dose[reagent].setlb(
                options["min_dose"] / 1000
            )
            self.chemical_reactor.reagent_dose[reagent].setub(
                options["max_dose"] / 1000
            )
            self.chemical_reactor.flow_mol_reagent[reagent].setlb(0)

    def set_optimization_operation(self):
        """if we have reaktoro chemistry, we need to unfix reagent addition"""
        if self.config.add_reaktoro_chemistry:
            for reagent, _ in self.selected_reagents.items():
                self.chemical_reactor.reagent_dose[reagent].unfix()

    def scale_before_initialization(self, **kwargs):
        dose_scale = 1 / 0.001  # step size of ppm (or 0.001 kg/m3)
        for reagent, options in self.selected_reagents.items():
            # use mol flow, as thats what will be propagated by default via mcas
            mass_flow_scale = dose_scale / (
                self.config.default_property_package._default_scaling_factors[
                    "flow_mol_phase_comp", ("Liq", "H2O")
                ]
                / value(self.config.default_property_package.mw_comp["H2O"])
            )
            mol_scale = mass_flow_scale * value(
                pyunits.convert(
                    self.config.viable_reagents[reagent]["mw"],
                    pyunits.kg / pyunits.mol,
                )
            )
            iscale.set_scaling_factor(
                self.chemical_reactor.flow_mol_reagent[reagent],
                mol_scale,
            )
            iscale.set_scaling_factor(
                self.chemical_reactor.flow_mass_reagent[reagent],
                mass_flow_scale,
            )
            iscale.set_scaling_factor(
                self.chemical_reactor.reagent_dose[reagent], dose_scale
            )
        if self.chemical_reactor.find_component("pE") is not None:
            iscale.set_scaling_factor(self.chemical_reactor.pE, 1 / 10)
        iscale.set_scaling_factor(self.chemical_reactor.pH, 1 / 10)
        if self.config.add_reaktoro_chemistry:
            self.config.viable_reagents.scale_solvent_vars_and_constraints(
                self.chemical_reactor
            )
            if self.config.add_alkalinity:
                iscale.set_scaling_factor(self.chemical_reactor.alkalinity, 1)

    def initialize_unit(self, **kwargs):
        self.chemical_reactor.initialize()
        if self.config.add_reaktoro_chemistry:
            self.chemistry_block.initialize()
            self.chemistry_block.display_jacobian_scaling()
            # recalcualte state with updated mol flow values
            self.chemical_reactor.initialize()

    def get_model_state_dict(self):
        def get_ion_comp(stream, pH, pE=None):
            data_dict = OrderedDict()
            data_dict["Mass flow of H2O"] = stream.flow_mass_phase_comp["Liq", "H2O"]
            for phase, ion in stream.conc_mass_phase_comp:
                if ion != "H2O":
                    data_dict[ion] = stream.conc_mass_phase_comp[phase, ion]
            data_dict["pH"] = pH
            if pE is not None:
                data_dict["pE"] = pE
            data_dict["Temperature"] = stream.temperature
            data_dict["Pressure"] = stream.pressure
            return data_dict

        def get_pe(port):
            if self.config.track_pE:
                return self.chemical_reactor.pE[port]
            else:
                return None

        self.inlet.fix()
        unit_dofs = degrees_of_freedom(self)
        self.inlet.unfix()
        treated_state = get_ion_comp(
            self.chemical_reactor.dissolution_reactor.properties_out[0],
            self.chemical_reactor.pH["outlet"],
            get_pe("outlet"),
        )
        if (
            self.config.add_alkalinity != False
            and self.chemical_reactor.find_component("alkalinity") is not None
        ):
            treated_state["Alkalinity (mg/L as CaCO3)"] = (
                self.chemical_reactor.alkalinity
            )
        model_state = {
            "Model": {"DOFs": unit_dofs},
            "Inlet state": get_ion_comp(
                self.chemical_reactor.dissolution_reactor.properties_in[0],
                self.chemical_reactor.pH["inlet"],
                get_pe("inlet"),
            ),
            "Chemical dosing:": self.chemical_reactor.reagent_dose,
            "Chemical mass flow:": self.chemical_reactor.flow_mass_reagent,
            "Treated state": treated_state,
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

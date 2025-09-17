from watertap.core.util import scaling
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.reaktoro_utils import (
    ViablePrecipitantsBase,
    ViableReagents,
    ReaktoroOptionsContainer,
)
from idaes.core.util.model_statistics import degrees_of_freedom

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
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.core.util.math as idaesMath

__author__ = "Alexander Dudchenko"


class ViablePrecipitants(ViablePrecipitantsBase):
    def __init__(self):
        self.register_solid(
            "Calcite",
            100.09 * pyunits.g / pyunits.mol,
            {"Ca_2+": 1, "HCO3_-": 1},
            "Ca_2+",
            reaktoro_modifier={"Ca": -1, "C": -1, "O": -3},
        )
        self.register_solid(
            "Gypsum",
            172.17 * pyunits.g / pyunits.mol,
            {"Ca_2+": 1, "SO4_2-": 1},
            "Ca_2+",
            reaktoro_modifier={"Ca": -1, "S": -1, "O": -4},
        )
        self.register_solid(
            "Brucite",
            58.3197 * pyunits.g / pyunits.mol,
            {"Mg_2+": 1, "H2O": 2},
            "Mg_2+",
            reaktoro_modifier={"Mg": -1, "O": -2, "H": -2},
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
        "add_non_eq_reaktoro_chemistry",
        ConfigValue(
            default=False,
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
        "add_hardness",
        ConfigValue(
            default=False,
            description="Defines if to add Calcium and Magnesium hardness precipitate output",
            doc="""
            This will calculate Calcium (if present), Magnesium (if present) and Total hardness as CaCO3.
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
        self.precipitation_reactor.precipitation_reactor.properties_out[
            0
        ].conc_mass_phase_comp[...]
        self.precipitation_reactor.precipitation_reactor.properties_out[
            0
        ].flow_mass_phase_comp[...]

        self.precipitation_reactor.dissolution_reactor.properties_in[
            0
        ].conc_mass_phase_comp[...]
        self.precipitation_reactor.dissolution_reactor.properties_in[
            0
        ].flow_mass_phase_comp[...]
        self.precipitation_reactor.dissolution_reactor.properties_out[
            0
        ].flow_mass_phase_comp[...]
        self.precipitation_reactor.dissolution_reactor.properties_out[
            0
        ].conc_mass_phase_comp[...]
        self.precipitation_reactor.separator.waste_state[0].conc_mass_phase_comp[...]
        self.precipitation_reactor.separator.waste_state[0].flow_mass_phase_comp[...]

        self.precipitation_reactor.pH = Var(
            ["inlet", "outlet"],
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(0, 13),
        )
        if self.config.track_pE:
            self.precipitation_reactor.pE = Var(
                ["inlet", "outlet"],
                initialize=0,
                units=pyunits.dimensionless,
                bounds=(None, None),
            )
        if self.config.default_costing_package is not None:
            self.precipitation_reactor.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package,
                **self.config.default_costing_package_kwargs,
            )

            self.precipitation_reactor.reagent_cost = Var(
                list(self.selected_reagents.keys()),
                units=self.config.default_costing_package.base_currency / pyunits.kg,
                # mutable=True,
                domain=Reals,
            )
            for reagent, reagent_config in self.selected_reagents.items():
                self.precipitation_reactor.reagent_cost[reagent].fix(
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
        elif self.config.add_non_eq_reaktoro_chemistry:
            self.add_non_eq_reaktoro_chemistry()
        else:
            self.build_equality_ph_pe_constraints()
        if self.config.add_hardness:
            self.add_hardness()
        inlet_vars = {"pH": self.precipitation_reactor.pH["inlet"]}
        outlet_vars = {"pH": self.precipitation_reactor.pH["outlet"]}
        if self.config.track_pE:
            inlet_vars["pE"] = self.precipitation_reactor.pE["inlet"]
            outlet_vars["pE"] = self.precipitation_reactor.pE["outlet"]
        self.register_port(
            "inlet",
            self.precipitation_reactor.inlet,
            inlet_vars,
        )
        self.register_port(
            "outlet",
            self.precipitation_reactor.outlet,
            outlet_vars,
        )
        self.register_port(
            "sludge",
            self.precipitation_reactor.waste,
            outlet_vars,
        )

    def build_equality_ph_pe_constraints(self):
        self.eq_outlet_pH = Constraint(
            expr=self.precipitation_reactor.pH["inlet"]
            == self.precipitation_reactor.pH["outlet"]
        )
        if self.config.track_pE:
            self.eq_outlet_pE = Constraint(
                expr=self.precipitation_reactor.pE["inlet"]
                == self.precipitation_reactor.pE["outlet"]
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

    def add_hardness(self):

        self.precipitation_reactor.hardness = Var(
            ["Ca", "Mg", "Total"],
            initialize=0,
            units=pyunits.mg / pyunits.L,
            doc="Hardness (mg/L as CaCO3)",
        )

        @self.precipitation_reactor.Constraint(["Ca", "Mg", "Total"])
        def eq_hardness(b, ion):
            if ion == "Ca":
                return b.hardness["Ca"] == pyunits.convert(
                    b.precipitation_reactor.properties_out[0].flow_mol_phase_comp[
                        "Liq", "Ca_2+"
                    ]
                    * (100.09 * pyunits.g / pyunits.mol)
                    / b.precipitation_reactor.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.mg / pyunits.L,
                )
            elif ion == "Mg":
                return b.hardness["Mg"] == pyunits.convert(
                    b.precipitation_reactor.properties_out[0].flow_mol_phase_comp[
                        "Liq", "Mg_2+"
                    ]
                    * (100.09 * pyunits.g / pyunits.mol)
                    / b.precipitation_reactor.properties_out[0].flow_vol_phase["Liq"],
                    to_units=pyunits.mg / pyunits.L,
                )

            elif ion == "Total":
                return b.hardness["Total"] == b.hardness["Ca"] + b.hardness["Mg"]

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
        self.rkt_block_outputs = {("pH", None): self.precipitation_reactor.pH["outlet"]}
        if self.config.track_pE:
            self.rkt_block_outputs[("pE", None)] = self.precipitation_reactor.pE[
                "outlet"
            ]
        for (
            phase,
            obj,
        ) in self.precipitation_reactor.flow_mol_precipitate.items():
            self.rkt_block_outputs[("speciesAmount", phase)] = obj

        if self.config.add_alkalinity:
            self.precipitation_reactor.alkalinity = Var(
                initialize=1,
                units=pyunits.mg / pyunits.L,
                doc="Alkalinity (mg/L as CaCO3)",
            )
            self.rkt_block_outputs[("alkalinityAsCaCO3", None)] = (
                self.precipitation_reactor.alkalinity
            )

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
        if self.config.track_pE:
            self.reaktoro_options.system_state_option(
                "pE", self.precipitation_reactor.pE["inlet"]
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
        self.reaktoro_options["outputs"] = self.rkt_block_outputs
        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)

        self.precipitation_block = ReaktoroBlock(**self.reaktoro_options)

    def add_non_eq_reaktoro_chemistry(self):

        # self.precipitation_reactor.scaling_tendency = Var(
        #     ["non_eq_reactor", "maximum"],
        #     self.selected_precipitants,
        #     initialize=1,
        #     units=pyunits.dimensionless,
        # )
        self.precipitation_reactor.non_eq_flow_mol_precipitate = Var(
            self.selected_precipitants,
            initialize=1,
            bounds=(None, None),
            units=pyunits.mol / pyunits.s,
        )
        self.precipitation_reactor.non_eq_efficacy = Var(
            self.selected_precipitants,
            initialize=0.9,
            units=pyunits.dimensionless,
        )
        self.precipitation_reactor.eq_pH = Var(
            initialize=7,
            units=pyunits.dimensionless,
        )
        self.precipitation_reactor.non_eq_efficacy.fix()
        solvents = self.config.viable_reagents.create_solvent_constraint(
            self.precipitation_reactor, self.precipitation_reactor.flow_mol_reagent
        )
        reagents = {}
        reagents_and_solids = {}

        for r in self.selected_reagents:
            reagents[r] = self.precipitation_reactor.flow_mol_reagent[r]
            reagents_and_solids[r] = self.precipitation_reactor.flow_mol_reagent[r]

        for precip in self.selected_precipitants:
            reagents_and_solids[f"precipitation_{precip}"] = (
                self.precipitation_reactor.flow_mol_precipitate[precip]
            )
            self.precipitation_reactor.flow_mol_precipitate[precip].unfix()
        if solvents is not None:
            for solvent in solvents:
                reagents[solvent] = self.precipitation_reactor.flow_mol_solvent[solvent]
                reagents_and_solids[solvent] = (
                    self.precipitation_reactor.flow_mol_solvent[solvent]
                )

        self.non_eq_rkt_outputs = {
            (
                "pH",
                None,
            ): self.precipitation_reactor.pH["outlet"]
        }
        if self.config.track_pE:
            self.non_eq_rkt_outputs[("pE", None)] = self.precipitation_reactor.pE[
                "outlet"
            ]
        self.eq_rkt_outputs = {("pH", None): self.precipitation_reactor.eq_pH}
        if self.config.add_alkalinity:
            self.precipitation_reactor.alkalinity = Var(
                initialize=1,
                units=pyunits.mg / pyunits.L,
                doc="Alkalinity (mg/L as CaCO3)",
            )
            self.non_eq_rkt_outputs[("alkalinityAsCaCO3", None)] = (
                self.precipitation_reactor.alkalinity
            )
        chem_modifier_for_solids = (
            self.config.viable_reagents.get_reaktoro_chemistry_modifiers()
        )
        for (
            phase,
            obj,
        ) in self.precipitation_reactor.flow_mol_precipitate.items():
            self.eq_rkt_outputs[("speciesAmount", phase)] = (
                self.precipitation_reactor.non_eq_flow_mol_precipitate[phase]
            )
            # self.non_eq_rkt_outputs[("scalingTendency", phase)] = (
            #     self.precipitation_reactor.scaling_tendency["non_eq_reactor", phase]
            # )
            # self.precipitation_reactor.scaling_tendency["maximum", phase].fix(1)
            chem_modifier_for_solids[f"precipitation_{phase}"] = (
                self.selected_precipitants[phase]["reaktoro_modifier"]
            )
        chem_modifiers = self.config.viable_reagents.get_reaktoro_chemistry_modifiers()
        reaktoro_options = ReaktoroOptionsContainer()
        reaktoro_options.system_state_option(
            "temperature",
            self.precipitation_reactor.dissolution_reactor.properties_in[0].temperature,
        )
        reaktoro_options.system_state_option(
            "pressure",
            self.precipitation_reactor.dissolution_reactor.properties_in[0].pressure,
        )
        reaktoro_options.system_state_option(
            "pH", self.precipitation_reactor.pH["inlet"]
        )
        if self.config.track_pE:
            reaktoro_options.system_state_option(
                "pE", self.precipitation_reactor.pE["inlet"]
            )
        reaktoro_options["mineral_phase"] = {
            "phase_components": list(self.selected_precipitants.keys())
        }
        reaktoro_options.aqueous_phase_option(
            "composition",
            self.precipitation_reactor.dissolution_reactor.properties_in[
                0
            ].flow_mol_phase_comp,
        )
        reaktoro_options["register_new_chemistry_modifiers"] = chem_modifiers
        reaktoro_options["chemistry_modifier"] = reagents
        reaktoro_options["outputs"] = self.eq_rkt_outputs
        reaktoro_options.update_with_user_options(self.config.reaktoro_options)

        self.eq_precipitation_block = ReaktoroBlock(**reaktoro_options)

        non_eq_reaktoro_options = ReaktoroOptionsContainer()
        non_eq_reaktoro_options.system_state_option(
            "temperature",
            self.precipitation_reactor.dissolution_reactor.properties_in[0].temperature,
        )
        non_eq_reaktoro_options.system_state_option(
            "pressure",
            self.precipitation_reactor.dissolution_reactor.properties_in[0].pressure,
        )
        non_eq_reaktoro_options.system_state_option(
            "pH", self.precipitation_reactor.pH["inlet"]
        )
        if self.config.track_pE:
            reaktoro_options.system_state_option(
                "pE", self.precipitation_reactor.pE["inlet"]
            )
        non_eq_reaktoro_options.aqueous_phase_option(
            "composition",
            self.precipitation_reactor.dissolution_reactor.properties_in[
                0
            ].flow_mol_phase_comp,
        )
        non_eq_reaktoro_options["register_new_chemistry_modifiers"] = (
            chem_modifier_for_solids
        )
        non_eq_reaktoro_options["chemistry_modifier"] = reagents_and_solids
        non_eq_reaktoro_options["outputs"] = self.non_eq_rkt_outputs
        non_eq_reaktoro_options.update_with_user_options(self.config.reaktoro_options)

        self.non_eq_precipitation_block = ReaktoroBlock(**non_eq_reaktoro_options)

        @self.precipitation_reactor.Constraint(self.selected_precipitants)
        def precipitation_limited_reaction(fs, phase):
            return (
                fs.non_eq_flow_mol_precipitate[phase] * fs.non_eq_efficacy[phase]
                == fs.flow_mol_precipitate[phase]
            )

    def set_fixed_operation(self):
        for reagent, options in self.selected_reagents.items():
            self.precipitation_reactor.reagent_dose[reagent].fix(10 / 1000)
            if options["max_dose"] is not None:
                self.precipitation_reactor.reagent_dose[reagent].setub(
                    options["max_dose"] / 1000
                )
            if options["min_dose"] is not None:
                self.precipitation_reactor.reagent_dose[reagent].setlb(
                    options["min_dose"] / 1000
                )

            self.precipitation_reactor.flow_mol_reagent[reagent].setlb(0)

        self.precipitation_reactor.waste_mass_frac_precipitate.fix(0.2)
        for phase, ion in self.precipitation_reactor.separator.waste_state[
            0.0
        ].flow_mass_phase_comp:
            self.precipitation_reactor.separator.waste_state[0.0].flow_mass_phase_comp[
                phase, ion
            ].setlb(None)
            self.precipitation_reactor.separator.waste_state[0.0].flow_mol_phase_comp[
                phase, ion
            ].setlb(0)
        for precip in self.selected_precipitants.keys():
            self.precipitation_reactor.flow_mol_precipitate[precip].setlb(0)
            self.precipitation_reactor.flow_mass_precipitate[precip].setlb(None)
            if self.config.add_reaktoro_chemistry == False:
                self.precipitation_reactor.flow_mol_precipitate[precip].fix(1e-5)
            if self.config.add_non_eq_reaktoro_chemistry:
                self.precipitation_reactor.precipitation_limited_reaction.deactivate()

    def set_optimization_operation(self):
        """if we have reaktoro chemistry, we need to unfix the precipitate flow and reagent addition"""
        if (
            self.config.add_reaktoro_chemistry
            or self.config.add_non_eq_reaktoro_chemistry
        ):
            for phase, data in self.selected_precipitants.items():
                self.precipitation_reactor.flow_mol_precipitate[phase].unfix()
                self.precipitation_reactor.flow_mass_precipitate[phase].unfix()
            for reagent, options in self.selected_reagents.items():
                self.precipitation_reactor.reagent_dose[reagent].unfix()

    def scale_before_initialization(self, **kwargs):
        # TODO: Identify better auto scaling factors for dose/precipitation
        # generally it appears using the smallest does step is a good scaling factor
        # but it might change depending on chemicals being added.
        dose_scale = 1 / 0.001  # step size of ppm (or 0.001 kg/m3)
        reagent_ions = {}
        for reagent, options in self.selected_reagents.items():
            # use mol flow, as thats what will be propagated by default via mcas
            mass_flow_scale = dose_scale / (
                self.config.default_property_package._default_scaling_factors[
                    "flow_mol_phase_comp", ("Liq", "H2O")
                ]
                / value(self.config.default_property_package.mw_comp["H2O"])
            )
            iscale.set_scaling_factor(
                self.precipitation_reactor.flow_mass_reagent[reagent],
                mass_flow_scale,
            )
            mol_scale = mass_flow_scale * value(
                pyunits.convert(
                    self.config.viable_reagents[reagent]["mw"],
                    pyunits.kg / pyunits.mol,
                )
            )
            for ion, stoich in self.config.viable_reagents[reagent][
                "dissolution_stoichiometric"
            ].items():
                if ion in reagent_ions and reagent_ions[ion] < mol_scale / stoich:
                    # if we already have this ion, we want to use the largest scale
                    reagent_ions[ion] = mol_scale / stoich
                else:
                    reagent_ions[ion] = mol_scale / stoich
            iscale.set_scaling_factor(
                self.precipitation_reactor.flow_mol_reagent[reagent],
                mol_scale,
            )
            iscale.set_scaling_factor(
                self.precipitation_reactor.reagent_dose[reagent], dose_scale
            )

        for precip in self.selected_precipitants.keys():
            scales = []
            # use scale for ion and add scale from added chemicals

            for ion, stoich in self.selected_precipitants[precip][
                "precipitation_stoichiometric"
            ].items():
                sc = (
                    self.config.default_property_package._default_scaling_factors[
                        "flow_mol_phase_comp", ("Liq", ion)
                    ]
                    / stoich
                )
                if ion in reagent_ions:
                    # the dose scale is for ppm, assume we adding 1000 ppm on average
                    # this says that system has x amount of ion + 1000 ppm from reagent ion
                    # This could also be related to maximum chemical dose, but we assume its about 1000 ppm
                    sc = (sc**2 + (reagent_ions[ion] / 1000) ** 2) ** 0.5
                scales.append(sc)
            # want maximum scale for limiting ion
            precip_scale = max(scales)
            mol_precip_scale = precip_scale
            mass_precip_scale = precip_scale / value(
                value(
                    pyunits.convert(
                        self.config.viable_precipitants[precip]["mw"],
                        pyunits.kg / pyunits.mol,
                    )
                )
            )
            iscale.set_scaling_factor(
                self.precipitation_reactor.flow_mass_precipitate[precip],
                mass_precip_scale,
            )

            iscale.set_scaling_factor(
                self.precipitation_reactor.flow_mol_precipitate[precip],
                mol_precip_scale,
            )
            if self.config.add_non_eq_reaktoro_chemistry:
                iscale.set_scaling_factor(
                    self.precipitation_reactor.non_eq_flow_mol_precipitate[precip],
                    mol_precip_scale,
                )
                iscale.set_scaling_factor(
                    self.precipitation_reactor.non_eq_efficacy[precip], 1 / 10
                )
                iscale.constraint_scaling_transform(
                    self.precipitation_reactor.precipitation_limited_reaction[precip],
                    mol_precip_scale,
                )

        iscale.set_scaling_factor(self.precipitation_reactor.pH, 1 / 10)
        if self.config.track_pE:
            iscale.set_scaling_factor(self.precipitation_reactor.pE, 1 / 10)
        if self.config.add_non_eq_reaktoro_chemistry:
            iscale.set_scaling_factor(self.precipitation_reactor.eq_pH, 1 / 10)
            # for phase in self.selected_precipitants:
            #     iscale.set_scaling_factor(
            #         self.precipitation_reactor.scaling_tendency[
            #             "non_eq_reactor", phase
            #         ],
            #         1,
            #     )
            #     iscale.set_scaling_factor(
            #         self.precipitation_reactor.scaling_tendency["maximum", phase],
            #         1,
            #     )
        if (
            self.config.add_reaktoro_chemistry
            or self.config.add_non_eq_reaktoro_chemistry
        ):
            self.config.viable_reagents.scale_solvent_vars_and_constraints(
                self.precipitation_reactor
            )
            if self.config.add_alkalinity:
                iscale.set_scaling_factor(self.precipitation_reactor.alkalinity, 1)
        else:
            iscale.constraint_scaling_transform(self.eq_outlet_pH, 1 / 10)
            if self.config.track_pE:
                iscale.constraint_scaling_transform(self.eq_outlet_pE, 1 / 10)
        if self.config.add_hardness:
            iscale.set_scaling_factor(self.precipitation_reactor.hardness["Ca"], 1)
            iscale.set_scaling_factor(self.precipitation_reactor.hardness["Mg"], 1)
            iscale.set_scaling_factor(self.precipitation_reactor.hardness["Total"], 1)
            iscale.constraint_scaling_transform(
                self.precipitation_reactor.eq_hardness["Ca"], 1
            )
            iscale.constraint_scaling_transform(
                self.precipitation_reactor.eq_hardness["Mg"], 1
            )
            iscale.constraint_scaling_transform(
                self.precipitation_reactor.eq_hardness["Total"], 1
            )

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
            self.precipitation_block.initialize()
            self.precipitation_block.display_jacobian_scaling()
            self.precipitation_reactor.initialize()
            for phase, data in self.selected_precipitants.items():
                self.precipitation_reactor.flow_mol_precipitate[phase].unfix()
                self.precipitation_reactor.flow_mass_precipitate[phase].unfix()

        elif self.config.add_non_eq_reaktoro_chemistry:
            self.eq_precipitation_block.initialize()
            self.eq_precipitation_block.display_jacobian_scaling()

            for phase in self.selected_precipitants:

                self.precipitation_reactor.flow_mol_precipitate[phase].value = (
                    self.precipitation_reactor.non_eq_flow_mol_precipitate[phase].value
                    * self.precipitation_reactor.non_eq_efficacy[phase].value
                )
            self.non_eq_precipitation_block.initialize()
            self.non_eq_precipitation_block.display_jacobian_scaling()
            self.precipitation_reactor.initialize()
            for phase, data in self.selected_precipitants.items():
                self.precipitation_reactor.flow_mol_precipitate[phase].unfix()
                self.precipitation_reactor.flow_mass_precipitate[phase].unfix()
                self.precipitation_reactor.non_eq_flow_mol_precipitate[phase].unfix()
                self.precipitation_reactor.precipitation_limited_reaction.activate()

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
                return self.precipitation_reactor.pE[port]
            else:
                return None

        unit_dofs = degrees_of_freedom(self)

        treated_state = get_ion_comp(
            self.precipitation_reactor.precipitation_reactor.properties_out[0],
            self.precipitation_reactor.pH["outlet"],
            get_pe("outlet"),
        )

        if (
            self.config.add_alkalinity != False
            and self.precipitation_reactor.find_component("alkalinity") is not None
        ):
            treated_state["Alkalinity (mg/L as CaCO3)"] = (
                self.precipitation_reactor.alkalinity
            )
        if self.config.add_hardness:
            treated_state["Ca hardness (mg/L as CaCO3)"] = (
                self.precipitation_reactor.hardness["Ca"]
            )
            treated_state["Mg hardness (mg/L as CaCO3)"] = (
                self.precipitation_reactor.hardness["Mg"]
            )
            treated_state["Total hardness (mg/L as CaCO3)"] = (
                self.precipitation_reactor.hardness["Total"]
            )
        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "Inlet state": get_ion_comp(
                self.precipitation_reactor.dissolution_reactor.properties_in[0],
                self.precipitation_reactor.pH["inlet"],
            ),
            "Chemical dosing:": self.precipitation_reactor.reagent_dose,
            "Solids formed:": self.precipitation_reactor.flow_mass_precipitate,
            "Solids formed (mol basis):": self.precipitation_reactor.flow_mol_precipitate,
            "Sludge formed": self.precipitation_reactor.total_sludge_product,
            "Treated state": treated_state,
            "Waste state": get_ion_comp(
                self.precipitation_reactor.separator.waste_state[0],
                self.precipitation_reactor.pH["outlet"],
            ),
        }
        if hasattr(self, "rkt_block_outputs"):
            model_state_dict["Reaktoro outputs"] = self.rkt_block_outputs
        if hasattr(self, "eq_rkt_outputs"):
            model_state_dict["Equilibrium rkt outputs"] = self.eq_rkt_outputs
            model_state_dict["Non-equilibrium rkt outputs"] = self.non_eq_rkt_outputs
            model_state_dict["Reaktoro precipitation efficiency"] = (
                self.precipitation_reactor.non_eq_efficacy
            )
            # for key, value in self.precipitation_block.outputs.items():
            #         model_state_dict[key] = value
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

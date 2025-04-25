__author__ = "Alexander Dudchenko"
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.chemical_utils import (
    ViableReagentsBase,
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
        )
        self.register_reagent(
            "H2SO4",
            98.08 * pyunits.g / pyunits.mol,
            {"SO4_2-": 1, "H2O": 1},
            min_dose=0.1,
            max_dose=3000,
            purity=0.93,
            solvent=("H2O", 18.01 * pyunits.g / pyunits.mol),
        )


@declare_process_block_class("AcidificationUnit")
class AcidificationUnitUnitData(WaterTapFlowsheetBlockData):
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
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="To use Reaktoro-PSE for estimating amount of solids formed",
            doc="""
            If True, builds a reaktoro block and useses it to calculate how much solids form based on feed
            composition adn amount of added reagent. 
            """,
        ),
    )
    CONFIG.declare(
        "reaktoro_options",
        ConfigValue(
            default={
                "activity_model": "ActivityModelPitzer",
                "database": "PhreeqcDatabase",
                "database_file": "pitzer.dat",
                "reaktoro_block_manager": None,
            },
            description="Options for configuring Reaktoro-PSE",
            doc="""
            Options for configuring Reaktoro-PSE
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

        self.acidification_reactor = StoichiometricReactor(
            property_package=self.config.default_property_package,
            reagent=self.selected_reagents,
        )
        self.acidification_reactor.pH = Var(
            ["inlet", "outlet"],
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(1, 12),
        )
        if self.config.default_costing_package is not None:
            self.acidification_reactor.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package
            )

            for reagent, reagent_config in self.selected_reagents.items():
                self.acidification_reactor.add_component(
                    reagent,
                    Param(units=self.config.default_costing_package / pyunits.kg),
                    mutable=True,
                )
                self.acidification_reactor.find_component(reagent).set(
                    reagent_config["cost"]
                )
                self.acidification_reactor.costing.register_flow_type(
                    f"softening_reagent_{reagent}",
                    self.softening_reactor.find_component(reagent),
                )

        if self.config.add_reaktoro_chemistry:
            self.add_reaktoro_chemistry()

        self.register_port(
            "inlet",
            self.softening_reactor.inlet,
            {"pH": self.softening_reactor.pH["inlet"]},
        )
        self.register_port(
            "outlet",
            self.softening_reactor.outlet,
            {"pH": self.softening_reactor.pH["outlet"]},
        )

    def add_reaktoro_chemistry(self):
        reagents = {}
        for r in self.selected_reagents:
            reagents[r] = self.softening_reactor.flow_mol_reagent[r]

        outputs = {("pH", None): self.softening_reactor.pH["outlet"]}
        for phase, obj in self.softening_reactor.flow_mol_precipitate.items():
            outputs[("speciesAmount", phase)] = obj

        self.precipitation_block = ReaktoroBlock(
            system_state={
                "temperature": self.softening_reactor.dissolution_reactor.properties_in[
                    0
                ].temperature,
                "pressure": self.softening_reactor.dissolution_reactor.properties_in[
                    0
                ].pressure,
                "pH": self.softening_reactor.pH["inlet"],
            },
            aqueous_phase={
                "composition": self.softening_reactor.dissolution_reactor.properties_in[
                    0
                ].flow_mol_phase_comp,
                "activity_model": self.config.reaktoro_options["activity_model"],
                "fixed_solvent_specie": "H2O",
                "convert_to_rkt_species": True,
            },
            mineral_phase={"phase_components": list(self.selected_precipitants.keys())},
            chemistry_modifier=reagents,
            outputs=outputs,
            database=self.config.reaktoro_options["database"],
            database_file=self.config.reaktoro_options["database_file"],
            dissolve_species_in_reaktoro=True,
            build_speciation_block=True,
            assert_charge_neutrality=True,
            reaktoro_block_manager=self.config.reaktoro_options[
                "reaktoro_block_manager"
            ],
        )

    def fix_operation(self):
        for reagent, options in self.selected_reagents.items():
            self.softening_reactor.reagent_dose[reagent].setlb(
                options["min_dose"] / 1000
            )
            self.softening_reactor.reagent_dose[reagent].setub(
                options["max_dose"] / 1000
            )
            self.softening_reactor.flow_mol_reagent[reagent].setlb(None)

        self.softening_reactor.waste_mass_frac_precipitate.fix(0.2)

        for precip in self.selected_precipitants.keys():
            self.softening_reactor.flow_mol_precipitate[precip].setlb(None)

    def scale_before_initialization(self, **kwargs):
        max_dose = []
        for reagent, options in self.selected_reagents.items():
            dose_scale = options["max_dose"]
            max_dose.append(dose_scale)
            # use mol flow, as thats what will be propagated by default via mcas
            mass_flow_scale = dose_scale / value(
                self.softening_reactor.inlet.flow_mol_phase_comp[0.0, "Liq", "H2O"]
                * (18.015 / 1000)
            )
            iscale.set_scaling_factor(
                self.softening_reactor.flow_mass_reagent[reagent],
                mass_flow_scale,
            )
            iscale.set_scaling_factor(
                self.softening_reactor.reagent_dose[reagent], dose_scale
            )
        precip_scale = max(max_dose) / value(
            self.softening_reactor.inlet.flow_mol_phase_comp[0.0, "Liq", "H2O"]
            * (18.015 / 1000)
        )
        for precip in self.selected_precipitants.keys():
            iscale.set_scaling_factor(
                self.softening_reactor.flow_mass_precipitate[precip],
                precip_scale,
            )
        iscale.set_scaling_factor(self.softening_reactor.pH, 1)

    def initialize_unit(self, **kwargs):

        for phase, data in self.selected_precipitants.items():
            # assume that only fraction of ions will actual preciptaitate
            flow = (
                self.softening_reactor.inlet.flow_mol_phase_comp[
                    0.0, "Liq", data["primary_ion"]
                ].value
                * 0.0001
            )
            self.softening_reactor.flow_mol_precipitate[phase].fix(flow)
        for reagent, _ in self.selected_reagents.items():
            self.softening_reactor.reagent_dose[reagent].fix(10 / 1000)

        self.softening_reactor.initialize()
        if self.config.add_reaktoro_chemistry:
            # get intial mole flows
            self.precipitation_block.initialize()
            self.precipitation_block.display_jacobian_scaling()
            # recalcualte state with updated mol flow values
            self.softening_reactor.initialize()
            for phase, data in self.selected_precipitants.items():
                self.softening_reactor.flow_mol_precipitate[phase].unfix()

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
                self.softening_reactor.dissolution_reactor.properties_in[0],
                self.softening_reactor.pH["inlet"],
            ),
            "Chemical dosing:": self.softening_reactor.reagent_dose,
            "Solids formed:": self.softening_reactor.flow_mass_precipitate,
            "Treated state": get_ion_comp(
                self.softening_reactor.precipitation_reactor.properties_out[0],
                self.softening_reactor.pH["outlet"],
            ),
            "Waste state": get_ion_comp(
                self.softening_reactor.precipitation_reactor.properties_out[0],
                self.softening_reactor.pH["inlet"],
            ),
        }

        return self.name, model_state

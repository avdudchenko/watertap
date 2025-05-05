__author__ = "Alexander Dudchenko"
__author__ = "Alexander Dudchenko"

from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)

from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.reaktoro_utils import (
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
    units as pyunits,
)
from pyomo.common.config import ConfigValue

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock

import idaes.core.util.scaling as iscale

from idaes.models.unit_models import (
    Translator,
)
from reaktoro_pse.reaktoro_block import ReaktoroBlock


@declare_process_block_class("MultiCompPumpUnit")
class MultiCompPumpUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "initialization_pressure",
        ConfigValue(
            default="osmotic_pressure",
            description="Pressure to use for initial guess",
            doc="""
            Can be:
                - a value with pyomo units (e.g. 1*pyunits.bar)
                - 'osmotic_pressure' - will use osmotic pressure entering pump to initialize its outlet with
            """,
        ),
    )
    CONFIG.declare(
        "ro_property_package",
        ConfigValue(
            default=None,
            description="Property package to use in RO model",
            doc="""
            This defines which property package should be used in RO model, either NaCl_prop_pack or seawater_prop_pack,
            if non is provided, will use seawater_prop_pack
            """,
        ),
    )
    CONFIG.declare(
        "selected_scalants",
        ConfigValue(
            default={"Calcite": 1, "Gypsum": 1},
            description="Dict of scalants to track",
            doc="""
            Defines which scalants to track, and their maximum potential during optimization. 
            Provide a dictionary with following structure:
            {'Scalant': maximum scaling potential}
            if add_reaktoro_chemistry == True:
                Will add maximum_scaling_potential variable, and use int eq_maximum_scaling_potential constraint to define scaling limits
                The eq_maximum_scaling_potential is deactivated until set_optimization_operation is called.
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
    CONFIG.declare(
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="To use Reaktoro-PSE for estimating scaling limits in RO",
            doc="""
            If True, builds a reaktoro block and uses it to calculate scaling potential, will also create scaling constraints. 
            """,
        ),
    )
    CONFIG.declare(
        "ro_options_dict",
        default=None,
        description="Options for RO, will override the defaults",
        doc="""
            Provide dict with options to change defaults in RO model,
            {'has_pressure_change:True} etc. 
            This will update default dictionary. 
            """,
    )
    CONFIG.declare(
        "build_monotonic_cp_constraint",
        default=True,
        description="Defines if monotonic concentration polarization constraint is added",
        doc="""
                   Builds a monotonic concentration polarization constraint to ensure CP is always highest at the end of the
                   module, this alow construction of a single Reaktoro block for monitoring scaling, if its not built during
                   optimization model might design system to operate with maximum CP in middle of the module or way from where
                   Reaktoro gets its composition information""",
    )

    def build(self):
        super().build()
        self.get_ro_solute_type()
        # define translator blocks
        self.ro_feed = Translator(
            inlet_property_package=self.config.default_config_block,
            outlet_property_package=self.config.ro_property_package,
        )
        self.ro_retentate = Translator(
            inlet_property_package=self.config.ro_property_package,
            outlet_property_package=self.config.default_config_block,
        )
        self.ro_product = Translator(
            inlet_property_package=self.config.ro_property_package,
            outlet_property_package=self.config.default_config_block,
        )
        # set them up for translating input prop pack to outlet prop pack
        self.setup_inlet_translator_block(self.ro_feed)

        # TODO: define outlet blocks to gether, so we can include pseudo rejection of ions
        self.setup_outlet_translator_block(self.ro_retentate)
        self.setup_outlet_translator_block(self.ro_product)

        # build ro unit, we will grab ro options, and redfine them with user provided overrides
        self.ro_unit = ReverseOsmosis1D(**self.get_ro_options())

        if self.config.default_costing_package is not None:
            self.ro_unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package
            )
        self.ro_feed.pH = Var(initialize=self.config.pH, units=pyunits.dimensionless)
        self.ro_retentate.pH = Var(
            initialize=self.config.pH, units=pyunits.dimensionless
        )
        self.ro_product.pH = Var(initialize=self.config.pH, units=pyunits.dimensionless)

        self.register_port("feed", self.ro_feed.inlet, self.ro_feed.pH)
        self.register_port("feed", self.ro_retentate.inlet, self.ro_retentate.pH)
        self.register_port("feed", self.ro_product.inlet, self.ro_product.pH)

        if self.config.build_monotone_cp_constraint:
            self.build_monotonic_cp_constraint()

        if self.add_reaktoro_chemistry:
            self.build_scaling_constraints()
            self.add_reaktoro_chemistry()

    def get_ro_options(self):
        """defines ro defaults and overrides them with user config options if provided"""
        default_ro_dict = {
            "property_package": self.config.ro_property_package,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 10,
        }
        if self.config.ro_options_dict is not None:
            default_ro_dict.update(self.config.ro_options_dict)
        return default_ro_dict

    def setup_inlet_translator_block(
        self,
        translator_block,
    ):
        """defines inlet translator block, will sum up all mcas massflow and translate them to mass of single solute
        in ro property package
        Args:
            translator_block -- Inlet translator block (should be feed)"""
        ions = []
        for index in translator_block.properties_in[0].flow_mol_phase_comp:
            if "H2O" not in index:
                ions.append(
                    translator_block.properties_in[0].flow_mol_phase_comp[index]
                    * translator_block.properties_in[0].mw_comp[index[1]]
                )

        @translator_block.Constraint(["H2O", self.prop_pack_ion])
        def eq_flow_mass_phase_comp(blk, ion):
            if ion == "H2O":
                return (
                    translator_block.properties_in[0].flow_mol_phase_comp["Liq", "H2O"]
                    * translator_block.properties_in[0].mw_comp["H2O"]
                    == translator_block.properties_out[0].flow_mass_phase_comp[
                        "Liq", "H2O"
                    ]
                )
            else:
                return translator_block.properties_out[0].flow_mass_phase_comp[
                    "Liq", self.ro_solute_type
                ] == sum(ions)

        translator_block.eq_pressure_equality = Constraint(
            expr=translator_block.properties_in[0].pressure
            == translator_block.properties_out[0].pressure
        )
        iscale.constraint_scaling_transform(translator_block.pressures_constraint, 1e-5)
        translator_block.eq_temperature_equality = Constraint(
            expr=translator_block.properties_in[0].temperature
            == translator_block.properties_out[0].temperature
        )
        iscale.constraint_scaling_transform(translator_block.temp_constraint, 1e-2)
        translator_block.properties_out[0].pressure_osm_phase[...]

    def setup_outlet_translator_block(self, translator_block, inlet_composition):
        """defines outlet translator block, will convert single outlet solute to multi solutes, assuming they are
        ratiometricly related to changes between inlet and outlet properties. e.g.
         out_ion=in_total_ion_mass/out_total_ion_mass*in_ion_mass
        Args:
            translator_block -- Outlet translator block (should be feed)
            inlet_composition -- Inlet composition used to get original mass flow of ions entering system
        """
        tds_in = []
        for index in inlet_composition:
            if "H2O" not in index:
                tds_in.append(
                    inlet_composition[index]
                    * translator_block.properties_out[0].mw_comp[index[-1]]
                )

        @translator_block.Constraint(
            list(translator_block.properties_out[0].flow_mass_phase_comp)
        )
        def eq_flow_mass_phase_comp(blk, liq, ion):
            if "H2O" in ion:
                return (
                    blk.properties_out[0].flow_mol_phase_comp[liq, ion]
                    * blk.properties_out[0].mw_comp[ion]
                    == blk.properties_in[0].flow_mass_phase_comp[liq, ion]
                )
            else:
                return (
                    blk.properties_out[0].flow_mol_phase_comp[liq, ion] * sum(tds_in)
                    == inlet_composition[liq, ion]
                    * blk.properties_in[0].flow_mass_phase_comp[
                        "Liq", self.ro_solute_type
                    ]
                )

        translator_block.eq_pressure_equality = Constraint(
            expr=translator_block.properties_in[0].pressure
            == translator_block.properties_out[0].pressure
        )

        iscale.constraint_scaling_transform(translator_block.pressures_constraint, 1e-5)
        translator_block.eq_temperature_equality = Constraint(
            expr=translator_block.properties_in[0].temperature
            == translator_block.properties_out[0].temperature
        )

        iscale.constraint_scaling_transform(translator_block.temp_constraint, 1e-2)

    def get_ro_solute_type(self):
        if len(self.config.ro_property_package.solute_set) > 1:
            raise TypeError(
                "current multi_comp_ro model expects a single solute RO property package"
            )
        self.ro_solute_type = self.config.ro_property_package.solute_set[0]

    def build_monotonic_cp_constraint(self):
        """builds monotone concentration polarization constraint"""
        domain = list(self.ro_unit.length_domain)[1:]
        length_index = list(range(len(domain) - 1))

        @self.ro_unit.Constraint(domain)
        def monotone_cp_constraint(fs, d):
            return (
                self.ro_unit.feed_side.properties_interface[
                    0.0, domain[d]
                ].conc_mass_phase_comp["Liq", self.ro_solute_type]
                <= self.ro_unit.feed_side.properties_interface[
                    0.0, domain[d + 1]
                ].conc_mass_phase_comp["Liq", self.ro_solute_type]
            )

        for eq in self.unit_block.monotone_cp_constraint:
            iscale.constraint_scaling_transform(
                self.ro_unit.monotone_cp_constraint[eq], 1 / 100
            )

    def build_scaling_constraints(self):
        """builds scaling constraints"""
        self.ro_unit.scaling_tendency = Var(
            list(self.config.selected_scalants.keys()),
            initialize=1,
            units=pyunits.dimensionless,
        )

        self.ro_unit.maximum_scaling_tendency = Var(
            list(self.config.selected_scalants.keys()),
            initialize=lambda scalant: self.config.selected_scalants[scalant],
            units=pyunits.dimensionless,
        )
        self.ro_unit.maximum_scaling_tendencies.fix()

        @self.ro_unit.Constraint
        def eq_max_scaling_tendency(blk, scalant):
            return blk.scaling_tendency <= blk.maximum_scaling_tendency

        for scalant, max_tendency in self.config.selected_scalants:
            iscale.set_scaling_factor(
                self.ro_unit.scaling_tendency[scalant], 1 / max_tendency
            )
            iscale.set_scaling_factor(
                self.ro_unit.maximum_scaling_tendencies[scalant], 1 / max_tendency
            )
            iscale.constraint_scaling_transform(
                self.ro_unit.eq_max_scaling_tendency[scalant], 1 / max_tendency
            )

    def add_reaktoro_chemistry(self):
        """add water removal constraint, and relevant reaktoro block"""
        self.ro_unit.water_removed_at_interface = Var(
            initialize=1, units=pyunits.mol / pyunits.s
        )

        ro_cp_interface = self.ro_unit.feed_side.properties_interface[0, 1]
        self.ro_unit.eq_water_removed_at_interface = Constraint(
            expr=(
                self.ro_unit.water_removed_at_interface
                * self.config.default_property_package.mw_comp["H2O"]
            )
            == self.ro_unit.inlet.flow_mass_phase_comp["Liq", "H2O"]
            - self.ro_unit.inlet.flow_mass_phase_comp["Liq", self.ro_solute_type]
            / ro_cp_interface.mass_frac_phase_comp["Liq", self.ro_solute_type]
        )

        outputs = {("pH", None): self.ro_retentate.pH}
        for scalant in self.ro_unit.scaling_tendency:
            outputs[("scalingTendency", scalant)] = self.ro_unit.scaling_tendency[
                scalant
            ]

        self.reaktoro_options = ReaktoroOptionsContainer()
        self.reaktoro_options.system_state_option(
            "temperature",
            self.ro_feed.dissolution_reactor.properties_in[0]
            .properties_in[0]
            .temperature,
        )
        self.reaktoro_options.system_state_option(
            "pressure",
            self.ro_feed.dissolution_reactor.properties_in[0].properties_in[0].pressure,
        )
        self.reaktoro_options.system_state_option("pH", self.ro_feed.pH)
        self.reaktoro_options.aqueous_phase_option(
            "composition",
            self.ro_feed.dissolution_reactor.properties_in[0].flow_mol_phase_comp,
        )
        self.reaktoro_options["register_new_chemistry_modifiers"] = (
            self.config.viable_reagents.get_reaktoro_chemistry_modifiers()
        )
        self.reaktoro_options["chemistry_modifier"] = {
            "H2O_evaporation": self.ro_unit.water_removed_at_interface
        }
        self.reaktoro_options.system_state_modifier_option(
            "pressure", ro_cp_interface.pressure
        )

        self.reaktoro_options["outputs"] = outputs

        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)

        self.precipitation_block = ReaktoroBlock(**self.reaktoro_options)

    def set_fixed_operation(self):
        pass

    def set_optimization_operation(self):
        pass

    def scale_before_initialization(self, **kwargs):
        iscale.set_scaling_factor(self.pump.outlet.pressure, 1e-5)
        iscale.set_scaling_factor(self.pump.control_volume.work, 1e-4)

    def initialize_unit(self):
        self.set_fixed_operation()
        self.pump.initialize()

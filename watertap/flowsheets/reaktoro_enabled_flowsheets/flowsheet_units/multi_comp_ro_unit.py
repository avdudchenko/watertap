__author__ = "Alexander Dudchenko"
__author__ = "Alexander Dudchenko"

from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from idaes.core.util.initialization import propagate_state
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
    value,
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
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from idaes.models.unit_models import (
    Translator,
)
from reaktoro_pse.reaktoro_block import ReaktoroBlock
from pyomo.network import Arc


@declare_process_block_class("MultiCompROUnit")
class MultiCompROUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
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
        "add_reaktoro_chemistry",
        ConfigValue(
            default=True,
            description="To use Reaktoro-PSE for estimating scaling potential",
            doc="""
            If True, builds a reaktoro block and uses it to calculate scaling potential at membrane water interface in
            final RO node.
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
        "ro_options_dict",
        ConfigValue(
            default=None,
            description="Options for RO, will override the defaults",
            doc="""
            Provide dict with options to change defaults in RO model,
            {'has_pressure_change:True} etc. 
            This will update default dictionary. 
            """,
        ),
    )
    CONFIG.declare(
        "build_monotonic_cp_constraint",
        ConfigValue(
            default=True,
            description="Defines if monotonic concentration polarization constraint is added",
            doc="""
                   Builds a monotonic concentration polarization constraint to ensure CP is always highest at the end of the
                   module, this alow construction of a single Reaktoro block for monitoring scaling, if its not built during
                   optimization model might design system to operate with maximum CP in middle of the module or way from where
                   Reaktoro gets its composition information""",
        ),
    )

    def build(self):
        super().build()
        self.get_ro_solute_type()
        # define translator blocks
        self.ro_feed = Translator(
            inlet_property_package=self.config.default_property_package,
            outlet_property_package=self.config.ro_property_package,
        )
        self.ro_feed.properties_out[0].flow_mol_phase_comp[...]
        self.ro_feed.properties_out[0].conc_mass_phase_comp[...]
        self.ro_retentate = Translator(
            inlet_property_package=self.config.ro_property_package,
            outlet_property_package=self.config.default_property_package,
        )
        self.ro_retentate.properties_in[0].flow_mol_phase_comp[...]
        self.ro_retentate.properties_in[0].conc_mass_phase_comp[...]
        self.ro_product = Translator(
            inlet_property_package=self.config.ro_property_package,
            outlet_property_package=self.config.default_property_package,
        )
        self.ro_product.properties_in[0].flow_mol_phase_comp[...]
        self.ro_product.properties_in[0].conc_mass_phase_comp[...]
        # set them up for translating input prop pack to outlet prop pack
        self.setup_inlet_translator_block(self.ro_feed)

        # TODO: define outlet blocks to gether, so we can include pseudo rejection of ions
        self.setup_outlet_translator_block(
            self.ro_retentate, self.ro_feed.inlet.flow_mol_phase_comp
        )
        self.setup_outlet_translator_block(
            self.ro_product, self.ro_feed.inlet.flow_mol_phase_comp
        )

        # build ro unit, we will grab ro options, and redfine them with user provided overrides
        self.ro_unit = ReverseOsmosis1D(**self.get_ro_options())

        if self.config.default_costing_package is not None:
            self.ro_unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package,
                **self.config.default_costing_package_kwargs
            )
        self.ro_feed.pH = Var(initialize=7, bounds=(1, 12), units=pyunits.dimensionless)
        self.ro_retentate.pH = Var(
            initialize=7, bounds=(1, 12), units=pyunits.dimensionless
        )
        self.ro_product.pH = Var(
            initialize=7, bounds=(1, 12), units=pyunits.dimensionless
        )

        self.register_port("feed", self.ro_feed.inlet, {"pH": self.ro_feed.pH})
        self.register_port(
            "retentate", self.ro_retentate.outlet, {"pH": self.ro_retentate.pH}
        )
        self.register_port(
            "product", self.ro_product.outlet, {"pH": self.ro_product.pH}
        )

        if self.config.build_monotonic_cp_constraint:
            self.build_monotonic_cp_constraint()

        self.build_water_removal_constraint()

        if self.config.add_reaktoro_chemistry:
            self.build_scaling_constraints()
            self.add_reaktoro_chemistry()
        else:
            self.add_retentate_ph_constraint()
        self.add_permeate_ph_constraint()

        # connecting translator blocks to ro unit
        self.feed_to_ro = Arc(
            source=self.ro_feed.outlet,
            destination=self.ro_unit.inlet,
        )
        self.ro_to_retentate = Arc(
            source=self.ro_unit.retentate,
            destination=self.ro_retentate.inlet,
        )
        self.ro_to_product = Arc(
            source=self.ro_unit.permeate,
            destination=self.ro_product.inlet,
        )

    def add_retentate_ph_constraint(self):
        """adds retentate pH constraint"""
        self.ro_retentate.eq_ph_equality = Constraint(
            expr=self.ro_retentate.pH == self.ro_feed.pH
        )

    def add_permeate_ph_constraint(self):
        """adds permeate pH constraint, which we assume is average of feed and retentate"""
        self.ro_product.eq_average_permeate_ph = Constraint(
            expr=self.ro_product.pH
            == 0.5 * self.ro_feed.pH + 0.5 * self.ro_retentate.pH
        )

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

        @translator_block.Constraint(["H2O", self.ro_solute_type])
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
        translator_block.eq_temperature_equality = Constraint(
            expr=translator_block.properties_in[0].temperature
            == translator_block.properties_out[0].temperature
        )

    def setup_outlet_translator_block(self, translator_block, inlet_composition):
        """defines outlet translator block, will convert single outlet solute to multi solutes, assuming they are
        ratiometrically related to changes between inlet and outlet properties. e.g.
         out_ion=in_total_ion_mass/out_total_ion_mass*in_ion_mass
        Args:
            translator_block -- Outlet translator block (should be feed)
            inlet_composition -- Inlet composition used to get original mass flow of ions entering system
        """
        tds_in = []
        ions = []
        for index in inlet_composition:
            if "H2O" not in index:
                tds_in.append(
                    inlet_composition[index]
                    * translator_block.properties_out[0].mw_comp[index[-1]]
                )
            ions.append(index[-1])

        @translator_block.Constraint(ions)
        def eq_flow_mass_phase_comp(blk, ion):
            liq = "Liq"
            if "H2O" in ion:
                return (
                    blk.properties_out[0].flow_mol_phase_comp[liq, ion]
                    * blk.properties_out[0].mw_comp[ion]
                    == blk.properties_in[0].flow_mass_phase_comp[liq, ion]
                )
            else:
                return (
                    blk.properties_out[0].flow_mol_phase_comp[liq, ion] * sum(tds_in)
                    == inlet_composition[0.0, liq, ion]
                    * blk.properties_in[0].flow_mass_phase_comp[
                        "Liq", self.ro_solute_type
                    ]
                )

        translator_block.eq_pressure_equality = Constraint(
            expr=translator_block.properties_in[0].pressure
            == translator_block.properties_out[0].pressure
        )

        translator_block.eq_temperature_equality = Constraint(
            expr=translator_block.properties_in[0].temperature
            == translator_block.properties_out[0].temperature
        )

    def get_ro_solute_type(self):
        if len(self.config.ro_property_package.solute_set) > 1:
            raise TypeError(
                "current multi_comp_ro model expects a single solute RO property package"
            )
        self.ro_solute_type = list(self.config.ro_property_package.solute_set)[0]

    def build_monotonic_cp_constraint(self):
        """builds monotone concentration polarization constraint"""
        domain = list(self.ro_unit.length_domain)[1:]
        length_index = list(range(len(domain) - 1))

        @self.ro_unit.Constraint(length_index)
        def monotone_cp_constraint(fs, d):
            return (
                self.ro_unit.feed_side.properties_interface[
                    0.0, domain[d]
                ].conc_mass_phase_comp["Liq", self.ro_solute_type]
                <= self.ro_unit.feed_side.properties_interface[
                    0.0, domain[d + 1]
                ].conc_mass_phase_comp["Liq", self.ro_solute_type]
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
            initialize=lambda blk, idx: self.config.selected_scalants[idx],
            units=pyunits.dimensionless,
        )
        self.ro_unit.maximum_scaling_tendency.fix()

        @self.ro_unit.Constraint(list(self.config.selected_scalants.keys()))
        def eq_max_scaling_tendency(blk, scalant):
            return (
                blk.scaling_tendency[scalant] <= blk.maximum_scaling_tendency[scalant]
            )

        self.deactivate_scaling_constraints()

    def activate_scaling_constraints(self):
        """activates scaling constraints"""
        if self.config.add_reaktoro_chemistry:
            for scalant in self.config.selected_scalants.keys():
                self.ro_unit.eq_max_scaling_tendency[scalant].activate()
                self.ro_unit.scaling_tendency[scalant].unfix()
                self.ro_unit.maximum_scaling_tendency[scalant].fix()

    def deactivate_scaling_constraints(self):
        """deactivates scaling constraints"""
        if self.config.add_reaktoro_chemistry:
            for scalant in self.config.selected_scalants.keys():
                self.ro_unit.eq_max_scaling_tendency[scalant].deactivate()
                self.ro_unit.maximum_scaling_tendency[scalant].fix()

    def build_water_removal_constraint(self):
        """builds water removal constraint"""
        self.ro_unit.water_removed_at_interface = Var(
            initialize=1, units=pyunits.mol / pyunits.s
        )

        ro_cp_interface = self.ro_unit.feed_side.properties_interface[0, 1]
        self.ro_unit.eq_water_removed_at_interface = Constraint(
            expr=(
                self.ro_unit.water_removed_at_interface
                * self.config.default_property_package.mw_comp["H2O"]
            )
            == self.ro_unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
            - self.ro_unit.inlet.flow_mass_phase_comp[0, "Liq", self.ro_solute_type]
            / ro_cp_interface.mass_frac_phase_comp["Liq", self.ro_solute_type]
        )

    def add_reaktoro_chemistry(self):
        """add water removal constraint, and relevant reaktoro block"""

        outputs = {("pHDirect", None): self.ro_retentate.pH}
        for scalant in self.ro_unit.scaling_tendency:
            outputs[("scalingTendencyDirect", scalant)] = self.ro_unit.scaling_tendency[
                scalant
            ]

        ro_cp_interface = self.ro_unit.feed_side.properties_interface[0, 1]

        self.reaktoro_options = ReaktoroOptionsContainer()
        self.reaktoro_options.system_state_option(
            "temperature",
            self.ro_feed.properties_in[0].temperature,
        )
        self.reaktoro_options.system_state_option(
            "pressure",
            self.ro_feed.properties_in[0].pressure,
        )
        self.reaktoro_options.system_state_option("pH", self.ro_feed.pH)
        self.reaktoro_options.aqueous_phase_option(
            "composition",
            self.ro_feed.properties_in[0].flow_mol_phase_comp,
        )

        self.reaktoro_options["chemistry_modifier"] = {
            "H2O_evaporation": self.ro_unit.water_removed_at_interface
        }
        self.reaktoro_options.system_state_modifier_option(
            "pressure", ro_cp_interface.pressure
        )

        self.reaktoro_options["outputs"] = outputs
        self.reaktoro_options.update_with_user_options(self.config.reaktoro_options)
        self.scaling_block = ReaktoroBlock(**self.reaktoro_options)

    def set_fixed_operation(self):
        """fixes operation point for pump unit model"""
        self.ro_unit = self.ro_unit
        # Default membrane performance metrics
        self.ro_unit.A_comp[0, "H2O"].fix(2 / (3600 * 1000 * 1e5))
        self.ro_unit.B_comp[0, self.ro_solute_type].fix(0.15 / (3600 * 1000))

        # Default module design metrics
        self.ro_unit.area.fix(100)
        self.ro_unit.length.unfix()
        self.ro_unit.width.unfix()
        # guess for initial operation, this gets length and width
        self.ro_unit.feed_side.velocity[0, 0].fix(0.1)

        # module design
        self.ro_unit.feed_side.channel_height.fix(1 / 10 / 100)
        self.ro_unit.feed_side.spacer_porosity.fix(0.9)
        self.ro_unit.permeate.pressure[0].fix(101325)

        self.ro_unit.recovery_vol_phase[0.0, "Liq"].setlb(0.05)
        self.ro_unit.recovery_vol_phase[0.0, "Liq"].setub(0.98)
        self.ro_unit.feed_side.velocity[0, 0].fix(0.1)
        for d in list(self.ro_unit.length_domain):
            self.ro_unit.feed_side.velocity[0, d].setub(0.3)
            self.ro_unit.feed_side.velocity[0, d].setlb(0.0001)

        # Deal with low flows
        self.ro_unit.feed_side.cp_modulus.setub(50)
        self.ro_unit.feed_side.cp_modulus.setlb(0.1)
        # Deals with low feed salinity
        for e in self.ro_unit.permeate_side:
            if e[-1] != 0:
                self.ro_unit.permeate_side[e].pressure_osm_phase["Liq"].setlb(200)
                self.ro_unit.permeate_side[e].molality_phase_comp[
                    "Liq", self.ro_solute_type
                ].setlb(1e-8)

        self.ro_unit.feed_side.K.setlb(1e-6)
        self.ro_unit.feed_side.friction_factor_darcy.setub(200)
        self.ro_unit.flux_mass_phase_comp.setub(1)
        print("DOFS", degrees_of_freedom(self))
        self.report()

    def set_optimization_operation(self):
        """unfixes operation for optimization"""
        self.ro_unit.area.unfix()
        self.ro_unit.width.unfix()
        self.ro_unit.length.unfix()
        # self.ro_unit.area.setub(1000)
        # self.ro_unit.area.setub(None)
        self.ro_unit.width.setub(None)
        self.ro_unit.length.setub(None)
        self.ro_unit.area.setlb(1)
        self.ro_unit.width.setlb(0.1)
        self.ro_unit.length.setlb(0.1)
        self.ro_unit.feed_side.velocity[0, 0].unfix()
        self.ro_unit.recovery_vol_phase[0.0, "Liq"].unfix()
        self.ro_unit.recovery_vol_phase[0.0, "Liq"].setlb(0.2)
        self.ro_unit.recovery_vol_phase[0.0, "Liq"].setub(0.95)
        self.activate_scaling_constraints()

    def scale_before_initialization(self, **kwargs):
        """scales the unit before initialization"""
        iscale.set_scaling_factor(self.ro_feed.pH, 1 / 1000)
        iscale.set_scaling_factor(self.ro_retentate.pH, 1 / 1000)
        iscale.set_scaling_factor(self.ro_product.pH, 1 / 1000)
        # scale monotone cp constraint if it exists
        if self.ro_unit.find_component("monotone_cp_constraint") is not None:
            for eq in self.ro_unit.monotone_cp_constraint:
                iscale.constraint_scaling_transform(
                    self.ro_unit.monotone_cp_constraint[eq], 1
                )

        # Get the scaling factors for the RO property package and default property package
        prop_scaling = self.get_default_scaling_factors()
        # Configure default scaling factors for the RO property package
        self.config.ro_property_package.set_default_scaling(
            "flow_mass_phase_comp",
            prop_scaling["H2O"],
            index=("Liq", "H2O"),
        )
        self.config.ro_property_package.set_default_scaling(
            "flow_mass_phase_comp",
            prop_scaling[self.ro_solute_type],
            index=("Liq", self.ro_solute_type),
        )
        # Scale the translator blocks
        for block in [self.ro_feed, self.ro_retentate, self.ro_product]:
            iscale.constraint_scaling_transform(block.eq_pressure_equality, 1e-5)
            iscale.constraint_scaling_transform(block.eq_temperature_equality, 1e-2)
            iscale.constraint_scaling_transform(
                block.eq_flow_mass_phase_comp["H2O"], prop_scaling["H2O"]
            )
            # Scale the mass flow of the solute in the translator block,
            # for inlet, we use scale from ro_prop_scaling, for outlet we use default scaling
            # for default property package, Water scaling is same for both

            for ion in block.eq_flow_mass_phase_comp:
                sf = prop_scaling[ion]
                iscale.constraint_scaling_transform(
                    block.eq_flow_mass_phase_comp[ion],
                    sf,
                )

        # scale water removal constraint if it exists
        if self.ro_unit.find_component("eq_water_removed_at_interface") is not None:
            iscale.constraint_scaling_transform(
                self.ro_unit.eq_water_removed_at_interface, prop_scaling["H2O"]
            )
            iscale.set_scaling_factor(
                self.ro_unit.water_removed_at_interface, prop_scaling["H2O"]
            )
        # scale ph constraints
        if self.ro_retentate.find_component("eq_ph_equality") is not None:
            iscale.constraint_scaling_transform(self.ro_retentate.eq_ph_equality, 1)
        iscale.constraint_scaling_transform(self.ro_product.eq_average_permeate_ph, 1)
        # scale scaling constraints if they exist
        if self.config.add_reaktoro_chemistry:
            for scalant, max_tendency in self.config.selected_scalants.items():
                iscale.set_scaling_factor(
                    self.ro_unit.scaling_tendency[scalant], 1 / max_tendency
                )
                iscale.set_scaling_factor(
                    self.ro_unit.maximum_scaling_tendency[scalant], 1 / max_tendency
                )
                iscale.constraint_scaling_transform(
                    self.ro_unit.eq_max_scaling_tendency[scalant], 1 / max_tendency
                )

        # scale RO unit
        iscale.set_scaling_factor(self.ro_unit.area, 1 / self.ro_unit.area.value)
        iscale.constraint_scaling_transform(
            self.ro_unit.eq_area, 1 / self.ro_unit.area.value
        )
        iscale.set_scaling_factor(self.ro_unit.width, 1 / self.ro_unit.width.value)
        iscale.set_scaling_factor(self.ro_unit.length, 1 / self.ro_unit.length.value)

    def get_default_scaling_factors(self):
        """returns scale for ro property package and default property package"""
        scale_factors = {}

        scale_factors["H2O"] = (
            self.config.default_property_package._default_scaling_factors[
                "flow_mol_phase_comp", ("Liq", "H2O")
            ]
            / value(self.config.default_property_package.mw_comp["H2O"])
        )
        tds_sf = []
        for ion in self.config.default_property_package.solute_set:
            sf = self.config.default_property_package._default_scaling_factors[
                "flow_mol_phase_comp", ("Liq", ion)
            ]

            scale_factors[ion] = sf
            tds_sf.append(sf * value(self.config.default_property_package.mw_comp[ion]))

        tds_sf = sum(tds_sf)
        scale_factors[self.ro_solute_type] = tds_sf
        return scale_factors

    def init_translator_block(self, block, additonal_var_to_fix=None):
        """initializes translator block"""
        block.initialize()
        if additonal_var_to_fix is not None:
            if not isinstance(additonal_var_to_fix, list):
                additonal_var_to_fix = [additonal_var_to_fix]
            for var in additonal_var_to_fix:
                var.fix()
        flags = fix_state_vars(block.properties_in)
        solver = get_solver()
        results = solver.solve(block, tee=False)
        assert_optimal_termination(results)
        revert_state_vars(block.properties_in, flags)
        if additonal_var_to_fix is not None:
            for var in additonal_var_to_fix:
                var.unfix()

    def initialize_unit(self):
        self.init_translator_block(self.ro_feed)
        propagate_state(self.feed_to_ro)
        self.ro_unit.initialize()
        propagate_state(self.ro_to_product)
        propagate_state(self.ro_to_retentate)
        self.ro_feed.properties_in[
            0
        ].flow_mol_phase_comp.fix()  # need to fix for initialization
        self.init_translator_block(
            self.ro_retentate, self.ro_feed.properties_in[0].flow_mol_phase_comp
        )
        self.init_translator_block(
            self.ro_product, self.ro_feed.properties_in[0].flow_mol_phase_comp
        )
        self.ro_feed.properties_in[0].flow_mol_phase_comp.unfix()
        if self.config.add_reaktoro_chemistry:
            self.scaling_block.initialize()

    def get_model_state_dict(self):
        """Returns a dictionary with the model state"""
        unit_dofs = degrees_of_freedom(self)
        ro_domains = list(self.ro_unit.length_domain)
        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "RO design": {
                "Area": self.ro_unit.area,
                "Length": self.ro_unit.length,
                "Width": self.ro_unit.width,
            },
            "RO operation": {
                "Average flux": pyunits.convert(
                    self.ro_unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"],
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.hr),
                ),
                "Inlet flux": pyunits.convert(
                    self.ro_unit.flux_mass_phase_comp[0, ro_domains[1], "Liq", "H2O"],
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.hr),
                ),
                "Outelt flux": pyunits.convert(
                    self.ro_unit.flux_mass_phase_comp[0, ro_domains[-1], "Liq", "H2O"],
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.hr),
                ),
                "Inlet velocity": self.ro_unit.feed_side.velocity[0, 0],
                "Recovery": self.ro_unit.recovery_vol_phase[0.0, "Liq"],
                "Rejection": self.ro_unit.rejection_phase_comp[
                    0, "Liq", self.ro_solute_type
                ],
            },
            "Inlet": {
                "pH": self.ro_feed.pH,
                "Pressure": pyunits.convert(
                    self.ro_unit.inlet.pressure[0], to_units=pyunits.bar
                ),
                "Temperature": self.ro_feed.inlet.temperature[0],
                "Flow rate": self.ro_feed.properties_out[0].flow_vol_phase["Liq"],
                self.ro_solute_type: self.ro_feed.properties_out[
                    0
                ].conc_mass_phase_comp["Liq", self.ro_solute_type],
            },
            "Retentate": {
                "pH": self.ro_retentate.pH,
                "Pressure": pyunits.convert(
                    self.ro_retentate.outlet.pressure[0], to_units=pyunits.bar
                ),
                "Temperature": self.ro_retentate.outlet.temperature[0],
                "Flow rate": self.ro_retentate.properties_in[0].flow_vol_phase["Liq"],
                self.ro_solute_type: self.ro_retentate.properties_in[
                    0
                ].conc_mass_phase_comp["Liq", self.ro_solute_type],
            },
            "Permeate": {
                "pH": self.ro_product.pH,
                "Pressure": pyunits.convert(
                    self.ro_product.outlet.pressure[0], to_units=pyunits.bar
                ),
                "Temperature": self.ro_product.outlet.temperature[0],
                "Flow rate": self.ro_product.properties_in[0].flow_vol_phase["Liq"],
                self.ro_solute_type: self.ro_product.properties_in[
                    0
                ].conc_mass_phase_comp["Liq", self.ro_solute_type],
            },
        }
        if self.config.add_reaktoro_chemistry:
            model_state_dict["Scaling potential"] = {}
            model_state_dict["Maximum scaling potential"] = {}
            for scalant in self.ro_unit.scaling_tendency:
                model_state_dict["Scaling potential"][scalant] = (
                    self.ro_unit.scaling_tendency[scalant]
                )
                model_state_dict["Maximum scaling potential"][scalant] = (
                    self.ro_unit.maximum_scaling_tendency[scalant]
                )
        if self.ro_unit.find_component("water_removed_at_interface") is not None:
            model_state_dict["RO operation"][
                "Water Removed"
            ] = self.ro_unit.water_removed_at_interface
        if self.config.default_costing_package is not None:
            model_state_dict["Captial cost"] = self.ro_unit.costing.capital_cost
            model_state_dict["Fixed operating cost"] = (
                self.ro_unit.costing.fixed_operating_cost
            )
        return self.name, model_state_dict

from idaes.core import (
    declare_process_block_class,
)

from pyomo.common.config import ConfigValue
from idaes.models.unit_models import (
    Feed,
)
from idaes.core import (
    FlowsheetBlock,
)

from pyomo.environ import (
    ConcreteModel,
    Var,
    value,
    Constraint,
    units as pyunits,
)
from watertap.core.solvers import get_solver
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.util.initialization import interval_initializer
import idaes.core.util.scaling as iscale
from reaktoro_pse.reaktoro_block import ReaktoroBlock
import idaes.logger as idaeslog
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.reaktoro_utils import (
    ReaktoroOptionsContainer,
)

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander Dudchenko"


@declare_process_block_class("MultiCompFeed")
class MultiCompFeedData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "ion_concentrations",
        ConfigValue(
            default=None,
            description="Feed composition used for MCAS",
            doc="""
            Provide a dictionary that contains apparent species, and their concentrations, each value should include a corresponding pyomo unit (etc.):
            {'Ca':50*pyunits.mg/pyunits.L,
            'Na':50*pyunits.mg/pyunits.L}
            """,
        ),
    )
    CONFIG.declare(
        "mass_flowrate",
        ConfigValue(
            default=1,
            description="mass flow rate of feed",
            doc="""
                Provide mass flowrate of feed with pyomo units
            """,
        ),
    )
    CONFIG.declare(
        "volumetric_flowrate",
        ConfigValue(
            default=None,
            description="volumetric flow rate of feed",
            doc="""
                Provide mass flowrate of feed with pyomo units
            """,
        ),
    )
    CONFIG.declare(
        "temperature",
        ConfigValue(
            default=293.15,
            description="Temperature of feed (must be in Kelvin)",
            doc="""
                Temperature of feed (must be in Kelvin)
            """,
        ),
    )
    CONFIG.declare(
        "pressure",
        ConfigValue(
            default=1 * pyunits.atm,
            description="Pressure of feed (must be in Pa or with py units)",
            doc="""
                Pressure of feed (must be in Pa or with py units)
            """,
        ),
    )
    CONFIG.declare(
        "pH",
        ConfigValue(
            default=7,
            description="pH of feed",
            doc="""
                Provide feed pH
            """,
        ),
    )
    CONFIG.declare(
        "pE",
        ConfigValue(
            default=None,
            description="pE of feed",
            doc="""
                Provide feed pE, specify either a float value, or True (bool) to have Reaktoro reconcile it for you. 
            """,
        ),
    )
    CONFIG.declare(
        "alkalinity_as_CaCO3",
        ConfigValue(
            default=None,
            description="Alkalinity target for the feed",
            doc="""
                If value is not None, will use it to reconcile alkalinity, should be provided 
                on basis of CaCO3 (mg/L)
            """,
        ),
    )
    CONFIG.declare(
        "reconcile_using_reaktoro",
        ConfigValue(
            default=True,
            description="Will reconcile feed using reaktoro",
            doc="""
                Reconciles feed for charge and alkalinity using reaktoro.
            """,
        ),
    )
    CONFIG.declare(
        "charge_balance_ion",
        ConfigValue(
            default="Cl_-",
            description="Ion to use for charge balancing",
            doc="""
            Ion that will be adjusted to reach charge balance of zero
            """,
        ),
    )
    CONFIG.declare(
        "alkalinity_balance_ions",
        ConfigValue(
            default="HCO3_-",
            description="Apparent species that should be adjusted to reach target alkalinity",
            doc="""
            Apparent species that should be adjusted to reach target alkalinity
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

        self.feed = Feed(property_package=self.config.default_property_package)
        self.feed.properties[0].flow_vol_phase["Liq"]

        self.feed.pH = Var(
            initialize=self.config.pH, bounds=(0, 13), units=pyunits.dimensionless
        )
        fs_vars = {"pH": self.feed.pH}
        if self.config.pE is not None and self.config.pE != False:
            self.feed.pE = Var(
                initialize=1,
                bounds=(None, None),
                units=pyunits.dimensionless,
            )
            fs_vars["pE"] = self.feed.pE
        self.register_port("outlet", self.feed.outlet, fs_vars)
        self.feed.alkalinity_as_CaCO3 = Var(units=pyunits.mg / pyunits.L)
        self.feed.alkalinity_as_CaCO3.fix()
        if self.config.alkalinity_as_CaCO3 is not None:
            self.feed.alkalinity_as_CaCO3.fix(self.config.alkalinity_as_CaCO3)

    def set_fixed_operation(self, solver=None, alt_block=None):
        """sets fixed unit operation"""
        if alt_block is None:
            block = self.feed
        else:
            block = alt_block
        for ion, value in self.config.ion_concentrations.items():
            block.properties[0].conc_mass_phase_comp["Liq", ion].fix(value)
        if self.config.volumetric_flowrate is not None:
            block.properties[0].flow_vol_phase["Liq"].fix(
                self.config.volumetric_flowrate
            )
        else:
            block.total_mass_flow = Constraint(
                expr=self.config.mass_flowrate
                == sum(
                    [
                        block.properties[0].flow_mass_phase_comp[phase, ion]
                        for phase, ion in block.properties[0].flow_mass_phase_comp
                    ]
                )
            )
            block.properties[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ] = self.config.mass_flowrate  # approximate
        if self.config.pE is not None and isinstance(self.config.pE, bool) == False:
            block.pE.fix(self.config.pE)
        elif self.config.pE == True:
            block.pE.fix()
        block.properties[0].pressure.fix(self.config.pressure)
        block.properties[0].temperature.fix(self.config.temperature)
        block.pH.fix(self.config.pH)

        # Doing a solve with out scaling -this in general works
        # as we have a simple problem, but solution might be sub optimal
        # this will still give a good initial state however to move forward

        interval_initializer(block)
        if solver is None:
            solver = get_solver()
        result = solver.solve(block, tee=False)
        assert_optimal_termination(result)
        # self.feed.display()
        assert degrees_of_freedom(block) == 0

    def reaktoro_reconciliation(self):
        sub_model = ConcreteModel()
        sub_model.fs = FlowsheetBlock()

        # create sub model for feed
        sub_model.fs.feed = Feed(property_package=self.config.default_property_package)
        sub_model.fs.feed.pH = Var(units=pyunits.dimensionless)
        iscale.set_scaling_factor(sub_model.fs.feed.pH, 1 / 10)
        if self.config.pE is not None and self.config.pE is True:
            sub_model.fs.feed.pE = Var(
                initialize=0,
                units=pyunits.dimensionless,
            )
            sub_model.fs.feed.pE.fix()
            iscale.set_scaling_factor(sub_model.fs.feed.pE, 1 / 10)
        sub_model.fs.feed.alkalinity_as_CaCO3 = Var(units=pyunits.mg / pyunits.L)
        iscale.set_scaling_factor(sub_model.fs.feed.alkalinity_as_CaCO3, 1 / 10)
        # set its conditions
        self.set_fixed_operation(alt_block=sub_model.fs.feed)

        iscale.calculate_scaling_factors(sub_model)
        # make sure we scale before solve

        sub_model.fs.feed.charge = Var(units=pyunits.dimensionless)
        iscale.set_scaling_factor(sub_model.fs.feed.charge, 1)
        if self.feed.find_component("charge") is None:
            self.feed.charge = Var(units=pyunits.dimensionless)
        iscale.calculate_scaling_factors(sub_model.fs.feed)
        initial_con = value(
            pyunits.convert(
                sub_model.fs.feed.properties[0].conc_mass_phase_comp[
                    "Liq", self.config.charge_balance_ion
                ],
                to_units=pyunits.g / pyunits.L,
            )
        )

        outputs = {("charge", None): sub_model.fs.feed.charge}

        outputs[("alkalinityAsCaCO3", None)] = sub_model.fs.feed.alkalinity_as_CaCO3
        if self.config.pE is not None and self.config.pE is True:
            outputs[("pE", None)] = sub_model.fs.feed.pE
            sub_model.fs.feed.pE.unfix()
        sub_model.fs.feed.properties[0].conc_mass_phase_comp[
            "Liq", self.config.charge_balance_ion
        ].unfix()

        reaktoro_options = ReaktoroOptionsContainer()
        reaktoro_options.system_state_option(
            "temperature",
            sub_model.fs.feed.properties[0].temperature,
        )
        reaktoro_options.system_state_option(
            "pressure",
            sub_model.fs.feed.properties[0].pressure,
        )
        reaktoro_options.system_state_option("pH", sub_model.fs.feed.pH)
        if self.config.pE is not None and isinstance(self.config.pE, bool) == False:
            reaktoro_options.system_state_option("pE", sub_model.fs.feed.pE)
        reaktoro_options.aqueous_phase_option(
            "composition",
            sub_model.fs.feed.properties[0].flow_mol_phase_comp,
        )
        reaktoro_options["build_speciation_block"] = False
        reaktoro_options["assert_charge_neutrality"] = False
        reaktoro_options["outputs"] = outputs
        reaktoro_options.update_with_user_options(self.config.reaktoro_options)
        reaktoro_options["reaktoro_block_manager"] = (
            None  # ensure we dont add this to parallel manager!
        )
        sub_model.fs.feed.charge_balance_block = ReaktoroBlock(**reaktoro_options)
        sub_model.fs.feed.charge_balance_block.initialize()

        sub_model.fs.feed.charge.fix(0)
        if self.config.alkalinity_as_CaCO3 is not None:
            sub_model.fs.feed.alkalinity_as_CaCO3.fix(self.config.alkalinity_as_CaCO3)
            if isinstance(self.config.alkalinity_balance_ions, str):
                self.config.alkalinity_balance_ions = [
                    self.config.alkalinity_balance_ions
                ]
            initial_alk_ions = {}
            for ion in self.config.alkalinity_balance_ions:
                sub_model.fs.feed.properties[0].conc_mass_phase_comp["Liq", ion].unfix()
                sub_model.fs.feed.properties[0].conc_mass_phase_comp["Liq", ion].setub(
                    5
                )
                initial_alk_ions[ion] = (
                    sub_model.fs.feed.properties[0]
                    .conc_mass_phase_comp["Liq", ion]
                    .value
                )
        else:
            sub_model.fs.feed.alkalinity_as_CaCO3.unfix()
        assert degrees_of_freedom(sub_model.fs.feed) == 0
        solver = get_cyipopt_watertap_solver()
        _log.info("Reconciling feed with reaktoro")
        result = solver.solve(sub_model.fs.feed, tee=True)
        if self.config.pE is not None and self.config.pE is True:
            sub_model.fs.feed.pE.fix()
        assert_optimal_termination(result)

        def replace_name(name, old_model, new_model):
            return name.replace(old_model, new_model)

        for v in sub_model.fs.component_data_objects(Var):
            v_n = replace_name(v.name, "fs.feed", "feed")
            self.find_component(v_n).value = v.value
            if v.fixed:
                self.find_component(v_n).fix()
            else:
                self.find_component(v_n).unfix()

        _log.info("Reconciliation complete")
        balanced_con = value(
            pyunits.convert(
                sub_model.fs.feed.properties[0].conc_mass_phase_comp[
                    "Liq", self.config.charge_balance_ion
                ],
                to_units=pyunits.g / pyunits.L,
            )
        )
        _log.info("Starting alkalinity reconciliation report")
        _log.info(
            f"Charge balanced, current charge is {sub_model.fs.feed.charge.value}"
        )
        _log.info(
            f"Increased {self.config.charge_balance_ion} from {initial_con} to {balanced_con} g/L)"
        )
        _log.info(f"Reconciled alkalinity is {self.feed.alkalinity_as_CaCO3.value}")
        if self.config.pE is not None and self.config.pE is True:
            _log.info(f"Reconciled pE is {self.feed.pE.value}")
        self.feed.properties[0].conc_mass_phase_comp[
            "Liq", self.config.charge_balance_ion
        ].fix()
        if self.config.alkalinity_as_CaCO3 is not None:
            for ion in self.config.alkalinity_balance_ions:
                self.feed.properties[0].conc_mass_phase_comp["Liq", ion].fix()
                _log.info(
                    f"Changed {ion} from {initial_alk_ions[ion]} to {self.feed.properties[0].conc_mass_phase_comp['Liq', ion].value} g/L"
                )
        self.feed.alkalinity_as_CaCO3.fix()
        solver.solve(self.feed, tee=False)
        _log.info(f"Report complete: DOFs {degrees_of_freedom(self)}")
        assert_optimal_termination(result)
        assert degrees_of_freedom(self) == 0

    def scale_before_initialization(self, **kwargs):
        iscale.set_scaling_factor(self.feed.pH, 1 / 10)
        if self.find_component("pE") is not None:
            iscale.set_scaling_factor(self.feed.pE, 1 / 10)
        self.scale_feed()

    def scale_feed(self):
        for idx in self.feed.properties[0].flow_mol_phase_comp:
            scale = 1 / self.feed.properties[0].flow_mol_phase_comp[idx].value
            self.config.default_property_package.set_default_scaling(
                "flow_mol_phase_comp", scale, index=idx
            )
            _log.info(f"Applied scaling factor to {idx} of {scale}")

    def initialize_unit(self, solver=None, tee=True):
        if solver is None:
            solver = get_solver()
        result = solver.solve(self.feed, tee=tee)
        assert_optimal_termination(result)
        if self.config.reconcile_using_reaktoro:
            self.reaktoro_reconciliation()
            if self.feed.find_component("total_mass_flow") != None:
                self.feed.total_mass_flow.deactivate()
                self.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix()
        assert degrees_of_freedom(self) == 0

    def get_model_state_dict(self):
        model_state = {
            "Composition": {},
            "Physical state": {},
        }
        model_state["Composition"]["Mass flow of H2O"] = self.feed.properties[
            0
        ].flow_mass_phase_comp["Liq", "H2O"]
        for phase, ion in self.feed.properties[0].conc_mass_phase_comp:
            if ion != "H2O":
                model_state["Composition"][ion] = self.feed.properties[
                    0
                ].conc_mass_phase_comp[phase, ion]
        model_state["Composition"]["pH"] = self.feed.pH
        if self.config.pE is not None and self.config.pE != False:
            model_state["Composition"]["pE"] = self.feed.pE
        model_state["Composition"]["Alkalinity"] = self.feed.alkalinity_as_CaCO3
        model_state["Physical state"]["Temperature"] = self.feed.properties[
            0
        ].temperature
        model_state["Physical state"]["Pressure"] = self.feed.properties[0].pressure
        model_state["Physical state"]["Volumetric flowrate"] = self.feed.properties[
            0
        ].flow_vol_phase["Liq"]
        return self.name, model_state

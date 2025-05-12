from idaes.core import (
    declare_process_block_class,
)

from pyomo.common.config import ConfigValue
from idaes.models.unit_models import (
    Feed,
)

from pyomo.environ import (
    Var,
    value,
    Constraint,
    units as pyunits,
)
from watertap.core.solvers import get_solver
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.cyipot_solver import (
    get_cyipopt_solver,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.util.initialization import interval_initializer
import idaes.core.util.scaling as iscale
from reaktoro_pse.reaktoro_block import ReaktoroBlock
import idaes.logger as idaeslog

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
        "charge_balance_with_reaktoro",
        ConfigValue(
            default=False,
            description="To charge balance with Reaktoro",
            doc="""
                To use Reaktoro to charge balance feed during initialization
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

        self.feed = Feed(property_package=self.config.default_property_package)

        self.feed.pH = Var(initialize=self.config.pH, units=pyunits.dimensionless)
        self.register_port("outlet", self.feed.outlet, {"pH": self.feed.pH})

    def set_fixed_operation(self, solver=None):
        """sets fixed unit operation"""

        for ion, value in self.config.ion_concentrations.items():
            self.feed.properties[0].conc_mass_phase_comp["Liq", ion].fix(value)

        if self.config.volumetric_flowrate is not None:
            self.feed.properties[0].flow_vol_phase["Liq"].fix(
                self.config.volumetric_flowrate
            )
        else:
            self.feed.total_mass_flow = Constraint(
                expr=self.config.mass_flowrate
                == sum(
                    [
                        self.feed.properties[0].flow_mass_phase_comp[phase, ion]
                        for phase, ion in self.feed.properties[0].flow_mass_phase_comp
                    ]
                )
            )
            self.feed.properties[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ] = self.config.mass_flowrate  # aproximate
        self.feed.properties[0].pressure.fix(self.config.pressure)
        self.feed.properties[0].temperature.fix(self.config.temperature)
        self.feed.pH.fix(self.config.pH)

        # Doing a solve with out scaling -this in general works
        # as we have a simple problem, but soluton might be sub optimal
        # this will still give a good intial state however to move forward
        interval_initializer(self.feed)
        if solver is None:
            solver = get_solver()
        result = solver.solve(self.feed, tee=False)
        assert_optimal_termination(result)
        assert degrees_of_freedom(self.feed) == 0

    def balance_charge_with_reaktoro(self):
        self.feed.charge = Var(units=pyunits.dimensionless)
        initial_con = value(
            pyunits.convert(
                self.feed.properties[0].conc_mass_phase_comp[
                    "Liq", self.config.charge_balance_ion
                ],
                to_units=pyunits.g / pyunits.L,
            )
        )
        self.feed.properties[0].conc_mass_phase_comp[
            "Liq", self.config.charge_balance_ion
        ].unfix()

        self.feed.charge_balance_block = ReaktoroBlock(
            system_state={
                "temperature": self.feed.properties[0].temperature,
                "pressure": self.feed.properties[0].pressure,
                "pH": self.feed.pH,
            },
            aqueous_phase={
                "composition": self.feed.properties[
                    0
                ].flow_mol_phase_comp,  # This is the spices mass flow
                "convert_to_rkt_species": True,  # We can use default converter as its defined for default database (Phreeqc and pitzer)
                "activity_model": self.config.reaktoro_options[
                    "activity_model"
                ],  # Can provide a string, or Reaktoro initialized class
                "fixed_solvent_specie": "H2O",  # We need to define our aqueous solvent as we have to speciate the block
            },
            outputs={("charge", None): self.feed.charge},
            assert_charge_neutrality=False,
            database=self.config.reaktoro_options[
                "database"
            ],  # can also be reaktoro.PhreeqcDatabase('pitzer.dat')
            database_file=self.config.reaktoro_options["database_file"],
        )
        self.feed.charge_balance_block.initialize()

        self.feed.charge.fix(0)

        assert degrees_of_freedom(self) == 0
        solver = get_cyipopt_solver()
        result = solver.solve(self.feed)
        _log.info(f"Charge balanced, current charge is {self.feed.charge.value}")
        balanced_con = value(
            pyunits.convert(
                self.feed.properties[0].conc_mass_phase_comp[
                    "Liq", self.config.charge_balance_ion
                ],
                to_units=pyunits.g / pyunits.L,
            )
        )
        _log.info(
            f"Increased {self.config.charge_balance_ion} from {initial_con} to {balanced_con} g/L)"
        )
        for v in self.feed.charge_balance_block.component_data_objects(Constraint):
            v.deactivate()
        self.feed.charge_balance_block.reaktoro_model.deactivate()
        self.feed.charge_balance_block.deactivate()

        self.feed.properties[0].conc_mass_phase_comp[
            "Liq", self.config.charge_balance_ion
        ].fix()
        assert_optimal_termination(result)
        assert degrees_of_freedom(self) == 0

    def scale_before_initialization(self, **kwargs):
        iscale.set_scaling_factor(self.feed.pH, 1)
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
        if self.config.charge_balance_with_reaktoro:
            self.balance_charge_with_reaktoro()
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
        model_state["Physical state"]["Temperature"] = self.feed.properties[
            0
        ].temperature
        model_state["Physical state"]["Pressure"] = self.feed.properties[0].pressure
        return self.name, model_state

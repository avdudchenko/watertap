__author__ = "Alexander Dudchenko"

from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import (
    assert_optimal_termination,
)

from pyomo.environ import (
    Var,
    units as pyunits,
)
from pyomo.common.config import ConfigValue
from watertap.unit_models.pressure_changer import Pump
from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock

import idaes.core.util.scaling as iscale


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
        "osmotic_over_pressure",
        ConfigValue(
            default=2,
            description="Value to multiply osmotic pressure by",
            doc="""
            Value to multiply osmotic pressure by if initialization_pressure is set to osmotic_pressure:
            guess pressure = osmotic pressure * osmotic overpressure + osmotic pressure offset
            """,
        ),
    )
    CONFIG.declare(
        "osmotic_pressure_offset",
        ConfigValue(
            default=2e5,
            description="Offset for osmotic pressure ",
            doc="""
            Sets an fixed offset of guessed initial pressure
            """,
        ),
    )
    CONFIG.declare(
        "maximum_pressure",
        ConfigValue(
            default=300e5,
            description="Maximum pressure for pump unit",
            doc="""
            Sets the maximum pressure for the pump unit
            """,
        ),
    )
    CONFIG.declare(
        "pump_efficiency",
        ConfigValue(
            default=0.8,
            description="default pump efficiency",
            doc="""
            default pump efficiency
            """,
        ),
    )

    def build(self):
        """Build a multi-component pump unit model"""
        super().build()
        self.pump = Pump(property_package=self.config.default_property_package)
        if self.config.default_costing_package is not None:
            self.pump.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package
            )
        # Add pump flow rate
        self.pump.control_volume.properties_in[0].flow_vol_phase[...]
        self.pump.pH = Var(initialize=7, bounds=(1, 12), units=pyunits.dimensionless)

        self.register_port("inlet", self.pump.inlet, {"pH": self.pump.pH})
        self.register_port("outlet", self.pump.outlet, {"pH": self.pump.pH})

    def set_fixed_operation(self):
        """fixes operation point for pump unit model
        Uses osmotic pressure to initialize pump outlet pressure or user defined pressure
        """
        if self.config.initialization_pressure == "osmotic_pressure":
            init_flags = self.pump.control_volume.initialize()
            self.pump.control_volume.release_state(init_flags)
            self.pump.outlet.pressure[0].fix(
                self.pump.control_volume.properties_in[0]
                .pressure_osm_phase["Liq"]
                .value
                * self.config.osmotic_over_pressure
                + self.config.osmotic_pressure_offset
            )
        else:
            self.pump.outlet.pressure[0].fix(self.config.initialization_pressure)
        self.pump.outlet.pressure[0].setub(self.config.maximum_pressure)
        self.pump.efficiency_pump[0].fix(self.config.pump_efficiency)
        self.inlet.fix()
        assert degrees_of_freedom(self) == 0
        self.inlet.unfix()

    def set_optimization_operation(self):
        self.pump.outlet.pressure[0].unfix()

    def scale_before_initialization(self, **kwargs):
        iscale.set_scaling_factor(self.pump.outlet.pressure, 1e-5)
        iscale.set_scaling_factor(self.pump.control_volume.work, 1e-4)
        iscale.set_scaling_factor(self.pump.pH, 1 / 1000)

    def initialize_unit(self):
        self.set_fixed_operation()
        self.pump.initialize()

    def get_model_state_dict(self):
        """Returns a dictionary with the model state"""

        unit_dofs = degrees_of_freedom(self)

        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "Overall": {
                "pH": self.pump.pH,
                "Temperature": self.pump.inlet.temperature[0],
                "Flow rate": self.pump.control_volume.properties_in[0].flow_vol_phase[
                    "Liq"
                ],
            },
            "Inlet": {
                "Pressure": self.pump.inlet.pressure[0],
            },
            "Outlet": {
                "Pressure": self.pump.outlet.pressure[0],
            },
        }
        if self.config.default_costing_package is not None:
            model_state_dict["Costs"] = {"Capital cost": self.pump.costing.capital_cost}
        return self.name, model_state_dict

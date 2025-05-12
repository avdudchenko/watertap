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
from watertap.unit_models.pressure_changer import EnergyRecoveryDevice
from idaes.core import (
    declare_process_block_class,
)
from idaes.core import UnitModelCostingBlock

import idaes.core.util.scaling as iscale


@declare_process_block_class("MultiCompERDUnit")
class MultiCompERDUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "erd_efficiency",
        ConfigValue(
            default=0.9,
            description="default energy recovery device efficiency",
            doc="""
            default energy recovery device efficiency
            """,
        ),
    )
    CONFIG.declare(
        "erd_outlet_pressure",
        ConfigValue(
            default=1 * pyunits.atm,
            description="ERD outlet pressure",
            doc="""
            Energy recovery device pressure

            """,
        ),
    )

    def build(self):
        """Build a multi-component ERD unit model"""
        super().build()
        self.ERD = EnergyRecoveryDevice(
            property_package=self.config.default_property_package
        )
        if self.config.default_costing_package is not None:
            self.ERD.costing = UnitModelCostingBlock(
                flowsheet_costing_block=self.config.default_costing_package
            )
        # Add ERD flow rate
        self.ERD.control_volume.properties_in[0].flow_vol_phase[...]
        self.ERD.pH = Var(initialize=7, units=pyunits.dimensionless)

        self.register_port("inlet", self.ERD.inlet, {"pH": self.ERD.pH})
        self.register_port("outlet", self.ERD.outlet, {"pH": self.ERD.pH})

    def set_fixed_operation(self):
        """fixes operation point for ERD unit model
        Uses osmotic pressure to initialize ERD outlet pressure or user defined pressure
        """
        self.ERD.efficiency_pump[0].fix(self.config.erd_efficiency)
        self.ERD.control_volume.properties_out[0].pressure.fix(
            self.config.erd_outlet_pressure
        )
        self.inlet.fix()
        assert degrees_of_freedom(self) == 0
        self.inlet.unfix()

    def scale_before_initialization(self, **kwargs):
        iscale.set_scaling_factor(self.ERD.inlet.pressure, 1e-5)
        iscale.set_scaling_factor(self.ERD.outlet.pressure, 1e-5)
        iscale.set_scaling_factor(self.ERD.control_volume.work, 1e-4)

    def initialize_unit(self):
        self.ERD.initialize()

    def get_model_state_dict(self):
        """Returns a dictionary with the model state"""
        self.inlet.fix()
        unit_dofs = degrees_of_freedom(self)
        self.inlet.unfix()

        model_state_dict = {
            "Model": {"DOFs": unit_dofs},
            "Overall": {
                "pH": self.ERD.pH,
                "Temperature": self.ERD.inlet.temperature[0],
                "Flow rate": self.ERD.control_volume.properties_in[0].flow_vol_phase[
                    "Liq"
                ],
            },
            "Inlet": {
                "Pressure": self.ERD.inlet.pressure[0],
            },
            "Outlet": {
                "Pressure": self.ERD.outlet.pressure[0],
            },
        }
        self.config.default_costing_package.display()
        if self.config.default_costing_package is not None:
            model_state_dict["Costs"] = {
                "Capital cost": self.ERD.costing.capital_cost,
                # "Operating cost": self.ERD.costing.operating_cost,
            }
        return self.name, model_state_dict

from idaes.core import (
    declare_process_block_class,
)

from idaes.models.unit_models import (
    Product,
)

from pyomo.environ import (
    Var,
    units as pyunits,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)

from pyomo.common.config import ConfigValue
import idaes.core.util.scaling as iscale

__author__ = "Alexander Dudchenko"


@declare_process_block_class("MultiCompProduct")
class MultiCompProductData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
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

        self.product = Product(property_package=self.config.default_property_package)

        self.product.pH = Var(initialize=7, bounds=(0, 13), units=pyunits.dimensionless)
        inlet_vars = {"pH": self.product.pH}
        if self.config.track_pE:
            self.product.pE = Var(
                initialize=0,
                units=pyunits.dimensionless,
                bounds=(None, None),
            )
            inlet_vars["pE"] = self.product.pE

        self.register_port("inlet", self.product.inlet, inlet_vars)
        self.product.properties[0].conc_mass_phase_comp[...]
        self.product.properties[0].flow_mass_phase_comp[...]

    def scale_before_initialization(self, **kwargs):
        iscale.set_scaling_factor(self.product.pH, 1 / 10)
        if self.config.track_pE:
            iscale.set_scaling_factor(self.product.pE, 1 / 10)

    def initialize_unit(self, solver=None, tee=True):
        self.product.initialize()

    def get_model_state_dict(self):
        model_state = {
            "Composition": {},
            "Physical state": {},
        }
        model_state["Composition"]["Mass flow of H2O"] = self.product.properties[
            0
        ].flow_mass_phase_comp["Liq", "H2O"]
        for phase, ion in self.product.properties[0].conc_mass_phase_comp:
            if ion != "H2O":
                model_state["Composition"][ion] = self.product.properties[
                    0
                ].conc_mass_phase_comp[phase, ion]
        model_state["Composition"]["pH"] = self.product.pH
        model_state["Physical state"]["Temperature"] = self.product.properties[
            0
        ].temperature
        model_state["Physical state"]["Pressure"] = self.product.properties[0].pressure
        return self.name, model_state

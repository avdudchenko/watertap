__author__ = "Alexander Dudchenko"
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)

from idaes.core.util.model_statistics import degrees_of_freedom


from pyomo.environ import (
    Var,
    units as pyunits,
)
from pyomo.common.config import ConfigValue

from idaes.models.unit_models import Separator
from idaes.core import (
    declare_process_block_class,
)

import idaes.core.util.scaling as iscale
from collections import OrderedDict


__author__ = "Alexander Dudchenko"


@declare_process_block_class("SplitterPhUnit")
class SplitterPhUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "outlet_ports",
        ConfigValue(
            default=["inlet_1", "inlet_2"],
            description="inlet ports to build on the splitter",
            doc="""
                Define a list of inlet ports to build on the splitter.
            """,
        ),
    )
    CONFIG.declare(
        "splitter_initialization_guess",
        ConfigValue(
            default=0.5,
            description="Sets split ratio",
            doc="""
               Sets initialization guess, if single value applies it to all outlet ports,except last one
               if dict, applies split ratio teach named port and fixes it.
               Always leave one port split ratio free to ensure DOFs=0.
        """,
        ),
    )

    CONFIG.declare(
        "splitter_options",
        ConfigValue(
            default=None,
            description="Options to pass to the splitter during build ",
            doc="""
               This will pass options to splitter and override default options:
                       defaults_splitter_props = {
            "energy_mixing_type": MixingType.none,
            "momentum_mixing_type": MomentumMixingType.minimize,
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
        defaults_splitter_props = {}
        if self.config.splitter_options is not None:
            defaults_splitter_props.update(self.config.splitter_options)
        self.splitter = Separator(
            property_package=self.config.default_property_package,
            outlet_list=self.config.outlet_ports,
            **defaults_splitter_props,
        )

        self.splitter.pH = Var(
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(0, 13),
        )
        if self.config.track_pE:
            self.splitter.pE = Var(
                initialize=0,
                units=pyunits.dimensionless,
                bounds=(None, None),
            )
        for port in self.config.outlet_ports:
            self.splitter.find_component(f"{port}_state")[0].flow_mass_phase_comp[...]
            self.splitter.find_component(f"{port}_state")[0].conc_mass_phase_comp[...]
            outlet_vars = {"pH": self.splitter.pH}
            if self.config.track_pE:
                outlet_vars["pE"] = self.splitter.pE
            self.register_port(
                port,
                self.splitter.find_component(port),
                outlet_vars,
            )
        inlet_vars = {"pH": self.splitter.pH}
        if self.config.track_pE:
            inlet_vars["pE"] = self.splitter.pE
        self.register_port(
            "inlet",
            self.splitter.inlet,
            inlet_vars,
        )
        self.splitter.find_component(f"mixed_state")[0].flow_mass_phase_comp[...]
        self.splitter.find_component(f"mixed_state")[0].conc_mass_phase_comp[...]

    def scale_before_initialization(self, **kwargs):
        iscale.set_scaling_factor(self.splitter.pH, 1 / 10)

    def set_fixed_operation(self):
        """fixes operation point for pump unit model
        Uses osmotic pressure to initialize pump outlet pressure or user defined pressure
        """
        if isinstance(self.config.splitter_initialization_guess, dict):
            for port, frac in self.config.splitter_initialization_guess.items():
                self.splitter.split_fraction[0, port].fix(frac)
        else:
            for port in self.config.outlet_ports[:-1]:
                self.splitter.split_fraction[0, port].fix(
                    self.config.splitter_initialization_guess
                )
        self.splitter.split_fraction[0, self.config.outlet_ports[0]].setlb(0.01)
        self.splitter.split_fraction[0, self.config.outlet_ports[0]].setub(0.99)

    def set_optimization_operation(self):
        for port in self.config.outlet_ports:
            self.splitter.split_fraction[0, port].unfix()

    def initialize_unit(self, **kwargs):
        self.splitter.initialize()

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

        unit_dofs = degrees_of_freedom(self)

        model_state = {
            "Model": {"DOFs": unit_dofs},
        }

        def get_pe():
            if self.config.track_pE:
                return self.splitter.pE
            else:
                return None

        model_state["Inlet"] = get_ion_comp(
            self.splitter.find_component("mixed_state")[0],
            self.splitter.pH,
            get_pe(),
        )
        for outlet in self.config.outlet_ports:
            model_state[outlet] = get_ion_comp(
                self.splitter.find_component(f"{outlet}_state")[0],
                self.splitter.pH,
                get_pe(),
            )

        return self.name, model_state

__author__ = "Alexander Dudchenko"
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.watertap_flowsheet_block import (
    WaterTapFlowsheetBlockData,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.reaktoro_utils import (
    ReaktoroOptionsContainer,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import (
    Var,
    Constraint,
    units as pyunits,
)
import math
from pyomo.common.config import ConfigValue
from idaes.models.unit_models import Mixer
from idaes.models.unit_models.mixer import MomentumMixingType, MixingType
from idaes.core import (
    declare_process_block_class,
)
import idaes.core.util.scaling as iscale
from reaktoro_pse.reaktoro_block import ReaktoroBlock
from collections import OrderedDict

from pyomo.util.calc_var_value import calculate_variable_from_constraint

__author__ = "Alexander Dudchenko"


@declare_process_block_class("MixerPhUnit")
class MixerPhUnitData(WaterTapFlowsheetBlockData):
    CONFIG = WaterTapFlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "inlet_ports",
        ConfigValue(
            default=["inlet_1", "inlet_2"],
            description="inlet ports to build on the mixer",
            doc="""
                Define a list of inlet ports to build on the mixer.
            """,
        ),
    )
    CONFIG.declare(
        "mixer_options",
        ConfigValue(
            default=None,
            description="Options to pass to the mixer during build ",
            doc="""
               This will pass options to mixer and override default options:
                       defaults_mixer_props = {
            "energy_mixing_type": MixingType.none,
            "momentum_mixing_type": MomentumMixingType.minimize,
        """,
        ),
    )
    CONFIG.declare(
        "guess_secondary_inlet_composition",
        ConfigValue(
            default=False,
            description="Set if secondary inlet composition should be guessed",
            doc="""
               Set this to True if you are using mixer in a recycle loop and don't have 
               the recycled stream composition. Otherwise to False, if composition from 
               secondary inlets is known and is propagated.
        """,
        ),
    )
    CONFIG.declare(
        "add_reaktoro_chemistry",
        ConfigValue(
            default=False,
            description="To use Reaktoro-PSE for estimate pH change",
            doc="""
            If True, builds a reaktoro block and uses it to calculate change in pH due to the addition of given chemical. 
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
        defaults_mixer_props = {
            "energy_mixing_type": MixingType.none,
            "momentum_mixing_type": MomentumMixingType.minimize,
        }
        if self.config.mixer_options is not None:
            defaults_mixer_props.update(self.config.mixer_options)
        self.mixer = Mixer(
            property_package=self.config.default_property_package,
            inlet_list=self.config.inlet_ports,
            **defaults_mixer_props,
        )
        all_ports = self.config.inlet_ports + ["outlet"]
        self.mixer.pH = Var(
            all_ports,
            initialize=7,
            units=pyunits.dimensionless,
            bounds=(0, 13),
        )

        if self.config.add_reaktoro_chemistry:
            self.add_reaktoro_chemistry()
        else:
            # flow average pH of all inlets
            self.mixer.eq_pH = Constraint(
                expr=sum(
                    self.mixer.pH[port]
                    * self.mixer.find_component(f"{port}_state")[0].flow_vol_phase[
                        "Liq"
                    ]
                    for port in self.config.inlet_ports
                )
                / self.mixer.find_component("mixed_state")[0].flow_vol_phase["Liq"]
                == self.mixer.pH["outlet"]
            )
        for port in self.config.inlet_ports:
            self.mixer.find_component(f"{port}_state")[0].flow_mass_phase_comp[...]
            self.mixer.find_component(f"{port}_state")[0].conc_mass_phase_comp[...]

            self.register_port(
                port,
                self.mixer.find_component(port),
                {"pH": self.mixer.pH[port]},
            )
        self.register_port(
            "outlet",
            self.mixer.outlet,
            {"pH": self.mixer.pH["outlet"]},
        )
        self.mixer.find_component(f"mixed_state")[0].flow_mass_phase_comp[...]
        self.mixer.find_component(f"mixed_state")[0].conc_mass_phase_comp[...]
        # make sure we pass temperature to outlet
        self.mixer.temp_constraint = Constraint(
            expr=self.mixer.find_component(f"{self.config.inlet_ports[0]}_state")[
                0
            ].temperature
            == self.mixer.mixed_state[0].temperature
        )
        self.mixer_initialized = False

    def add_reaktoro_chemistry(self):
        # for tracking all true species in the inlet streams
        mixing_blocks = []
        for port in self.config.inlet_ports[1:]:
            mixer_inlet = self.mixer.find_component(f"{port}_state")[0]
            reaktoro_options = ReaktoroOptionsContainer()
            reaktoro_options.system_state_option(
                "temperature",
                mixer_inlet.temperature,
            )
            reaktoro_options.system_state_option(
                "pressure",
                mixer_inlet.pressure,
            )
            reaktoro_options.system_state_option("pH", self.mixer.pH[port])
            reaktoro_options.aqueous_phase_option(
                "composition",
                mixer_inlet.flow_mol_phase_comp,
            )
            reaktoro_options["build_speciation_block"] = False

            reaktoro_options["outputs"] = {"speciesAmount": True}
            reaktoro_options["build_graybox_model"] = False
            reaktoro_options.update_with_user_options(self.config.reaktoro_options)

            self.add_component(
                f"{port}_speciation_block", ReaktoroBlock(**reaktoro_options)
            )
            mixing_blocks.append(self.find_component(f"{port}_speciation_block"))
        feed_port = self.mixer.find_component(f"{self.config.inlet_ports[0]}_state")[0]
        reaktoro_options = ReaktoroOptionsContainer()
        reaktoro_options.system_state_option(
            "temperature",
            feed_port.temperature,
        )
        reaktoro_options.system_state_option(
            "pressure",
            feed_port.pressure,
        )
        reaktoro_options.aqueous_phase_option(
            "composition",
            feed_port.flow_mol_phase_comp,
        )
        reaktoro_options.system_state_option(
            "pH", self.mixer.pH[self.config.inlet_ports[0]]
        )
        reaktoro_options["build_speciation_block"] = True
        reaktoro_options["outputs"] = {("pH", None): self.mixer.pH["outlet"]}
        reaktoro_options["external_speciation_reaktoro_blocks"] = mixing_blocks
        reaktoro_options.update_with_user_options(self.config.reaktoro_options)
        self.mixer_speciation_block = ReaktoroBlock(**reaktoro_options)

    def scale_before_initialization(self, **kwargs):
        iscale.constraint_scaling_transform(self.mixer.temp_constraint, 1e-2)
        for ph in self.mixer.pH:
            iscale.set_scaling_factor(self.mixer.pH[ph], 1 / 10)
        if self.config.add_reaktoro_chemistry == False:
            iscale.constraint_scaling_transform(self.mixer.eq_pH, 1 / 10)

    def initialize_streams(self, **kwargs):
        self.fixed_streams = []
        if (
            self.config.guess_secondary_inlet_composition
            and self.mixer_initialized == False
        ):
            stream_init = {}
            for i, inlet in enumerate(self.config.inlet_ports):
                inlet_var = self.mixer.find_component(f"{inlet}_state")[0]

                if i > 0:
                    stream_init[inlet] = True
                else:
                    stream_init[inlet] = False
                    ref_stream = inlet_var
            for inlet, init_state in stream_init.items():
                if init_state:
                    inlet_var = self.mixer.find_component(f"{inlet}_state")[0]
                    for idx, obj in inlet_var.flow_mol_phase_comp.items():
                        obj.fix(ref_stream.flow_mol_phase_comp[idx].value * 1)
                        self.fixed_streams.append(obj)
                    self.mixer.pH[inlet].value = self.mixer.pH[
                        f"{self.config.inlet_ports[0]}"
                    ].value
            self.mixer.mixed_state[0].temperature.value = self.mixer.find_component(
                f"{self.config.inlet_ports[0]}_state"
            )[0].temperature.value
            self.mixer_initialized = True

    def initialize_unit(self, **kwargs):
        self.mixer.pH.display()
        self.initialize_streams()
        self.mixer.pH.display()
        self.mixer.pH.fix()
        self.mixer.pH["outlet"].unfix()
        self.mixer.initialize()
        self.mixer.pH.unfix()
        self.mixer.pH.display()
        for obj in self.fixed_streams:
            obj.unfix()
        self.fixed_streams = []
        calculate_variable_from_constraint(
            self.mixer.mixed_state[0].temperature,
            self.mixer.temp_constraint,
        )
        if self.config.add_reaktoro_chemistry:
            self.mixer_speciation_block.initialize()

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

        unit_dofs = degrees_of_freedom(self)

        model_state = {
            "Model": {"DOFs": unit_dofs},
        }
        for inlet in self.config.inlet_ports:
            model_state[inlet] = get_ion_comp(
                self.mixer.find_component(f"{inlet}_state")[0],
                self.mixer.pH[inlet],
            )

        model_state["Outlet"] = get_ion_comp(
            self.mixer.find_component("mixed_state")[0],
            self.mixer.pH["outlet"],
        )
        return self.name, model_state

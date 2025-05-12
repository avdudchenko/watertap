import pytest
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_feed_unit import (
    MultiCompFeed,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_ph_mixer_unit import (
    MixerPhUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.water_sources.source_water_importer import (
    get_source_water_data,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.cyipot_solver import (
    get_cyipopt_solver,
)

from pyomo.environ import ConcreteModel
from idaes.core import (
    FlowsheetBlock,
)

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
from pyomo.environ import (
    TransformationFactory,
    units as pyunits,
)

__author__ = "Alexander Dudchenko"


@pytest.mark.component
def test_mixing_sea_brackish_water():
    mcas_props, USDA_feed_specs = get_source_water_data(
        f"../../water_sources/USDA_brackish.yaml"
    )
    _, sea_water_feed_specs = get_source_water_data(
        f"../../water_sources/Seawater.yaml"
    )
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.usda_feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        charge_balance_with_reaktoro=True,
        **USDA_feed_specs,
    )
    m.fs.sea_water_feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        charge_balance_with_reaktoro=True,
        **sea_water_feed_specs,
    )
    m.fs.usda_feed.fix_and_scale()
    m.fs.sea_water_feed.set_fixed_operation()

    m.fs.mixer = MixerPhUnit(
        default_property_package=m.fs.properties,
        inlet_ports=["usda_feed", "sea_water_feed"],
    )
    m.fs.mixer.fix_and_scale()
    m.fs.usda_feed.outlet.connect_to(m.fs.mixer.usda_feed)
    m.fs.sea_water_feed.outlet.connect_to(m.fs.mixer.sea_water_feed)
    TransformationFactory("network.expand_arcs").apply_to(m)
    assert degrees_of_freedom(m) == 0
    iscale.calculate_scaling_factors(m)

    m.fs.usda_feed.initialize()

    m.fs.sea_water_feed.initialize()
    m.fs.mixer.initialize()
    m.fs.mixer.report()
    assert degrees_of_freedom(m) == 0
    solver = get_cyipopt_solver(10)
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.mixer.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.mixer.mixer.pH["outlet"].value,
            1e-5,
        )
        == 7.1033
    )
    assert (
        pytest.approx(
            m.fs.mixer.mixer.mixed_state[0].pressure.value,
            1e-5,
        )
        == 101325
    )
    assert (
        pytest.approx(
            m.fs.mixer.mixer.mixed_state[0].temperature.value,
            1e-5,
        )
        == 293.15
    )

    assert (
        pytest.approx(
            m.fs.mixer.mixer.mixed_state[0].flow_mass_phase_comp["Liq", "H2O"].value,
            1e-5,
        )
        == 1.9622
    )

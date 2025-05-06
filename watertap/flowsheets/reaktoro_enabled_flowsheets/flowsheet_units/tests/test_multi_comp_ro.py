import pytest
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_pump_unit import (
    MultiCompPumpUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_ro import (
    MultiCompROUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.tests.test_multicomponent_feed import (
    build_case,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.cyipot_solver import (
    get_cyipopt_solver,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
from pyomo.environ import (
    TransformationFactory,
    value,
    units as pyunits,
)

from watertap.costing import WaterTAPCosting
import watertap.property_models.seawater_prop_pack as sea_water_props

__author__ = "Alexander Dudchenko"


@pytest.mark.component
def test_osmotic_init_pressure():
    m = build_case("USDA_brackish", True)
    m.fs.sea_water_prop_pack = sea_water_props.SeawaterParameterBlock()
    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        initialization_pressure="osmotic_pressure",
    )
    m.fs.ro_unit = MultiCompROUnit(
        default_property_package=m.fs.properties,
        ro_property_package=m.fs.sea_water_prop_pack,
    )

    m.fs.feed.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ro_unit.feed)

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.pump_unit.fix_and_scale()
    m.fs.ro_unit.fix_and_scale()

    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.pump_unit.initialize()

    m.fs.pump_unit.report()
    m.fs.ro_unit.initialize()
    m.fs.ro_unit.report()
    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.ro_unit.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.ro_unit.ro_unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"],
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.hr),
                )
            ),
            1e-2,
        )
        == 4.4490
    )

    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.ro_unit.ro_unit.inlet.pressure[0], to_units=pyunits.bar
                )
            ),
            1e-2,
        )
        == 6.1710
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.scaling_tendency["Calcite"]),
            1e-2,
        )
        == 1.4353
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.scaling_tendency["Gypsum"]),
            1e-2,
        )
        == 0.35895951780356383
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_product.pH),
            1e-2,
        )
        == 7.0674
    )


@pytest.mark.component
def test_user_options():
    m = build_case("USDA_brackish", True)
    m.fs.sea_water_prop_pack = sea_water_props.SeawaterParameterBlock()
    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        initialization_pressure="osmotic_pressure",
    )
    m.fs.ro_unit = MultiCompROUnit(
        default_property_package=m.fs.properties,
        ro_property_package=m.fs.sea_water_prop_pack,
        selected_scalants={"Brucite": 1, "Calcite": 1},
        build_monotonic_cp_constraint=False,
        ro_options_dict={"finite_elements": 5},
    )

    m.fs.feed.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ro_unit.feed)

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.pump_unit.fix_and_scale()
    m.fs.ro_unit.fix_and_scale()

    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.pump_unit.initialize()

    m.fs.pump_unit.report()
    m.fs.ro_unit.initialize()
    m.fs.ro_unit.report()
    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.ro_unit.report()
    assert degrees_of_freedom(m) == 0
    assert (
        len(m.fs.ro_unit.ro_unit.length_domain) == 6
    )  # we have 5 finite elements + 1 for 0th node

    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.ro_unit.ro_unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"],
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.hr),
                )
            ),
            1e-2,
        )
        == 4.374655086
    )

    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.ro_unit.ro_unit.inlet.pressure[0], to_units=pyunits.bar
                )
            ),
            1e-2,
        )
        == 6.1710
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.scaling_tendency["Calcite"]),
            1e-2,
        )
        == 1.4353
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.scaling_tendency["Brucite"]),
            1e-2,
        )
        == 8.15764261602e-07
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_product.pH),
            1e-2,
        )
        == 7.0674
    )


@pytest.mark.component
def test_costing():
    m = build_case("USDA_brackish", True)
    m.fs.costing = WaterTAPCosting()
    m.fs.sea_water_prop_pack = sea_water_props.SeawaterParameterBlock()
    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        initialization_pressure="osmotic_pressure",
    )
    m.fs.ro_unit = MultiCompROUnit(
        default_property_package=m.fs.properties,
        ro_property_package=m.fs.sea_water_prop_pack,
        default_costing_package=m.fs.costing,
    )

    m.fs.feed.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ro_unit.feed)

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.pump_unit.fix_and_scale()
    m.fs.ro_unit.fix_and_scale()

    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.pump_unit.initialize()

    m.fs.pump_unit.report()
    m.fs.ro_unit.initialize()
    m.fs.ro_unit.report()
    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.ro_unit.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.ro_unit.ro_unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"],
                    to_units=pyunits.kg / (pyunits.m**2 * pyunits.hr),
                )
            ),
            1e-2,
        )
        == 4.4490
    )

    assert (
        pytest.approx(
            value(
                pyunits.convert(
                    m.fs.ro_unit.ro_unit.inlet.pressure[0], to_units=pyunits.bar
                )
            ),
            1e-2,
        )
        == 6.1710
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.scaling_tendency["Calcite"]),
            1e-2,
        )
        == 1.4353
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.scaling_tendency["Gypsum"]),
            1e-2,
        )
        == 0.35895951780356383
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_product.pH),
            1e-2,
        )
        == 7.0674
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.costing.capital_cost),
            1e-2,
        )
        == 6000
    )
    assert (
        pytest.approx(
            value(m.fs.ro_unit.ro_unit.costing.fixed_operating_cost),
            1e-2,
        )
        == 600
    )

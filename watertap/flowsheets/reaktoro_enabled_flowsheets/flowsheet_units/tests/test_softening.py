import pytest
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.softening_unit import (
    SofteningUnit,
    ViablePrecipitants,
    ViableReagents,
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
    units as pyunits,
)

__author__ = "Alexander Dudchenko"


@pytest.mark.component
def test_softening_default():
    m = build_case("USDA_brackish", True)
    m.fs.softening = SofteningUnit(
        default_property_package=m.fs.properties,
    )
    m.fs.softening.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.softening.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.softening.initialize()
    m.fs.softening.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.softening.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.softening.softening_reactor.flow_mass_precipitate["Calcite"].value,
            1e-5,
        )
        == 4.2935e-05
    )


@pytest.mark.component
def test_softening_with_all_options():
    m = build_case("USDA_brackish", True)
    m.fs.softening = SofteningUnit(
        default_property_package=m.fs.properties,
        selected_precipitants=ViablePrecipitants().keys(),
        selected_reagents=ViableReagents().keys(),
    )
    m.fs.softening.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.softening.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.softening.initialize()
    m.fs.softening.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.softening.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.softening.softening_reactor.flow_mass_precipitate["Calcite"].value,
            1e-5,
        )
        == 4.2935e-05
    )
    assert (
        pytest.approx(
            m.fs.softening.softening_reactor.flow_mass_precipitate["Gypsum"].value,
            1e-5,
        )
        == 6.9149e-30
    )
    assert (
        pytest.approx(
            m.fs.softening.softening_reactor.flow_mass_precipitate["Brucite"].value,
            1e-5,
        )
        == 6.9149e-30
    )


@pytest.mark.component
def test_softening_with_custom_options():
    m = build_case("USDA_brackish", True)
    viable_precipitants = ViablePrecipitants()
    viable_precipitants.register_solid(
        "Anhydrite",
        172.17 * pyunits.g / pyunits.mol,
        {"Ca_2+": 1, "SO4_2-": 1},
        "Ca_2+",
    )
    viable_precipitants.register_solid(
        "Aragonite",
        100.09 * pyunits.g / pyunits.mol,
        {"Ca_2+": 1, "HCO3_-": 1},
        "Ca_2+",
    )
    viable_reagents = ViableReagents()
    viable_reagents.register_reagent(
        "NaOH", 39.9971 * pyunits.g / pyunits.mol, {"Na_+": 1, "H2O": 1}
    )
    m.fs.softening = SofteningUnit(
        default_property_package=m.fs.properties,
        selected_precipitants=["Anhydrite", "Aragonite"],
        viable_reagents=viable_reagents,
        viable_precipitants=viable_precipitants,
        selected_reagents=["NaOH"],
    )
    m.fs.softening.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.softening.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.softening.initialize()
    m.fs.softening.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.softening.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.softening.softening_reactor.flow_mass_precipitate["Anhydrite"].value,
            1e-5,
        )
        == 3.6981e-26
    )
    assert (
        pytest.approx(
            m.fs.softening.softening_reactor.flow_mass_precipitate["Aragonite"].value,
            1e-5,
        )
        == 5.5670e-06
    )

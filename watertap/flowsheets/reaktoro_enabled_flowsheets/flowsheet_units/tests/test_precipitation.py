import pytest
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.precipitation_unit import (
    PrecipitationUnit,
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
from watertap.costing import WaterTAPCosting

__author__ = "Alexander Dudchenko"


@pytest.mark.component
def test_precip_default():
    m = build_case("USDA_brackish", True)
    m.fs.precipitation = PrecipitationUnit(
        default_property_package=m.fs.properties,
    )
    m.fs.precipitation.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.precipitation.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.precipitation.initialize()
    m.fs.precipitation.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.precipitation.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mass_precipitate[
                "Calcite"
            ].value,
            1e-5,
        )
        == 4.2935e-05
    )


@pytest.mark.component
def test_costing():
    m = build_case("USDA_brackish", True)
    m.fs.costing = WaterTAPCosting()
    m.fs.precipitation = PrecipitationUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
    )
    m.fs.precipitation.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.precipitation.inlet)
    m.fs.costing.cost_process()
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)
    m.fs.costing.initialize()
    m.fs.feed.initialize()
    m.fs.precipitation.initialize()
    m.fs.precipitation.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.precipitation.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.costing.capital_cost.value,
            1e-1,
        )
        == 2433.2
    )
    assert (
        pytest.approx(
            m.fs.costing.aggregate_flow_costs["fs_precipitation_reagent_CaO"].value,
            1e-1,
        )
        == 48.914
    )

    assert (
        pytest.approx(
            m.fs.costing.aggregate_flow_costs["fs_precipitation_reagent_Na2CO3"].value,
            1e-1,
        )
        == 59.959
    )


@pytest.mark.component
def test_precipitation_with_all_options():
    m = build_case("USDA_brackish", True)
    m.fs.precipitation = PrecipitationUnit(
        default_property_package=m.fs.properties,
        selected_precipitants=ViablePrecipitants().keys(),
        selected_reagents=ViableReagents().keys(),
    )
    m.fs.precipitation.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.precipitation.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.precipitation.initialize()
    m.fs.precipitation.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.precipitation.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mass_precipitate[
                "Calcite"
            ].value,
            1e-5,
        )
        == 4.2935e-05
    )
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mass_precipitate[
                "Gypsum"
            ].value,
            1e-5,
        )
        == 6.9149e-30
    )
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mass_precipitate[
                "Brucite"
            ].value,
            1e-5,
        )
        == 6.9149e-30
    )


@pytest.mark.component
def test_precipitation_with_custom_options():
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
        "NaOH",
        39.9971 * pyunits.g / pyunits.mol,
        {"Na_+": 1, "H2O": 1},
        solvent=("H2O", 18.01 * pyunits.g / pyunits.mol),
        purity=0.8,
    )
    m.fs.precipitation = PrecipitationUnit(
        default_property_package=m.fs.properties,
        selected_precipitants=["Anhydrite", "Aragonite"],
        viable_reagents=viable_reagents,
        viable_precipitants=viable_precipitants,
        selected_reagents=["NaOH"],
    )
    m.fs.precipitation.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.precipitation.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.precipitation.initialize()
    m.fs.precipitation.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.precipitation.report()

    print(m.fs.precipitation.config.viable_reagents.solvents)

    mols_reagent = 0.8 / 39.9971
    mols_solvent = (1 - 0.8) / 18.01
    ratio = mols_solvent / mols_reagent
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mol_solvent["H2O"].value,
            1e-5,
        )
        == 0.00011104941
    )
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mol_reagent["NaOH"].value,
            1e-5,
        )
        == 0.00020001450105132612
    )
    assert (
        pytest.approx(
            m.fs.precipitation.config.viable_reagents.solvents["NaOH"]["solvent_ratio"],
            1e-2,
        )
        == ratio
    )
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mol_solvent["H2O"].value
            / m.fs.precipitation.precipitation_reactor.flow_mol_reagent["NaOH"].value,
            1e-5,
        )
        == m.fs.precipitation.config.viable_reagents.solvents["NaOH"]["solvent_ratio"]
    )
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mass_precipitate[
                "Anhydrite"
            ].value,
            1e-5,
        )
        == 3.6981e-26
    )
    assert (
        pytest.approx(
            m.fs.precipitation.precipitation_reactor.flow_mass_precipitate[
                "Aragonite"
            ].value,
            1e-5,
        )
        == 1.7223085987945725e-06
    )

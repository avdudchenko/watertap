import pytest
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.chemical_addition_unit import (
    ChemicalAdditionUnit,
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
def test_acid_default():
    m = build_case("USDA_brackish", True)
    m.fs.acidification = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
    )
    m.fs.acidification.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.acidification.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.acidification.initialize()
    m.fs.acidification.report()
    m.fs.acidification.report(use_default_units=True)
    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.acidification.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.acidification.chemical_reactor.pH["outlet"].value,
            1e-5,
        )
        == 7.0159
    )


@pytest.mark.component
def test_acid_without_reaktoro_default():
    m = build_case("USDA_brackish", True)
    m.fs.acidification = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        add_reaktoro_chemistry=False,
    )
    m.fs.acidification.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.acidification.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.acidification.initialize()
    m.fs.acidification.report()
    m.fs.acidification.report(use_default_units=True)
    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.acidification.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.acidification.chemical_reactor.pH["outlet"].value,
            1e-5,
        )
        == 7.07
    )


@pytest.mark.component
def test_costing():
    m = build_case("USDA_brackish", True)
    m.fs.costing = WaterTAPCosting()
    m.fs.acidification = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
    )
    m.fs.acidification.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.acidification.inlet)
    m.fs.costing.cost_process()
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)
    m.fs.costing.initialize()
    m.fs.feed.initialize()
    m.fs.acidification.initialize()
    m.fs.acidification.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.acidification.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.acidification.chemical_reactor.costing.capital_cost.value,
            1e-1,
        )
        == 42.114875938347524
    )
    assert (
        pytest.approx(
            m.fs.costing.aggregate_flow_costs["fs_acidification_reagent_HCl"].value,
            1e-1,
        )
        == 53.648
    )


@pytest.mark.component
def test_acidification_with_all_options():
    m = build_case("USDA_brackish", True)
    m.fs.acidification = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        selected_reagents=ViableReagents().keys(),
    )
    m.fs.acidification.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.acidification.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.acidification.initialize()
    m.fs.acidification.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.acidification.report()
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.acidification.chemical_reactor.pH["outlet"].value,
            1e-5,
        )
        == 6.9272470967105
    )


@pytest.mark.component
def test_acidification_with_custom_options():
    m = build_case("USDA_brackish", True)
    viable_reagents = ViableReagents()
    viable_reagents.register_reagent(
        "NaOH",
        39.9971 * pyunits.g / pyunits.mol,
        {"Na_+": 1, "H2O": 1},
        solvent=("H2O", 18.01 * pyunits.g / pyunits.mol),
        purity=0.8,
    )
    m.fs.basification = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        viable_reagents=viable_reagents,
        selected_reagents=["NaOH"],
    )
    m.fs.basification.fix_and_scale()

    m.fs.feed.outlet.connect_to(m.fs.basification.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.basification.initialize()
    m.fs.basification.report()

    assert degrees_of_freedom(m) == 0

    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.fs.basification.report()

    print(m.fs.basification.config.viable_reagents.solvents)

    mols_reagent = 0.8 / 39.9971
    mols_solvent = (1 - 0.8) / 18.01
    ratio = mols_solvent / mols_reagent
    assert (
        pytest.approx(
            m.fs.basification.chemical_reactor.flow_mol_solvent["H2O"].value,
            1e-5,
        )
        == 0.00011104941
    )
    assert (
        pytest.approx(
            m.fs.basification.chemical_reactor.flow_mol_reagent["NaOH"].value,
            1e-5,
        )
        == 0.00020001450105132612
    )
    assert (
        pytest.approx(
            m.fs.basification.config.viable_reagents.solvents["NaOH"]["solvent_ratio"],
            1e-2,
        )
        == ratio
    )
    assert degrees_of_freedom(m) == 0
    assert (
        pytest.approx(
            m.fs.basification.chemical_reactor.pH["outlet"].value,
            1e-5,
        )
        == 7.18918452782
    )

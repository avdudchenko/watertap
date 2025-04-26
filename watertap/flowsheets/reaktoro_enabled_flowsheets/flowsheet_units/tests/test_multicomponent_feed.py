import pytest
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_feed import (
    MultiCompFeed,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.water_sources.source_water_importer import (
    get_source_water_data,
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


__author__ = "Alexander Dudchenko"


def build_case(water, charge_balance_with_reaktoro=False):
    mcas_props, feed_specs = get_source_water_data(f"../../water_sources/{water}.yaml")
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        charge_balance_with_reaktoro=charge_balance_with_reaktoro,
        **feed_specs,
    )
    m.fs.feed.fix_and_scale()
    m.fs.feed.report()
    return m


@pytest.mark.component
def test_mc_feed():
    m = build_case("USDA_brackish", False)
    m.fs.feed.report(use_default_units=True)
    assert degrees_of_freedom(m) == 0
    iscale.calculate_scaling_factors(m)
    m.fs.feed.initialize()
    assert pytest.approx(m.fs.feed.feed.pH.value, 1e-1) == 7
    assert (
        pytest.approx(
            m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value, 1e-2
        )
        == 0.99663795
    )
    assert (
        pytest.approx(
            m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", "Ca_2+"].value,
            1e-2,
        )
        == 0.0002580
    )


@pytest.mark.component
def test_with_reaktoro_intialization_feed():
    test_results = {
        "USDA_brackish": {
            "pH": 7.07,
            "H2O": 0.9965786471316685,
            "Cl_-": 0.9293536926966417,
        },
        "Seawater": {
            "pH": 7.56,
            "H2O": 0.9656354752143796,
            "Cl_-": 18.97752498232188,
        },
    }
    for watertype, results in test_results.items():
        m = build_case(watertype, True)
        print(m)
        iscale.calculate_scaling_factors(m)
        assert degrees_of_freedom(m) == 0
        m.fs.feed.initialize()
        assert degrees_of_freedom(m) == 0
        assert pytest.approx(m.fs.feed.feed.pH.value, 1e-3) == results["pH"]
        assert (
            pytest.approx(
                m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value,
                1e-5,
            )
            == results["H2O"]
        )
        assert (
            pytest.approx(
                m.fs.feed.feed.properties[0].conc_mass_phase_comp["Liq", "Cl_-"].value,
                1e-3,
            )
            == results["Cl_-"]
        )

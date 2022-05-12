###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
Tests for zero-order a WWTP with secondary treatment
"""
import pytest

from io import StringIO
from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    value,
    Var,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.unit_models.zero_order import SecondaryTreatmentWWTPZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestSecondaryTreatmentWWTPZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["viruses_enteric", "tss", "foo"]}
        )

        m.fs.unit = SecondaryTreatmentWWTPZO(
            default={"property_package": m.fs.params, "database": m.db}
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("secondary_treatment_wwtp")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            if j not in data["removal_frac_mass_solute"]:
                assert v.value == data["default_removal_frac_mass_solute"]["value"]
            else:
                assert v.value == data["removal_frac_mass_solute"][j]["value"]

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert (
            model.fs.unit.energy_electric_flow_vol_inlet.value
            == data["energy_electric_flow_vol_inlet"]["value"]
        )

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(1.3e-2, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(76.9231, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(76.9231, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(76.9231, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(0.01115859, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(9.042361, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(5.170008, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(89.617057, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(0.0018414, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(488.267142, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(511.73285, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(4.344496876e-7, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(6.552, abs=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    - model.fs.unit.treated.flow_mass_comp[0, j]
                    - model.fs.unit.byproduct.flow_mass_comp[0, j]
                )
            )

    @pytest.mark.component
    def test_report(self, model):
        stream = StringIO()

        model.fs.unit.report(ostream=stream)

        output = """
====================================================================================
Unit : fs.unit                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                              : Value   : Fixed : Bounds
                  Electricity Demand :  6.5520 : False : (0, None)
               Electricity Intensity : 0.14000 :  True : (None, None)
                Solute Removal [foo] :  0.0000 :  True : (0, None)
                Solute Removal [tss] : 0.94231 :  True : (0, None)
    Solute Removal [viruses_enteric] : 0.89910 :  True : (0, None)
                      Water Recovery :  1.0000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                         Inlet   Treated  Byproduct
    Volumetric Flowrate                0.013000 0.011159  0.0018414
    Mass Concentration H2O               769.23   896.17 4.3445e-07
    Mass Concentration viruses_enteric   76.923   9.0424     488.27
    Mass Concentration tss               76.923   5.1700     511.73
    Mass Concentration foo               76.923   89.617 4.3445e-07
====================================================================================
"""

        assert output in stream.getvalue()


db = Database()
params = db._get_technology("secondary_treatment_wwtp")


class Test_SecondaryTreatmentWWTPZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(default={"solute_list": ["viruses_enteric"]})

        m.fs.unit = SecondaryTreatmentWWTPZO(
            default={"property_package": m.fs.params, "database": db}
        )

        return m

    @pytest.mark.parametrize("subtype", [params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters(
            "secondary_treatment_wwtp", subtype=subtype
        )

        model.fs.unit.load_parameters_from_database()

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_solute"][j]["value"]


@pytest.mark.parametrize("subtype", [k for k in params.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = WaterParameterBlock(default={"solute_list": ["sulfur", "toc", "tss"]})

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = SecondaryTreatmentWWTPZO(
        default={
            "property_package": m.fs.params,
            "database": m.db,
            "process_subtype": subtype,
        }
    )

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )

    assert isinstance(m.fs.costing.secondary_treatment_wwtp, Block)
    assert isinstance(m.fs.costing.secondary_treatment_wwtp.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.secondary_treatment_wwtp.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.secondary_treatment_wwtp.reference_state, Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]
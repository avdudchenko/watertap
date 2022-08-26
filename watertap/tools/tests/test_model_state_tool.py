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
import pytest
from pyomo.environ import ConcreteModel, Constraint
from watertap.tools.model_state_tool import modelStateStorage
from watertap.unit_models.pressure_changer import Pump

import watertap.property_models.seawater_prop_pack as props
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale


class TestModelStateStorage:
    @pytest.fixture(scope="class")
    def pump_model(self):
        # Note direct copy of pump unit test
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = props.SeawaterParameterBlock()
        m.fs.unit = Pump(default={"property_package": m.fs.properties})
        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_TDS = 0.035
        feed_pressure_in = 1e5
        feed_pressure_out = 5e5
        feed_temperature = 273.15 + 25
        efi_pump = 0.75

        feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
            feed_flow_mass * feed_mass_frac_TDS
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)
        m.fs.unit.efficiency_pump.fix(efi_pump)
        return m

    @pytest.mark.component
    def test_var_storage(self, pump_model):
        model = pump_model
        assert model.fs.unit.inlet.pressure[0].value == 1e5
        model.prior_state = modelStateStorage(
            model,
            restorVars=True,
            restorVarFixedState=False,
            restorConstraints=False,
            restoreScaling=False,
        )
        model.fs.unit.inlet.pressure[0] = 0
        assert model.fs.unit.inlet.pressure[0].value == 0
        model.prior_state.restore_state()
        assert model.fs.unit.inlet.pressure[0].value == 1e5

    @pytest.mark.component
    def test_var_fixed_storage(self, pump_model):
        model = pump_model
        assert model.fs.unit.inlet.pressure[0].fixed == True
        model.prior_state = modelStateStorage(
            model,
            restorVars=False,
            restorVarFixedState=True,
            restorConstraints=False,
            restoreScaling=False,
        )
        model.fs.unit.inlet.pressure[0].unfix()

        assert model.fs.unit.inlet.pressure[0].fixed == False
        model.prior_state.restore_state()
        assert model.fs.unit.inlet.pressure[0].fixed == True

    @pytest.mark.component
    def test_constraint_storage(self, pump_model):
        model = pump_model
        model.temp_change_constraint = Constraint(
            expr=model.fs.unit.outlet.temperature[0]
            == model.fs.unit.inlet.pressure[0] * model.fs.unit.inlet.temperature[0]
        )
        assert model.temp_change_constraint.active == True
        model.prior_state = modelStateStorage(
            model,
            restorVars=False,
            restorVarFixedState=False,
            restorConstraints=True,
            restoreScaling=False,
        )
        model.temp_change_constraint.deactivate()
        assert model.temp_change_constraint.active == False
        model.prior_state.restore_state()
        assert model.temp_change_constraint.active == True

    @pytest.mark.component
    def test_scale_storage(self, pump_model):
        model = pump_model
        iscale.set_scaling_factor(model.fs.unit.inlet.pressure[0], 1.11)
        assert iscale.get_scaling_factor(model.fs.unit.inlet.pressure[0]) == 1.11
        model.prior_state = modelStateStorage(
            model,
            restorVars=False,
            restorVarFixedState=False,
            restorConstraints=False,
            restoreScaling=True,
        )

        iscale.set_scaling_factor(model.fs.unit.inlet.pressure[0], 1)
        assert iscale.get_scaling_factor(model.fs.unit.inlet.pressure[0]) == 1
        model.prior_state.restore_state()
        assert iscale.get_scaling_factor(model.fs.unit.inlet.pressure[0]) == 1.11

    @pytest.mark.component
    def test_separate_model(self, pump_model):
        model = pump_model
        model.temp_change_constraint = Constraint(
            expr=model.fs.unit.outlet.temperature[0]
            == model.fs.unit.inlet.pressure[0] * model.fs.unit.inlet.temperature[0]
        )
        assert model.temp_change_constraint.active == True
        cloned_model = model.clone()

        iscale.set_scaling_factor(model.fs.unit.inlet.pressure[0], 1.11)
        assert iscale.get_scaling_factor(model.fs.unit.inlet.pressure[0]) == 1.11

        iscale.set_scaling_factor(cloned_model.fs.unit.inlet.pressure[0], 1)
        assert iscale.get_scaling_factor(cloned_model.fs.unit.inlet.pressure[0]) == 1
        cloned_model.fs.unit.inlet.pressure[0] = 0
        assert cloned_model.fs.unit.inlet.pressure[0].value == 0
        cloned_model.fs.unit.inlet.pressure[0].unfix()
        assert cloned_model.fs.unit.inlet.pressure[0].fixed == False
        cloned_model.temp_change_constraint.deactivate()
        assert cloned_model.temp_change_constraint.active == False

        model.prior_state = modelStateStorage(model)
        model.prior_state.restore_state(cloned_model)

        assert cloned_model.fs.unit.inlet.pressure[0].value == 1e5
        assert cloned_model.temp_change_constraint.active == True
        assert cloned_model.fs.unit.inlet.pressure[0].fixed == True
        assert iscale.get_scaling_factor(cloned_model.fs.unit.inlet.pressure[0]) == 1.11

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

from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.environ import (
    Var,
    Constraint,
)
import idaes.core.util.scaling as iscale
from idaes.core.util import to_json, from_json


class modelStateStorage:
    """Basic class for storing model vars, constraints, and scaling values and restoring them on the fly"""

    def __init__(
        self,
        model,
        restorVars=True,
        restorVarFixedState=True,
        restorConstraints=True,
        restoreScaling=True,
        ignoreList=[],
        **kwargs
    ):
        """
        Args:
            model: the pyomo/idaes model block which to store
            restorVars: flag to store and restore variable values and fixed/unfixed state
            restorConstraints: flag to store and restore constraints states
            restoreScaling: flag to store and restore scaling of variables and constraints
        Returns:
            N/A
        """
        self.model = model
        self.restorVars = restorVars
        self.restorVarFixedState = restorVarFixedState
        self.restorConstraints = restorConstraints
        self.restoreScaling = restoreScaling
        self.ignoreList = ignoreList
        self.store_state()

    def store_state(self):
        if self.restorVars:
            self._store_vars()
        if self.restorVarFixedState:
            self._store_var_fixed_state()
        if self.restorConstraints:
            self._store_constraint_state()
        if self.restoreScaling:
            self._store_scaling_state()

    def restore_state(self, alternative_model=None):
        """Restore state for inherited model, or alternative block with same struture
        Args:
            alternative_model (optional): if passed use this model block instead, must be same in structure as original model
        Returns:
            N/A
        """

        if alternative_model is None:
            model = self.model
        else:
            model = alternative_model

        if self.restorVars:
            self._restore_vars(model)
        if self.restorVarFixedState:
            self._restore_var_fixed_state(model)
        if self.restorConstraints:
            self._restore_constraint_state(model)
        if self.restoreScaling:
            self._restore_scaling_state(model)

    def _ignore_check(self, v):
        # print(v, self.ignoreList, str(v) in self.ignoreList)
        if str(v) in self.ignoreList:
            print("ignored ", v)
            return False
        else:
            return True

    def _store_vars(self):
        """Stores model variable states"""
        self.variableStates = ComponentMap()
        for v in self.model.component_data_objects(Var):
            self.variableStates[v] = (v.value, v.lb, v.ub)

    def _store_var_fixed_state(self):
        """Stores model variable fixed states"""
        self.variableFixedStates = ComponentMap()
        for v in self.model.component_data_objects(Var):
            self.variableFixedStates[v] = v.fixed

    def _store_constraint_state(self):
        """Stores model constraint states"""
        self.constraintStates = ComponentMap()
        for v in self.model.component_data_objects(Constraint):
            self.constraintStates[v] = v.active

    def _store_scaling_state(self):
        """Stores model scaling states"""
        self.variableScalingState = ComponentMap()
        self.constrainteScalingState = ComponentMap()
        for v in self.model.component_data_objects(Var):
            self.variableScalingState[v] = iscale.get_scaling_factor(v)
        for v in self.model.component_data_objects(Constraint):
            self.constrainteScalingState[
                v
            ] = iscale.get_constraint_transform_applied_scaling_factor(v)

    def _restore_vars(self, model):
        """Stores model variable states"""
        for v, val in self.variableStates.items():
            if self._ignore_check(v):
                model.find_component(v).set_value(val[0])
                # print(v, val)
                model.find_component(v).setlb(val[1])
                model.find_component(v).setub(val[2])

    def _restore_var_fixed_state(self, model):
        """Stores model variable fixed states"""
        for v, fixed in self.variableFixedStates.items():
            if self._ignore_check(v):
                if fixed:
                    model.find_component(v).fix()
                else:
                    model.find_component(v).unfix()

    def _restore_constraint_state(self, model):
        """Stores model constraint states"""
        for v, active in self.constraintStates.items():
            if self._ignore_check(v):
                if active:
                    model.find_component(v).activate()
                else:
                    model.find_component(v).deactivate()

    def _restore_scaling_state(self, model):
        """Stores model scaling states"""
        for v, scale in self.variableScalingState.items():
            if self._ignore_check(v):
                if scale is None:
                    iscale.unset_scaling_factor(model.find_component(v))
                else:
                    iscale.set_scaling_factor(model.find_component(v), scale)
        for v, active in self.constrainteScalingState.items():
            if self._ignore_check(v):
                if scale is None:
                    iscale.constraint_scaling_transform_undo(model.find_component(v))
                else:
                    iscale.constraint_scaling_transform(model.find_component(v), scale)
        iscale.calculate_scaling_factors(model)

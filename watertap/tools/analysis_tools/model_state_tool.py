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
from pyomo.environ import Var, Constraint, Param, Expression
import idaes.core.util.scaling as iscale

import copy
from numpy.random import default_rng


class modelStateStorage:
    """Basic class for storing model vars, constraints, and scaling values and restoring them on the fly"""

    def __init__(
        self,
        model,
        restore_dict=None,
        restorVars=True,
        restorVarBounds=True,
        restorVarFixedState=True,
        restorConstraints=True,
        restoreScaling=True,
        restoreParams=True,
        evalExpressions=True,
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
        self.restorVarBounds = restorVarBounds
        self.restoreParams = restoreParams
        self.evalExpressions = evalExpressions
        self.ignoreList = ignoreList
        if restore_dict is not None:
            self.restore_from_dict(restore_dict)

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
        if self.restoreParams:
            self._store_params()

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
        if self.restoreParams:
            self._restore_params(model)
        if self.evalExpressions:
            self._eval_expressions(model)

    def get_dict_state(self):
        return_dict = {}
        if self.restorVars:
            return_dict["variableStatesDict"] = self.variableStatesDict
        if self.restorVarFixedState:
            return_dict["variableFixedStatesDict"] = self.variableFixedStatesDict
        if self.restorConstraints:
            return_dict["constraintStatesDict"] = self.constraintStatesDict
        if self.restoreParams:
            return_dict["paramStatesDict"] = self.paramStatesDict
        if self.restoreScaling:
            return_dict["variableScalingStateDict"] = self.variableScalingStateDict
            return_dict[
                "constrainteScalingStateDict"
            ] = self.constrainteScalingStateDict
        return copy.deepcopy(return_dict)

    def restore_from_dict(self, restore_dict):
        if "variableStatesDict" in restore_dict:
            self.variableStatesDict = restore_dict["variableStatesDict"]
            self.restorVars = True
        if "variableFixedStatesDict" in restore_dict:
            self.variableFixedStatesDict = restore_dict["variableFixedStatesDict"]
            self.restorVarFixedState = True
        if "constraintStatesDict" in restore_dict:
            self.constraintStatesDict = restore_dict["constraintStatesDict"]
            self.restorConstraints = True
        if "paramStatesDict" in restore_dict:
            self.paramStatesDict = restore_dict["paramStatesDict"]
            self.restorParams = True
        if "variableScalingStateDict" in restore_dict:
            self.variableScalingStateDict = restore_dict["variableScalingStateDict"]
            self.constrainteScalingStateDict = restore_dict[
                "constrainteScalingStateDict"
            ]
            self.restoreScaling = True
        self.restore_state()

    def _ignore_check(self, v):
        # print(v, self.ignoreList, str(v) in self.ignoreList)
        for val in self.ignoreList:
            if val in str(v):
                # if str(v) in self.ignoreList:
                print("State model tool ignored ", val, "==", v)
                return False
                break
        else:
            return True

    def _store_vars(self):
        """Stores model variable states"""
        self.variableStates = ComponentMap()
        self.variableStatesDict = {}
        for v in self.model.component_data_objects(Var):
            self.variableStates[v] = (v.value, v.lb, v.ub)
            self.variableStatesDict[str(v)] = {"value": v.value, "lb": v.lb, "ub": v.ub}

    def _store_params(self):
        """Stores model param states"""
        self.paramStates = ComponentMap()
        self.paramStatesDict = {}
        for v in self.model.component_data_objects(Param):
            self.paramStates[v] = v.value
            self.paramStatesDict[str(v)] = {"value": v.value}

    def _store_var_fixed_state(self):
        """Stores model variable fixed states"""
        self.variableFixedStates = ComponentMap()
        self.variableFixedStatesDict = {}
        for v in self.model.component_data_objects(Var):
            # p#rint(str(v))
            self.variableFixedStates[v] = v.fixed
            self.variableFixedStatesDict[str(v)] = v.fixed

    def _store_constraint_state(self):
        """Stores model constraint states"""
        self.constraintStates = ComponentMap()
        self.constraintStatesDict = {}
        for v in self.model.component_data_objects(Constraint):
            self.constraintStates[v] = v.active
            self.constraintStatesDict[str(v)] = v.active

    def _store_scaling_state(self):
        """Stores model scaling states"""
        self.variableScalingState = ComponentMap()
        self.variableScalingStateDict = {}

        self.constrainteScalingState = ComponentMap()
        self.constrainteScalingStateDict = {}
        for v in self.model.component_data_objects(Var):
            scale = iscale.get_scaling_factor(v)

            self.variableScalingState[v] = scale
            self.variableScalingStateDict[str(v)] = scale
        for v in self.model.component_data_objects(Constraint):
            scalin_value = iscale.get_constraint_transform_applied_scaling_factor(v)
            self.constrainteScalingState[v] = scalin_value
            self.constrainteScalingStateDict[str(v)] = scalin_value

    def _restore_vars(self, model):
        """Restores model variable states"""
        for v, val in self.variableStatesDict.items():
            if self._ignore_check(v):
                if self.restorVarBounds:
                    model.find_component(v).setlb(val["lb"])
                    model.find_component(v).setub(val["ub"])
                model.find_component(v).set_value(val["value"])
                # print(v, val)
                # model.find_component(v).setlb(val["lb"])
                # model.find_component(v).setub(val["ub"])

    def _restore_params(self, model):
        """Restores model param states"""
        for v, val in self.paramStatesDict.items():
            if self._ignore_check(v):
                # print(v, val)
                model.find_component(v).set_value(val["value"])

    def _restore_var_fixed_state(self, model):
        """Restores model variable fixed states"""
        for v, fixed in self.variableFixedStatesDict.items():
            if self._ignore_check(v):
                if fixed:
                    model.find_component(v).fix()
                else:
                    model.find_component(v).unfix()

    def _restore_constraint_state(self, model):
        """Restores model constraint states"""
        for v, active in self.constraintStatesDict.items():
            if self._ignore_check(v):
                try:
                    if active:
                        model.find_component(v).activate()
                    else:
                        model.find_component(v).deactivate()
                except:
                    print("faile seeting state", v)

    def _restore_scaling_state(self, model):
        """Restores model scaling states"""
        for v, scale in self.variableScalingStateDict.items():
            if self._ignore_check(v):
                if scale is None:
                    iscale.unset_scaling_factor(model.find_component(v))
                else:
                    iscale.set_scaling_factor(model.find_component(v), scale)
        for v, scale in self.constrainteScalingStateDict.items():
            if self._ignore_check(v):
                if scale is None:
                    iscale.constraint_scaling_transform_undo(model.find_component(v))
                else:
                    iscale.constraint_scaling_transform(model.find_component(v), scale)

    def _eval_expressions(self, model):
        """Evals model expressions"""
        for v in model.component_data_objects(Expression):
            try:
                value(v)
            except:
                pass

    def find_changed_vars(self, model):
        """Finds changed vars"""
        self.changedVariables = ComponentMap()
        for v, val in self.variableStatesDict.items():
            cur_val = model.find_component(v).value
            if cur_val != val["value"] and model.find_component(v).fixed:
                self.changedVariables[v] = (val["value"], cur_val)
        for v, val in self.paramStatesDict.items():
            cur_val = model.find_component(v).value
            if cur_val != val["value"]:
                self.changedVariables[v] = (val["value"], cur_val)
        return self.changedVariables

    def find_changed_scailing(self, model):
        """Finds changed vars"""
        self.changedVarScales = {}
        for v, scale in self.variableScalingStateDict.items():
            cur_scale = iscale.get_scaling_factor(model.find_component(v))
            if scale != cur_scale:
                self.changedVariables[v] = {scale, cur_scale}
        self.changedConstraintScales = {}
        for v, scale in self.constrainteScalingStateDict.items():
            cur_scale = iscale.get_constraint_transform_applied_scaling_factor(
                model.find_component(v)
            )
            if scale != cur_scale:
                self.changedConstraintScales[v] = {scale, cur_scale}
        return self.changedVarScales, self.changedConstraintScales

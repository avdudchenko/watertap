import numpy as np
from pyomo.common.collections import ComponentSet, ComponentMap
from watertap.tools.analysis_tools.model_state_tool import modelStateStorage
from pyomo.environ import (
    assert_optimal_termination,
)
from idaes.core.util import to_json, from_json
import time


def updateLocalState(model, state=None, returnState=False):
    if returnState is True or state is None:
        return to_json(model, return_json_string=True)
    else:
        from_json(model, s=state)


def step_optimize(
    model,
    solver,
    solve_func,
    init_state,
    steps=10,
    re_steps=3,
    try_final_first=True,
    step_back_threshold=None,
    perturb_unfixed_vars=None,
    test_succesfull_solve=True,
):
    solve_start = time.time()
    changed_vars = init_state.find_changed_vars(model)
    step_vars = gen_steps(changed_vars, steps)
    last_step = None
    model.last_solved_state = modelStateStorage(model, restorVarBounds=False)

    for k in range(re_steps):
        if last_step != None:
            step_vars = gen_steps(
                last_step,
                steps=steps * (k + 1),
                step_back_threshold=step_back_threshold,
            )
        results, success, last_step = try_stepping(
            model,
            solver,
            solve_func,
            steps * (k + 1),
            step_vars,
            try_final_first,
            perturb_unfixed_vars,
        )
        if success:
            break
    if test_succesfull_solve:
        assert_optimal_termination(results)
    print(
        "----step-solve-took: {} seconds------------".format(time.time() - solve_start)
    )
    return results


def fix_value(param, val):
    if param.is_variable_type():
        # Fix the single value to values[k]
        param.fix(val)

    elif param.is_parameter_type():
        # Fix the single value to values[k]
        param.set_value(val)


def try_stepping(
    model, solver, solve_func, steps, step_vars, try_final_first, perturb_unfixed_vars
):
    solved_successfull = False
    results = None
    last_step = None
    step_diff = ComponentMap()
    local_state = updateLocalState(model)
    if perturb_unfixed_vars is not None:
        model.last_solved_state.perturb_unfixed_vars(perturb_unfixed_vars)

    for v, vals in step_vars.items():
        step_diff[v] = (vals[0], vals[-1])
    if try_final_first:
        for v, vals in step_vars.items():
            # print(v, vals[i], vals[-1])
            fix_value(model.find_component(v), vals[-1])
        try:
            results = solve_func(model, solver)
            assert_optimal_termination(results)
            print("solved on first try!!!")
            # model.last_solved_state.store_state()
            # local_state = updateLocalState(model)
            solved_successfull = True
        except:
            updateLocalState(model, local_state)
            model.last_solved_state.store_state()
            if perturb_unfixed_vars is not None:
                model.last_solved_state.perturb_unfixed_vars(perturb_unfixed_vars)

    if solved_successfull == False:
        for i in range(steps):
            print("taking step ", i, steps)
            for v, vals in step_vars.items():
                print("taking step {} {} {}".format(i, v, vals[i]))
                # print(v, vals[i], vals[-1])
                # model.find_component(v).fix(vals[i])
                fix_value(model.find_component(v), vals[i])
            try:
                results = solve_func(model, solver)
                assert_optimal_termination(results)
                for v, vals in step_vars.items():
                    step_diff[v] = (vals[i], vals[-1])
                local_state = updateLocalState(model)
                model.last_solved_state.store_state()
                solved_successfull = True

            except:
                updateLocalState(model, local_state)
                model.last_solved_state.store_state()
                solved_successfull = False
                results = None
                last_step = step_diff
                if perturb_unfixed_vars is not None:
                    model.last_solved_state.perturb_unfixed_vars(perturb_unfixed_vars)
                break
            print("solved succesful: ", solved_successfull)

    return results, solved_successfull, last_step


def gen_steps(changed_vars, steps=5, step_back_threshold=None):
    step_vars = ComponentMap()
    for v, vals in changed_vars.items():
        # print("Generated step", v, vals, steps)
        # if vals[0] > vals[1]:
        #     val_0 = vals[0] - vals[1] * 0.01
        # else:
        #     val_0 = vals[0] + vals[1] * 0.01
        # if vals[1] > vals[0]:
        #     val_ini = vals[0] * 1.01
        # else:
        val_ini = vals[0]  # * 0.99
        # if step_back_threshold is not None:
        #     if abs(vals[1] - vals[0]) / abs(vals[0]) < step_back_threshold:

        #         if vals[1] > vals[0]:
        #             val_ini = vals[0] - vals[1] * step_back_threshold
        #         if vals[1] < vals[0]:
        #             val_ini = vals[0] + vals[1] * step_back_threshold
        #         print(
        #             "threshold reached stepping expanding gap!{} : {} expanded to {}, moving to {}".format(
        #                 v, vals[0], val_ini, vals[1]
        #             )
        #         )
        # print("Step generated {} : {} --> {}".format(v, val_ini, vals[1]))
        step_vars[v] = list(np.linspace(val_ini, vals[1], steps))
    return step_vars

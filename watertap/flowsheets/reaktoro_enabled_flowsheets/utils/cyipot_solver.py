from watertap_solvers import get_solver


def get_cyipopt_solver(max_iters=3000, limited_memory=True):
    """Get a cyipopt solver with custom options."""
    solver_name = "cyipopt-watertap"
    solver = get_solver(solver_name)
    solver.options["acceptable_dual_inf_tol"] = 1.0e-8
    solver.options["max_iter"] = max_iters
    solver.options["linear_solver"] = "ma27"
    # solver.options["print_level"] = 10
    # solver.options["recalc_y"] = "yes"
    # solver.options["recalc_y_feas_tol"] = 1e-2
    # solver.options["ma27_pivtol"] = 1e-6
    # solver.options["ma27_pivtolmax"] = 0.1
    # solver.options["min_refinement_steps"] = 2
    if limited_memory:
        solver.options["hessian_approximation"] = "limited-memory"
    return solver

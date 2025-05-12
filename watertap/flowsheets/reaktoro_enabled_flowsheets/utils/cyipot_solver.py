from watertap_solvers import get_solver


def get_cyipopt_solver(max_iters=300):
    """Get a cyipopt solver with custom options."""
    solver_name = "cyipopt-watertap"
    solver = get_solver(solver_name)
    solver.options["acceptable_dual_inf_tol"] = 1.0e-8
    solver.options["max_iter"] = max_iters
    solver.options["linear_solver"] = "ma27"

    return solver

from watertap_solvers import get_solver


def get_cyipopt_solver():
    solver_name = "cyipopt-watertap"
    solver = get_solver(solver_name)
    solver.options["acceptable_dual_inf_tol"] = 1.0e-8
    solver.options["max_iter"] = 300
    solver.options["linear_solver"] = "ma27"

    return solver

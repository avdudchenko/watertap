#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
from watertap.examples.flowsheets.nf_dspmde.nf import main


@pytest.fixture
def model():
    """
    # Example Usage:
    # Set the number of trials
    nn = 50
    results = np.zeros((nn, 2))

    for k in range(nn):
        # Set a random value for the parameter a
        m.fs.a.set_value(np.random.fs.rand())

        # Attempt to solve the model (infeasible if m.fs.a > success_prob)
        solver.solve(m)

        # Store the value of a and the solution x
        a = pyo.value(m.fs.a)
        x = pyo.value(m.fs.x)

        if np.abs(pyo.value(m.fs.err)) < 1e-6:
            results[k, :] = [a, x]
        else:
            results[k, :] = [a, np.nan]
    """

    # Initialize a model and solver
    m = pyo.ConcreteModel()

    m.fs = pyo.Block()

    # Declare decision variable and param
    m.fs.x = pyo.Var()
    m.fs.a = pyo.Param(mutable=True)
    m.fs.success_prob = pyo.Param(initialize=0.5, mutable=True)

    # Define expressions and constraints:
    # Numbers must sum to success_prob and x must be positive
    m.fs.err = pyo.Expression(expr=m.fs.a + m.fs.x - m.fs.success_prob)
    m.fs.sum = pyo.Constraint(expr=m.fs.err == 0.0)
    m.fs.pos = pyo.Constraint(expr=m.fs.x >= 0.0)

    return m


@pytest.mark.component
def test_main(model):
    m = model

    for key, (solve, testval) in test_dict.items():
        assert round(solve, 4) == round(testval, 4)

from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_feed_unit import (
    MultiCompFeed,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.water_sources.source_water_importer import (
    get_source_water_data,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_product_unit import (
    MultiCompProduct,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_ph_mixer_unit import (
    MixerPhUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_erd_unit import (
    MultiCompERDUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_pump_unit import (
    MultiCompPumpUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.multi_comp_ro_unit import (
    MultiCompROUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.precipitation_unit import (
    PrecipitationUnit,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.flowsheet_units.chemical_addition_unit import (
    ChemicalAdditionUnit,
)
from pyomo.environ import (
    TransformationFactory,
    units as pyunits,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.costing import WaterTAPCosting
from pyomo.environ import ConcreteModel, Var, Reals, Constraint, Objective
from idaes.core import (
    FlowsheetBlock,
)
import idaes.core.util.scaling as iscale

from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.report_util import (
    build_report_table,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.util.model_diagnostics.infeasible import *
import os
import tempfile
import re


def main():
    feed_water = "../water_sources/USDA_brackish.yaml"
    feed_water = "../water_sources/sample_500_hardness.yaml"
    feed_water = "../water_sources/sample_1500_hardness.yaml"
    # feed_water = "../water_sources/Seawater.yaml"
    m = build_model(
        feed_water,
        hpro=False,
        rkt_hessian_type="LBFGS",
        bfgs_initialization_type="GaussNewton",
    )
    initialize(m)
    if False:  # "USDA" in feed_water:
        m.fs.water_recovery.fix(60 / 100)
        solve_model(m)
        m.fs.water_recovery.fix(65 / 100)
        solve_model(m)
        start = 70
    else:
        start = 50
    for r in range(start, 91, 1):
        m.fs.water_recovery.fix(r / 100)
        print(f"\n\n------------Solving for water recovery: {r}%------------")
        solve_model(m, tee=True)
        report_all_units(m)
        print("------vars_close_to_bound-tests---------")
        print_variables_close_to_bounds(m)
        print("------constraints_close_to_bound-tests---------")

        print_constraints_close_to_bounds(m)


def enable_multi_process_reaktoro(
    m, rkt_hessian_type="LBFGS", bfgs_initialization_type="GaussNewton"
):
    """Enables use of parallel solves for reaktoro blocks,
    in RO mode there will be 3 reaktoro blocks
        1 for Softening,
        1 for Acidification,
        1 for RO
        requiring 3 logical cores
    in HPRO mode there will be 4 reaktoro blocks
        1 for Softening,
        1 for Acidification,
        1 for RO,
        1 for HPRO
        requiring 4 logical cores
    """
    from reaktoro_pse.parallel_tools.reaktoro_block_manager import (
        ReaktoroBlockManager,
    )

    rkt_options = {}
    opt = {
        "hessian_type": rkt_hessian_type,
        "bfgs_initialization_type": bfgs_initialization_type,
    }
    if bfgs_initialization_type == "constant":
        opt["bfgs_init_const_hessian_value"] = 1e-16

    m.reaktoro_manager = ReaktoroBlockManager(
        hessian_options=opt,
    )

    rkt_options["reaktoro_block_manager"] = m.reaktoro_manager
    return rkt_options


def build_model(
    water_case,
    multi_process_reaktoro=True,
    hpro=False,
    rkt_hessian_type="BFGS",
    bfgs_initialization_type="GaussNewton",
    # rkt_scaling_type="variable_oi_scaling_square_sum",
):
    """Builds the flowsheet model for the softening-acidification-RO process.
    Args:
        water_case (str): location of water source data file with yaml format.
        multi_process_reaktoro (bool): If True, enables parallel processing for reaktoro blocks.
        hpro (bool): If True, includes high-pressure RO unit in the flowsheet.
        rkt_hessian_type (str): Type of Hessian to use in reaktoro blocks.
            Options are:
                - limited-memory: Uses limited memory BFGS Hessian approximation within ipopt solver.
                - ZeroHessian - no hessian
                - GaussNewton - Naive Gauss-Newton Hessian approximation (Jacobian^T * Jacobian)
                - LBFGS - Limited Memory Broyden-Fletcher-Goldfarb-Shanno
                - BFGS - Broyden-Fletcher-Goldfarb-Shanno
                - CBFGS - conditional BFGS
                - BFGS_mod - modified BFGS
                - BFGS_damp - damped BFGS
                - BFGS_ipopt - BFGS with ipopt update step
    """

    mcas_props, feed_specs = get_source_water_data(water_case)
    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m = ConcreteModel()
    m.water_case = water_case
    if rkt_hessian_type == "limited-memory":
        rkt_hessian_type = "ZeroHessian"
        m.solver_limited_memory = True
    else:
        m.solver_limited_memory = False
    m.solver_limited_memory_scalar = bfgs_initialization_type
    rkt_options = {
        "hessian_options": {
            "hessian_type": rkt_hessian_type,
            "bfgs_initialization_type": bfgs_initialization_type,
        }
        # "jacobian_options": {
        #     "scaling_type": rkt_scaling_type,
        # },
    }
    if bfgs_initialization_type == "constant":
        rkt_options["hessian_options"]["bfgs_init_const_hessian_value"] = 1e-16
    if multi_process_reaktoro:
        rkt_options = enable_multi_process_reaktoro(
            m, rkt_hessian_type, bfgs_initialization_type
        )
        # rkt_options["jacobian_options"] = {
        #     "scaling_type": rkt_scaling_type,
        # }

    m.fs = FlowsheetBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.ro_properties = SeawaterParameterBlock()
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        reconcile_using_reaktoro=True,
        **feed_specs,
    )
    m.fs.softening_unit = PrecipitationUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        selected_precipitants=[
            "Calcite",
            "Brucite",
        ],
        selected_reagents=["Na2CO3", "CaO"],
        add_alkalinity=True,
        reaktoro_options=rkt_options,
    )

    m.fs.acidification_unit = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        selected_reagents=["HCl", "H2SO4"],
        reaktoro_options=rkt_options,
    )
    if "Seawater" in water_case:
        ovp = 1.5
    else:
        ovp = 6.5
    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        initialization_pressure="osmotic_pressure",
        osmotic_over_pressure=ovp,
        maximum_pressure=85 * pyunits.bar,
    )

    m.fs.ro_unit = MultiCompROUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        ro_property_package=m.fs.ro_properties,
        selected_scalants={"Calcite": 1, "Gypsum": 1},
        reaktoro_options=rkt_options,
    )

    m.fs.erd_unit = MultiCompERDUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
    )

    if hpro:
        m.fs.hp_pump_unit = MultiCompPumpUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            initialization_pressure="osmotic_pressure",
            maximum_pressure=400 * pyunits.bar,
        )
        m.fs.hpro_unit = MultiCompROUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            ro_property_package=m.fs.ro_properties,
            selected_scalants={"Calcite": 1, "Gypsum": 1},
            reaktoro_options=rkt_options,
            default_costing_package_kwargs={
                "costing_method_arguments": {"ro_type": "high_pressure"}
            },
        )
        m.fs.product_mixer = MixerPhUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            inlet_ports=["ro_inlet", "hpro_inlet"],
            add_reaktoro_chemistry=False,
        )

    m.fs.product = MultiCompProduct(
        default_property_package=m.fs.properties,
    )
    m.fs.brine = MultiCompProduct(
        default_property_package=m.fs.properties,
    )
    if multi_process_reaktoro:
        m.reaktoro_manager.build_reaktoro_blocks()
    # to simplify initialization and fixing model state
    m.flowsheet_unit_order = []
    m.flowsheet_unit_order.append(m.fs.feed)
    m.flowsheet_unit_order.append(m.fs.softening_unit)
    m.flowsheet_unit_order.append(m.fs.acidification_unit)
    m.flowsheet_unit_order.append(m.fs.pump_unit)
    m.flowsheet_unit_order.append(m.fs.ro_unit)

    if hpro:
        m.flowsheet_unit_order.append(m.fs.hp_pump_unit)
        m.flowsheet_unit_order.append(m.fs.hpro_unit)
        m.flowsheet_unit_order.append(m.fs.product_mixer)

    m.flowsheet_unit_order.append(m.fs.erd_unit)
    m.flowsheet_unit_order.append(m.fs.product)
    m.flowsheet_unit_order.append(m.fs.brine)

    m.fs.reaktoro_blocks = []
    m.fs.reaktoro_blocks.append(m.fs.softening_unit.precipitation_block)
    m.fs.reaktoro_blocks.append(m.fs.acidification_unit.chemistry_block)
    m.fs.reaktoro_blocks.append(m.fs.ro_unit.scaling_block)
    if hpro:
        m.fs.reaktoro_blocks.append(m.fs.hpro_unit.scaling_block)

    # build all connections
    m.fs.feed.outlet.connect_to(m.fs.softening_unit.inlet)
    m.fs.softening_unit.outlet.connect_to(m.fs.acidification_unit.inlet)
    m.fs.acidification_unit.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ro_unit.feed)

    if hpro:
        m.fs.ro_unit.retentate.connect_to(m.fs.hp_pump_unit.inlet)
        m.fs.hp_pump_unit.outlet.connect_to(m.fs.hpro_unit.feed)
        m.fs.hpro_unit.retentate.connect_to(m.fs.erd_unit.inlet)

        m.fs.ro_unit.product.connect_to(m.fs.product_mixer.ro_inlet)
        m.fs.hpro_unit.product.connect_to(m.fs.product_mixer.hpro_inlet)
        m.fs.product_mixer.outlet.connect_to(m.fs.product.inlet)
    else:
        m.fs.ro_unit.retentate.connect_to(m.fs.erd_unit.inlet)
        m.fs.ro_unit.product.connect_to(m.fs.product.inlet)
    m.fs.erd_unit.outlet.connect_to(m.fs.brine.inlet)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.product.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.product.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.product.product.properties[0].flow_vol
    )
    # brackish  water cost is low, so lets scale it up
    m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    TransformationFactory("network.expand_arcs").apply_to(m)
    add_global_constraints(m)
    fix_and_scale(m)
    add_perfoance_tracking_vars(m)
    return m


def add_global_constraints(m):
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0.0, 1),
        domain=Reals,
        units=pyunits.dimensionless,
        doc="System Water Recovery",
    )
    m.fs.eq_water_recovery = Constraint(
        expr=sum(m.fs.feed.feed.properties[0].flow_mass_phase_comp["Liq", :])
        * m.fs.water_recovery
        == m.fs.product.product.properties[0].flow_mass_phase_comp["Liq", "H2O"]
    )
    iscale.set_scaling_factor(m.fs.eq_water_recovery, 1)
    iscale.constraint_scaling_transform(m.fs.eq_water_recovery, 1)


def add_perfoance_tracking_vars(m):
    m.fs.ipopt_iterations = Var(
        [
            "Number of iterations",
            "Objective function evaluations",
            "Objective gradient evaluations",
            "Equality constraint evaluations",
            "Inequality constraint evaluations",
            "Equality constraint Jacobian evaluations",
            "Inequality constraint Jacobian evaluations",
            "Lagrangian Hessian evaluations",
        ],
        initialize=0,
        domain=Reals,
        doc="Number of ipopt iterations",
    )
    m.fs.ipopt_iterations.fix()
    m.fs.scaled_ipopt_result = Var(
        [
            "dual_infeasibility",
            "constraint_violation",
            "complementarity_error",
            "variable_bound_violation",
            "overall_nlp_error",
        ],
        initialize=0,
        domain=Reals,
        doc="Scaled dual infeasibility from ipopt result",
    )
    m.fs.unscaled_ipopt_result = Var(
        [
            "dual_infeasibility",
            "constraint_violation",
            "complementarity_error",
            "variable_bound_violation",
            "overall_nlp_error",
        ],
        initialize=0,
        domain=Reals,
        doc="unscaled dual infeasibility from ipopt result",
    )
    m.fs.scaled_ipopt_result.fix()


def report_global_state(m):
    data_dict = {"Global results": {}}
    data_dict["Global results"]["DOfs"] = int(degrees_of_freedom(m))
    data_dict["Global results"]["Water recovery"] = m.fs.water_recovery
    data_dict["Global results"]["LCOW"] = m.fs.costing.LCOW

    build_report_table("Global results", data_dict)


def fix_and_scale(m):
    for unit in m.flowsheet_unit_order:
        unit.fix_and_scale()
    # limiiting maximum pH during optimization for
    # stability - in all example waters, pH in softening does not really operate
    # above 10, but ability of softening unit to go above 10 can cuase instability
    # when alkalinitry drops bellow <20-25 ppm.
    max_ph = 11
    iscale.calculate_scaling_factors(m)
    m.fs.softening_unit.precipitation_reactor.alkalinity.setlb(10)
    m.fs.ro_unit.ro_feed.pH.setlb(6)
    m.fs.ro_unit.ro_feed.pH.setub(max_ph)
    if m.fs.find_component("hp_pump_unit") is not None:
        m.fs.hpro_unit.ro_feed.pH.setlb(6)
        m.fs.hpro_unit.ro_feed.pH.setub(max_ph)
    m.fs.softening_unit.precipitation_reactor.pH["outlet"].setlb(6)
    m.fs.softening_unit.precipitation_reactor.pH["outlet"].setub(max_ph)
    assert degrees_of_freedom(m) == 0


def initialize(m, **kwargs):
    for unit in m.flowsheet_unit_order:
        unit.initialize()
    m.fs.costing.initialize()
    # report_all_units(m)
    solve_model(m)
    set_optimization(m)

    if m.fs.water_recovery.value < 0.5:
        m.fs.water_recovery.fix()
        solve_model(m)
        m.fs.water_recovery.fix(0.5)
    else:
        m.fs.water_recovery.fix()
    solve_model(m)
    print("--------------Initialization complete--------")


def report_all_units(m):
    for unit in m.flowsheet_unit_order:
        unit.report()
    report_global_state(m)


def set_optimization(m):
    for unit in m.flowsheet_unit_order:
        unit.set_optimization_operation()
    report_all_units(m)


def test_func(m, **kwargs):
    """Test function to check if the model is suitable for solving
    For Seawater case, with out HPRO maximum recoveriy is 65% while, with HPRO it is 86%
    """
    print(
        f"Testing func: {m.water_case}, {m.fs.find_component('hpro_unit')}, {m.fs.water_recovery.value}"
    )

    if "Seawater" in m.water_case:
        if (
            m.fs.find_component("hpro_unit") is None
            and m.fs.water_recovery.value > 0.64
        ):
            return False
        elif m.fs.water_recovery.value > 0.86:
            return False
    return True


def solve_model(m, tee=True, **kwargs):
    solver = get_cyipopt_watertap_solver(
        ma27=True,
        max_iter=1000,
        limited_memory=m.solver_limited_memory,
        scalar_type=m.solver_limited_memory_scalar,
    )
    print(solver.options["linear_solver"])
    if m.fs.find_component("ipopt_iterations") is not None:
        tmp = tempfile.NamedTemporaryFile(suffix=".txt", delete=False)
        tmp.close()
        print(tmp.name)
        solver.options["output_file"] = tmp.name
    else:
        tmp = None
    result = solver.solve(m, tee=tee)
    if tmp is not None:
        with open(tmp.name) as f:
            lines = f.read()
        iters_match = re.search(
            r"""
            Number\ of\ objective\ function\ evaluations\s*=\s*([-+eE0-9.]+)\s*
            Number\ of\ objective\ gradient\ evaluations\s*=\s*([-+eE0-9.]+)\s*
            Number\ of\ equality\ constraint\ evaluations\s*=\s*([-+eE0-9.]+)\s*
            Number\ of\ inequality\ constraint\ evaluations\s*=\s*([-+eE0-9.]+)\s*
            Number\ of\ equality\ constraint\ Jacobian\ evaluations\s*=\s*([-+eE0-9.]+)\s*
            Number\ of\ inequality\ constraint\ Jacobian\ evaluations\s*=\s*([-+eE0-9.]+)\s*
            Number\ of\ Lagrangian\ Hessian\ evaluations\s*=\s*([-+eE0-9.]+)\s*
            """,
            lines,
            re.DOTALL | re.VERBOSE,
        )
        os.remove(tmp.name)
        parsed_output = _parse_ipopt_output(lines)

        iters_keys = list(m.fs.ipopt_iterations.keys())
        for i, k in enumerate(iters_match.groups()):
            if k is not None:
                k = int(float(k))
            else:
                k = 0

            m.fs.ipopt_iterations[iters_keys[i]] = k
        for key in m.fs.unscaled_ipopt_result:
            m.fs.unscaled_ipopt_result[key] = parsed_output.get(key, 0)
            m.fs.scaled_ipopt_result[key] = parsed_output["final_scaled_results"].get(
                key, 0
            )
        m.fs.ipopt_iterations["Number of iterations"] = int(parsed_output["iters"])
        # m.fs.ipopt_iterations.display()
        # m.fs.scaled_ipopt_result.display()
        # m.fs.unscaled_ipopt_result.display()
    # report_all_units(m)
    assert_optimal_termination(result)
    return result


def _parse_ipopt_output(output):
    _ALPHA_PR_CHARS = set("fFhHkKnNRwstTr")

    parsed_data = {}
    import io

    # Convert output to a string so we can parse it
    if isinstance(output, io.StringIO):
        output = output.getvalue()

    # Extract number of iterations
    iter_match = re.search(r"Number of Iterations.*:\s+(\d+)", output)
    if iter_match:
        parsed_data["iters"] = int(iter_match.group(1))
    # Gather all the iteration data
    iter_table = re.findall(r"^(?:\s*\d+.*?)$", output, re.MULTILINE)
    if iter_table:
        columns = [
            "iter",
            "objective",
            "inf_pr",
            "inf_du",
            "lg_mu",
            "d_norm",
            "lg_rg",
            "alpha_du",
            "alpha_pr",
            "ls",
        ]
        iterations = []

        for line in iter_table:
            tokens = line.strip().split()
            if len(tokens) != len(columns):
                continue
            iter_data = dict(zip(columns, tokens))

            # Extract restoration flag from 'iter'
            iter_data["restoration"] = iter_data["iter"].endswith("r")
            if iter_data["restoration"]:
                iter_data["iter"] = iter_data["iter"][:-1]

            # Separate alpha_pr into numeric part and optional tag
            iter_data["step_acceptance"] = iter_data["alpha_pr"][-1]
            if iter_data["step_acceptance"] in _ALPHA_PR_CHARS:
                iter_data["alpha_pr"] = iter_data["alpha_pr"][:-1]
            else:
                iter_data["step_acceptance"] = None

            # Attempt to cast all values to float where possible
            for key in columns:
                if iter_data[key] == "-":
                    iter_data[key] = None
                else:
                    try:
                        iter_data[key] = float(iter_data[key])
                    except (ValueError, TypeError):
                        print(
                            "Error converting Ipopt log entry to "
                            f"float:\n\t{sys.exc_info()[1]}\n\t{line}"
                        )

            # assert len(iterations) == iter_data.pop("iter"), (
            #     f"Parsed row in the iterations table\n\t{line}\ndoes not "
            #     f"match the next expected iteration number ({len(iterations)})"
            # )
            iterations.append(iter_data)

        parsed_data["iteration_log"] = iterations

    # Extract scaled and unscaled table
    scaled_unscaled_match = re.search(
        r"""
        Objective\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        Dual\ infeasibility\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        Constraint\ violation\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        (?:Variable\ bound\ violation:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*)?
        Complementarity\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        Overall\ NLP\ error\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)
        """,
        output,
        re.DOTALL | re.VERBOSE,
    )

    if scaled_unscaled_match:
        groups = scaled_unscaled_match.groups()
        all_fields = [
            "incumbent_objective",
            "dual_infeasibility",
            "constraint_violation",
            "variable_bound_violation",  # optional
            "complementarity_error",
            "overall_nlp_error",
        ]

        # Filter out None values and create final fields and values.
        # Nones occur in old-style IPOPT output (<= 3.13)
        zipped = [
            (field, scaled, unscaled)
            for field, scaled, unscaled in zip(all_fields, groups[0::2], groups[1::2])
            if scaled is not None and unscaled is not None
        ]

        scaled = {k: float(s) for k, s, _ in zipped}
        unscaled = {k: float(u) for k, _, u in zipped}

        parsed_data.update(unscaled)
        parsed_data["final_scaled_results"] = scaled

    # Newer versions of IPOPT no longer separate timing into
    # two different values. This is so we have compatibility with
    # both new and old versions
    parsed_data["cpu_seconds"] = {
        k.strip(): float(v)
        for k, v in re.findall(
            r"Total(?: CPU)? sec(?:ond)?s in ([^=]+)=\s*([0-9.]+)", output
        )
    }

    return parsed_data


if __name__ == "__main__":
    main()

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
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.cyipot_solver import (
    get_cyipopt_solver,
)

from watertap.costing import WaterTAPCosting
from pyomo.environ import ConcreteModel, Var, Reals, Constraint, Objective
from idaes.core import (
    FlowsheetBlock,
)
import idaes.core.util.scaling as iscale

from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.cyipot_solver import (
    get_cyipopt_solver,
)
from pyomo.environ import (
    assert_optimal_termination,
)
from watertap.flowsheets.reaktoro_enabled_flowsheets.utils.report_util import (
    build_report_table,
)
from idaes.core.util.model_statistics import degrees_of_freedom
import numpy as np

from watertap.core.util.model_diagnostics.infeasible import *


def test_softeining():
    m = build_lime_only("../water_sources/USDA_brackish.yaml")
    solver = get_cyipopt_solver()
    # initialize(m)'

    m.fs.softening_unit.report()
    m.fs.acidification_unit.report()
    solver.solve(m, tee=True)
    m.fs.softening_unit.report()
    m.fs.acidification_unit.report()
    # assert False
    m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"].unfix()
    m.fs.softening_unit.precipitation_reactor.pH["outlet"].fix(7.5)
    m.fs.acidification_unit.chemical_reactor.reagent_dose["HCl"].unfix()
    m.fs.acidification_unit.chemical_reactor.pH["outlet"].fix(6.5)
    m.fs.softening_unit.precipitation_block.update_block_scaling()
    m.fs.acidification_unit.chemistry_block.update_block_scaling()
    m.fs.softening_unit.precipitation_block.update_jacobian_scaling()
    m.fs.acidification_unit.chemistry_block.update_jacobian_scaling()
    # m.fs.ro_unit.deactivate_scaling_constraints()
    # m.fs.water_recovery.fix(0.9)
    # solve_model(m)
    solver.solve(m, tee=True)
    m.fs.softening_unit.report()
    for p in np.linspace(7.6, 11.5, 15):
        m.fs.softening_unit.precipitation_reactor.pH["outlet"].fix(p)
        print(f"\n\n------------Solving for pH: {p}------------")
        result = solver.solve(m, tee=True)
        assert_optimal_termination(result)
        m.fs.softening_unit.report()
        m.fs.acidification_unit.report()
        m.fs.softening_unit.precipitation_block.update_block_scaling()
        m.fs.acidification_unit.chemistry_block.update_block_scaling()
        print("------vars_close_to_bound-tests---------")
        print_variables_close_to_bounds(m)
        print("------constraints_close_to_bound-tests---------")

        print_constraints_close_to_bounds(m)
        # m.fs.softening_unit.precipitation_block.update_jacobian_scaling()
        # m.fs.acidification_unit.chemistry_block.update_jacobian_scaling()
    # for c in range(180, 220, 1):
    #     m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"].fix(c / 1000)
    #     print(f"\n\n------------Solving for CaO dose: {c} ppm------------")
    #     solver.solve(m, tee=True)
    #     m.fs.softening_unit.report()
    #     m.fs.acidification_unit.report()
    #     m.fs.softening_unit.precipitation_block.update_block_scaling()
    #     m.fs.acidification_unit.chemistry_block.update_block_scaling()
    #     m.fs.softening_unit.precipitation_block.update_jacobian_scaling()
    #     m.fs.acidification_unit.chemistry_block.update_jacobian_scaling()


def main():
    # test_softeining()
    main_ro_model()


def main_ro_model():
    m = build_model("../water_sources/USDA_brackish.yaml", hpro=False)
    initialize(m)

    # m.fs.water_recovery.fix(73 / 100)

    # solve_model(m)
    # for c in range(190, 220, 1):
    #     m.fs.softening_unit.precipitation_reactor.reagent_dose["CaO"].fix(c / 1000)
    #     print(f"\n\n------------Solving for CaO dose: {c} ppm------------")
    #     solve_model(m)
    #     process_reaktoro_blocks(m)
    #     print("------vars_close_to_bound-tests---------")
    #     print_variables_close_to_bounds(m)
    #     print("------constraints_close_to_bound-tests---------")

    #     print_constraints_close_to_bounds(m)

    for r in np.linspace(73, 90, 18):
        process_reaktoro_blocks(m)
        m.fs.water_recovery.fix(r / 100)
        print(f"\n\n------------Solving for water recovery: {r}%------------")
        solve_model(m)

        print("------vars_close_to_bound-tests---------")
        print_variables_close_to_bounds(m)
        print("------constraints_close_to_bound-tests---------")

        print_constraints_close_to_bounds(m)


def build_lime_only(water_case):
    mcas_props, feed_specs = get_source_water_data(water_case)
    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m = ConcreteModel()
    rkt_options = enable_multi_process_reaktoro(m)
    m.fs = FlowsheetBlock()
    # m.fs.costing = WaterTAPCosting()
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.ro_properties = SeawaterParameterBlock()
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        charge_balance_with_reaktoro=True,
        **feed_specs,
    )
    m.fs.softening_unit = PrecipitationUnit(
        default_property_package=m.fs.properties,
        # default_costing_package=m.fs.costing,
        selected_precipitants=["Calcite"],
        selected_reagents=[
            "CaO",
        ],
        reaktoro_options=rkt_options,
    )
    m.fs.acidification_unit = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        # default_costing_package=m.fs.costing,
        selected_reagents=["HCl"],
        reaktoro_options=rkt_options,
    )
    m.fs.feed.outlet.connect_to(m.fs.softening_unit.inlet)
    m.fs.softening_unit.outlet.connect_to(m.fs.acidification_unit.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    m.reaktoro_manager.build_reaktoro_blocks()
    m.fs.feed.fix_and_scale()
    m.fs.softening_unit.fix_and_scale()
    m.fs.acidification_unit.fix_and_scale()
    iscale.calculate_scaling_factors(m)

    m.fs.feed.initialize()
    m.fs.softening_unit.initialize()
    m.fs.acidification_unit.initialize()
    return m


def enable_multi_process_reaktoro(m):
    """Enables use of parallel solves for reaktoro blocks,
    in RO mode there will be 6 reaktoro blocks
        2 for Softening,
        2 for Acidification,
        2 for RO
        requiring 6 logical cores
    in HPRO mode there will be 8 reaktoro blocks
        2 for Softening,
        2 for Acidification,
        2 for RO,
        2 for HPRO
        requiring 8 logical cores
    """
    from reaktoro_pse.parallel_tools.reaktoro_block_manager import (
        ReaktoroBlockManager,
    )

    rkt_options = {}
    m.reaktoro_manager = ReaktoroBlockManager()
    rkt_options["reaktoro_block_manager"] = m.reaktoro_manager

    return rkt_options


def build_model(water_case, multi_process_reaktoro=True, hpro=False):
    rkt_options = None

    mcas_props, feed_specs = get_source_water_data(water_case)
    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m = ConcreteModel()
    if multi_process_reaktoro:
        rkt_options = enable_multi_process_reaktoro(m)
    m.fs = FlowsheetBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.ro_properties = SeawaterParameterBlock()
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        charge_balance_with_reaktoro=True,
        **feed_specs,
    )
    m.fs.softening_unit = PrecipitationUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        selected_precipitants=["Calcite"],
        selected_reagents=["CaO", "Na2CO3"],
        reaktoro_options=rkt_options,
    )

    m.fs.acidification_unit = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        selected_reagents=["HCl"],
        reaktoro_options=rkt_options,
    )

    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        initialization_pressure="osmotic_pressure",
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
            maximum_pressure=300 * pyunits.bar,
        )
        m.fs.hpro_unit = MultiCompROUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            ro_property_package=m.fs.ro_properties,
            selected_scalants={"Calcite": 1, "Gypsum": 1},
            reaktoro_options=rkt_options,
            default_costing_package_kwargs={
                "costing_method_arguments": {"RO_type": "high_pressure"}
            },
        )
        m.fs.product_mixer = MixerPhUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            inlet_ports=["ro", "hpro"],
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

        m.fs.ro_unit.product.connect_to(m.fs.product_mixer.inlet["ro"])
        m.fs.hpro_unit.product.connect_to(m.fs.product_mixer.inlet["hpro"])
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
    m.fs.lcow_objective = Objective(expr=(10 * m.fs.costing.LCOW) ** 2)
    TransformationFactory("network.expand_arcs").apply_to(m)
    add_global_constraints(m)
    fix_and_scale(m)
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


def report_global_state(m):
    data_dict = {"Global results": {}}
    data_dict["Global results"]["DOfs"] = int(degrees_of_freedom(m))
    data_dict["Global results"]["Water recovery"] = m.fs.water_recovery
    data_dict["Global results"]["LCOW"] = m.fs.costing.LCOW

    build_report_table("Global results", data_dict)


def fix_and_scale(m):
    for unit in m.flowsheet_unit_order:
        unit.fix_and_scale()

    iscale.calculate_scaling_factors(m)
    assert degrees_of_freedom(m) == 0


def initialize(m):
    for unit in m.flowsheet_unit_order:
        unit.initialize()
    m.fs.costing.initialize()
    report_all_units(m)

    solve_model(m)
    set_optimization(m)
    m.fs.water_recovery.fix(0.5)
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


def process_reaktoro_blocks(m):
    # pass
    for block in m.fs.reaktoro_blocks:
        block.update_block_scaling()
        block.update_jacobian_scaling()


def solve_model(m, solver=None):
    if solver is None:
        solver = get_cyipopt_solver(max_iters=100)
    result = solver.solve(m, tee=True)
    report_all_units(m)
    assert_optimal_termination(result)
    return result


if __name__ == "__main__":
    main()

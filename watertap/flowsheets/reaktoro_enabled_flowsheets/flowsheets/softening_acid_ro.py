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


def main():
    pass


def build_model(water_case, hpro=False):
    mcas_props, feed_specs = get_source_water_data(water_case)
    mcas_props["activity_coefficient_model"] = ActivityCoefficientModel.ideal
    mcas_props["density_calculation"] = DensityCalculation.constant

    m = ConcreteModel()

    m.fs = FlowsheetBlock()
    m.fs.costing = WaterTAPCosting()
    m.fs.properties = MCASParameterBlock(**mcas_props)
    m.fs.feed = MultiCompFeed(
        default_property_package=m.fs.properties,
        charge_balance_with_reaktoro=True,
        **feed_specs,
    )
    m.fs.softening_unit = PrecipitationUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        viable_precipitants=["Calcite"],
        viable_reagents=["CaO", "Na2CO3"],
    )

    m.fs.acidification_unit = ChemicalAdditionUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        viable_precipitants=["Calcite"],
        viable_reagents=["HCl", "H2SO4"],
    )

    m.fs.pump_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        initialization_pressure="osmotic_pressure",
        maximum_pressure=85 * pyunits.bar,
    )

    m.fs.ro_unit = MultiCompPumpUnit(
        default_property_package=m.fs.properties,
        default_costing_package=m.fs.costing,
        initialization_pressure="osmotic_pressure",
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
        m.fs.hpro_unit = MultiCompPumpUnit(
            default_property_package=m.fs.properties,
            default_costing_package=m.fs.costing,
            initialization_pressure="osmotic_pressure",
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

    # build all connections
    m.fs.feed.outlet.connect_to(m.fs.softening_unit.inlet)
    m.fs.softening_unit.outlet.connect_to(m.fs.acidification_unit.inlet)
    m.fs.acidification_unit.outlet.connect_to(m.fs.pump_unit.inlet)
    m.fs.pump_unit.outlet.connect_to(m.fs.ro_unit.inlet)

    if hpro:
        m.fs.ro_unit.retentate.connect_to(m.fs.hp_pump_unit.inlet)
        m.fs.hp_pump_unit.outlet.connect_to(m.fs.hpro_unit.inlet)
        m.fs.hpro_unit.outlet.connect_to(m.fs.erd_unit.inlet)

        m.fs.ro_unit.product.connect_to(m.fs.product_mixer.inlet["ro"])
        m.fs.hpro_unit.product.connect_to(m.fs.product_mixer.inlet["hpro"])
        m.fs.product_mixer.outlet.connect_to(m.fs.product.inlet)
    else:
        m.fs.ro_unit.outlet.connect_to(m.fs.erd_unit.inlet)
        m.fs.ro_unit.product.connect_to(m.fs.product.inlet)
    m.fs.erd_unit.outlet.connect_to(m.fs.brine.inlet)

    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.product.product.properties[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.product.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)
    m.fs.lcow_objective = Objective(expr=m.fs.costing.LCOW)
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def add_global_constraints(m):
    m.fs.water_recovery = Var(
        initialize=0.5,
        bounds=(0.2, 0.991),
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


def fix_and_scale(m):
    for unit in m.flowsheet_unit_order:
        unit.fix_and_scale()
    iscale.calculate_scaling_factors(m)


def initialize(m):
    for unit in m.flowsheet_unit_order:
        unit.initialize()
    report_all_units(m)
    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    report_all_units(m)


def report_all_units(m):
    for unit in m.flowsheet_unit_order:
        unit.report()


def set_optimization(m):
    for unit in m.flowsheet_unit_order:
        unit.set_optimization()
    report_all_units(m)
    m.fs.water_recovery.fix(0.5)
    solver = get_cyipopt_solver()
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)
    report_all_units(m)


if __name__ == "__main__":
    main()

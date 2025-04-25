from idaes.core.util.initialization import propagate_state
from pyomo.network import Arc
from pyomo.environ import (
    Constraint,
)
import idaes.core.util.scaling as iscale


class PortContainer:
    """Port container will aggregate ports and variables
    user wants to connect between unit models
    Although, the variable can be registered to the port directly, it can lead
    to issues when using UI visualization tools"""

    def __init__(self, name, port, var_dict, unit_block_reference):
        self.name = name
        self.port = port
        self.unit_block_reference = unit_block_reference
        if isinstance(var_dict, dict) == False:
            raise TypeError(
                "Var dict must be a dictionary with structure {'variable name': pyomo.var}"
            )
        self.var_dict = var_dict

    def connect_to(self, inlet):
        """connect current port to provided port, registering
        the generated connection container using provided function"""
        connection = ConnectionContainer(self, inlet)
        self.unit_block_reference.register_outlet_connection(connection)


class ConnectionContainer:
    """This serves as container for any ports and equality constraints for use between unit
    models"""

    def __init__(self, outlet, inlet):
        self.registered_equality_constraints = []
        if not isinstance(outlet, PortContainer) or not isinstance(
            inlet, PortContainer
        ):
            raise TypeError("Provided outlet and inlet must be PortContainer objects")
        self.build_arc(outlet, inlet)
        self.build_constraints(outlet, inlet)

    def build_arc(self, outlet, inlet):
        """
        builds a standard arc while naming it with outlet and inlet name, this should
        be a unique pair always (e.g. for a unit model you should only havbe single outlet-> inlet conenction)
        """
        self.registered_arc = outlet.unit_block_reference.add_component(
            f"{outlet.name}_to_{inlet.name}",
            Arc(source=outlet.port, destination=inlet.port),
        )
        self.registered_arc = outlet.unit_block_reference.find_component(
            f"{outlet.name}_to_{inlet.name}"
        )

    def build_constraints(self, outlet, inlet):
        """
        builds equality constraints for provided variables and scales them.
        We ensure that oulet and inlet vars are the same.
        """
        if set(outlet.var_dict) != set(inlet.var_dict):
            raise KeyError(
                f"Provided inlet keys: {outlet.var_dict} do not match outlet keys: {inlet.var_dict}"
            )

        for outlet_key in outlet.var_dict:
            outlet.unit_block_reference.add_component(
                f"eq_{outlet_key}_{outlet.name}_to_{inlet.name}",
                Constraint(
                    expr=outlet.var_dict[outlet_key] == inlet.var_dict[outlet_key]
                ),
            )
            constraint = outlet.unit_block_reference.find_component(
                f"eq_{outlet_key}_{outlet.name}_to_{inlet.name}"
            )
            # ensure we register it for propagation later on
            self.registered_equality_constraints.append(
                (constraint, outlet.var_dict[outlet_key], inlet.var_dict[outlet_key])
            )

            sf = iscale.get_scaling_factor(outlet.var_dict[outlet_key])
            if sf != None:
                iscale.constraint_scaling_transform(
                    constraint,
                    sf,
                )

    def propagate(self):
        """this should prop any arcs and also ensure all equality constraints are satisfied"""
        propagate_state(self.registered_arc)
        self.propagate_equality_constraints()

    def propagate_equality_constraints(self):
        for (
            constraint,
            outlet_var,
            inlet_var,
        ) in self.registered_equality_constraints:
            inlet_var.value = outlet_var.value

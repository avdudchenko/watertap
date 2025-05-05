from pyomo.environ import (
    value,
    Var,
    Constraint,
    units as pyunits,
)


import idaes.core.util.scaling as iscale


class ViableReagentsBase(dict):
    """class for tracking reagents and creating approriate constraints to handle non-pure reagents"""

    def register_reagent(
        self,
        reagent,
        mw,
        dissolution_stoichiometric,
        density_reagent=1 * pyunits.kg / pyunits.L,
        min_dose=0.1,
        max_dose=3000,
        cost=None,
        purity=1,
        solvent=None,
        reagent_modifier_dict=None,
        solvent_modifier_dict={"H2O": {"H": 2, "O": 1}},
    ):
        """
        Add a new reagent to default list:

        Args:
            reagent - name of reagent
            mw - molecular weight of reagent (include pyomo units)
            dissolution_stoichiometric - dictionary that defines which ions the reagent dissociates into {'ion':moles}
            density_reagent - density of reagent
            cost - (optional) if costing is used, provide cost for reagent in costing package base units per kg.
            min_dose (default: 0.1 PPM) - minimum reagent dose
            max_dose (default: 3000 PPM) - maximum reagent dose
            purity: Purity of reagent based on w/w basis
                if less then 1, it will add solvent as one of dissolution species, and estimate how many moles
                of solvent is added per mole of reagent based on purity.
            solvent: Solvent information as a tuple ('solvent',mw with pyomo units)
                example ('H2O':18.01*pyunts.g/pyunits.mol)
            reagent_modifier_dict: dictionary for defining a new reaktoro chemistry modifier should be in form of
                {"reagent":{'element':mols, element_2:mols}, e.g. "Ca(OH)2": {"Ca": 1, "O": 2, "H":2}
            solvent_modifier_dict: dictionary for defining a new reaktoro chemistry modifier should be in form of
                {"reagent":{'element':mols, element_2:mols}, e.g. "solvent_H2O": {"H": 2, "O": 1}
        """

        if purity < 1:
            # if reagent already dissocaites into solvent track amount
            # (eg NaOH dissociates into Na, and H2O, so if purity is 50% by weight
            # we need to track amount of solvent + original H2O in pure chemical)
            _solvent_adjust = 0
            if solvent[0] in dissolution_stoichiometric:
                _solvent_adjust = dissolution_stoichiometric[solvent[0]]
            mols_reagent = purity / value(pyunits.convert(mw, pyunits.g / pyunits.mol))
            mols_solvent = (1 - purity) / value(
                value(pyunits.convert(solvent[1], pyunits.g / pyunits.mol))
            )
            # calculate solution mw
            mw = mw + solvent[1] * mols_solvent / mols_reagent
            # ratio on mol to mol basis
            solvent_ratio = mols_solvent / mols_reagent
            dissolution_stoichiometric[solvent[0]] = solvent_ratio + _solvent_adjust
        self[reagent] = {
            "mw": mw,
            "dissolution_stoichiometric": dissolution_stoichiometric,
            "cost": cost,
            "min_dose": min_dose,  # ppm
            "max_dose": max_dose,  # ppm
            "purity": purity,
            "solvent": solvent,
            "density_reagent": density_reagent,
        }
        if reagent_modifier_dict is not None or solvent_modifier_dict is not None:
            if hasattr(self, "modifier") is False:
                self.modifier = {}
            if reagent_modifier_dict is not None:
                self.modifier.update(reagent_modifier_dict)
            if solvent_modifier_dict is not None:
                self.modifier.update(solvent_modifier_dict)
        if solvent is not None:
            if hasattr(self, "solvents") is False:
                self.solvents = {}

            self.solvents[reagent] = {
                "solvent": solvent[0],
                "reagent": reagent,
                "solvent_ratio": solvent_ratio,  # for reaktoro only - we need to know extra water, reaktoro handles dissociation for us
            }

    def get_reaktoro_chemistry_modifiers(self):
        """Return a dictionary of reaktoro chemistry modifiers"""
        if hasattr(self, "modifier"):
            return self.modifier
        else:
            return None

    def check_if_reagents_are_pure(self, reagent_var):
        """Check if any of the reagents are not pure"""
        for reagent in reagent_var:
            if reagent in self.solvents:
                return False
        return True

    def get_unqiue_solvents(self, reagent_var):
        """Get a list of unique solvents from the reagents"""
        solvents = []
        for reagent in reagent_var:
            if reagent in self.solvents:
                if self.solvents[reagent]["solvent"] not in solvents:
                    solvents.append(self.solvents[reagent]["solvent"])
        return solvents

    def create_solvent_constraint(self, block, reagent_var, time_unit=pyunits.s):
        """Create a constraint for the solvent amount based on the reagent amount

        This will create a new variable for the solvent amount (flow_mol_solvent) and a constraint (eq_flow_mol_solvent)
        that relates the solvent amount to the reagent amount.
        The solvent amount is calculated based on the solvent ratio and the reagent amount.
        """
        self.solvent_info = {}

        if (
            hasattr(self, "solvents")
            and self.check_if_reagents_are_pure(reagent_var) == False
        ):
            active_solvents = self.get_unqiue_solvents(reagent_var)
            block.add_component(
                "flow_mol_solvent",
                Var(
                    active_solvents,
                    units=pyunits.mol / time_unit,
                ),
            )

            @block.Constraint(active_solvents)
            def eq_flow_mol_solvent(blk, solvent):
                total_solvent = []
                self.solvent_info[solvent] = []
                for reagent in reagent_var:
                    if solvent == self.solvents[reagent]["solvent"]:
                        total_solvent.append(
                            reagent_var[reagent]
                            * self.solvents[reagent]["solvent_ratio"]
                        )
                        # grab refernce for each var so we can scale things later
                        self.solvent_info[solvent].append(reagent_var[reagent])
                if len(total_solvent) >= 1:
                    return block.find_component("flow_mol_solvent")[solvent] == sum(
                        total_solvent
                    )
                else:
                    return Constraint.skip()

            return block.find_component("flow_mol_solvent")
        else:
            return None

    def scale_solvent_vars_and_constraints(self, block, reagent_var):
        for solvent, info in self.solvent_info.items():
            sfs = []
            for var in info:
                sf = iscale.get_scaling_factor(var)
                sfs.append(sf)
            sff = max(sfs)  # lets grab maximum scaling factor
            print(solvent, sff, sfs)
            iscale.set_scaling_factor(
                block.find_component("flow_mol_solvent")[solvent], sff
            )
            iscale.constraint_scaling_transform(
                block.find_component("eq_flow_mol_solvent")[solvent], sff
            )


class ViablePrecipitantsBase(dict):
    def register_solid(
        self, precipitant, mw, precipitation_stoichiometric, primary_ion
    ):
        """
        Add a new reagent to default list:

        Args:
            reagent - name of reagent
            mw - molecular weight of reagent (include pyomo units)
            precipitation_stoichiometric - dictionary that contains what species form the solid {'ion':moles}
            primary_ion - a primary ion that defines formation of this precipitants

        """

        self[precipitant] = {
            "mw": mw,
            "precipitation_stoichiometric": precipitation_stoichiometric,
            "primary_ion": primary_ion,
        }


class ReaktoroOptionsContainer(dict):
    """Container for storing reaktoro options and providing simple methods for updating
    options dictionary with user provided options"""

    def __init__(self):
        self["system_state"] = {}
        self["aqueous_phase"] = {
            "activity_model": "ActivityModelPitzer",
            "fixed_solvent_specie": "H2O",
            "convert_to_rkt_species": True,
        }
        self["database"] = "PhreeqcDatabase"
        self["database_file"] = "pitzer.dat"
        self["reaktoro_block_manager"] = None
        self["build_speciation_block"] = True
        self["assert_charge_neutrality"] = True
        self["dissolve_species_in_reaktoro"] = True

    def system_state_option(self, option, value):
        self["system_state"][option] = value

    def system_state_modifier_option(self, option, value):
        self["system_state"][option] = value

    def aqueous_phase_option(self, option, value):
        self["aqueous_phase"][option] = value

    def update_with_user_options(self, options):
        if options is None:
            pass
        else:
            self.update(options)

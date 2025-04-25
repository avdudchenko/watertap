from pyomo.environ import (
    value,
    units as pyunits,
)


class ViableReagentsBase(dict):
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
            purity: Purity of reagent
            solvent: Solvent information as a tuple ('solvent',mw with pyomo units)
                example ('H2O':18.01*pyunts.g/pyunits.mol)

        """

        if purity < 1:
            _solvent_adjust = 0
            if solvent[0] in dissolution_stoichiometric:
                _solvent_adjust = dissolution_stoichiometric[solvent[0]]
            mols_reagent = purity / value(pyunits.convert(mw, pyunits.g / pyunits.mol))
            mols_solvent = (1 - purity) / value(
                value(pyunits.convert(solvent[1], pyunits.g / pyunits.mol))
            )
            # calculate solution mw
            mw = mw + solvent[1] * mols_solvent / mols_reagent
            dissolution_stoichiometric[solvent[0]] = mols_solvent / mols_reagent
        self[reagent] = {
            "mw": mw,
            "dissolution_stoichiometric": dissolution_stoichiometric,
            "cost": cost,
            "min_dose": min_dose,  # ppm
            "max_dose": max_dose,  # ppm
            "purity": purity,
            "solvent": solvent,
        }


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

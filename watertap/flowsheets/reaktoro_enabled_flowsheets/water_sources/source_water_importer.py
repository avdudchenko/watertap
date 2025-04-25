import yaml
from pyomo.environ import (
    units as pyunits,
)


def get_source_water_data(file_location):
    """simple function to load feed water compostion from yaml file"""
    with open(file_location, "r") as ymlfile:
        data_dict = yaml.safe_load(ymlfile)
    # Converts yaml structure to dict structure for use with MCAS
    mcas_param_dict = {}
    mcas_param_dict["solute_list"] = get_solute_dict(data_dict)
    mcas_param_dict["diffusivity_data"] = gen_diffusivity_dict(data_dict)
    mcas_param_dict["mw_data"] = gen_mw_dict(data_dict)
    mcas_param_dict["stokes_radius_data"] = gen_stoke_dict(data_dict)
    mcas_param_dict["charge"] = gen_charge_dict(data_dict)

    # Creats dict with feed properties to pass into multi_comp_feed
    mass_comp_dict = get_feed_comp(data_dict)
    pH = float(data_dict["pH"])
    mass_flowrate = data_dict.get("flow_mass") * pyunits.kg / pyunits.s
    feed_temperature = data_dict.get("temperature", 293.15)
    feed_spec_dict = {
        "ion_concentrations": mass_comp_dict,
        "pH": pH,
        "temperature": feed_temperature,
        "mass_flowrate": mass_flowrate,
    }
    return mcas_param_dict, feed_spec_dict


def get_solute_dict(data_dict):
    solute_list = list(data_dict["solute_list"].keys())
    return solute_list


def gen_diffusivity_dict(data_dict):
    diff_dict = {}
    for solute in data_dict["solute_list"].keys():
        diff_dict[("Liq", solute)] = float(
            data_dict["solute_list"][solute]["diffusivity"]
        )
    return diff_dict


def gen_mw_dict(data_dict):
    mw_dict = {}
    for solute in data_dict["solvent_list"].keys():
        mw_dict[solute] = float(
            data_dict["solvent_list"][solute]["molecular_weight (kg/mol)"]
        )
    for solute in data_dict["solute_list"].keys():
        mw_dict[solute] = float(
            data_dict["solute_list"][solute]["molecular_weight (kg/mol)"]
        )
    return mw_dict


def gen_stoke_dict(data_dict):
    stokes_dict = {}
    for solute in data_dict["solute_list"].keys():
        stokes_dict[solute] = float(
            data_dict["solute_list"][solute]["stokes_radius (m)"]
        )
    return stokes_dict


def gen_charge_dict(data_dict):
    charge_dict = {}
    for solute in data_dict["solute_list"].keys():
        charge_dict[solute] = float(
            data_dict["solute_list"][solute]["elemental charge"]
        )
    return charge_dict


def get_feed_comp(data_dict):
    mass_loading_dict = {}
    for solute in data_dict["solute_list"].keys():
        mass_loading_dict[solute] = (
            float(data_dict["solute_list"][solute]["concentration (mg/L)"])
            * pyunits.mg
            / pyunits.L
        )

    return mass_loading_dict

import datetime
import os

import yaml

import lephare as lp

__all__ = [
    "prepare",
    "overwrite_config",
    "read_yaml_config",
    "write_yaml_config",
    "write_para_config",
    "keymap_to_string_dict",
    "string_dict_to_keymap",
    "all_types_to_keymap",
]

_DEFAULT = object()


def prepare(config, star_config=_DEFAULT, gal_config=_DEFAULT, qso_config=_DEFAULT):
    """Run the prepare stages of LePHARE

    In order to run "zphota" we must create the filter files, run sedtolib to
    create the SED libraries, and finally run mag_gal, to create the
    magnitude libararies. We abstract these tasks into a single prepare stage.

    We overide the config for each type if distinct config set. If no overide
    configs are set we use the same for each type.

    Parameters
    ==========
    config : dict of lephare.keyword
        The config base to run all tasks
    star_config : dict of lephare.keyword or None
        Config values to override for stars. If None do not run.
    gal_config : dict of lephare.keyword or None
        Config values to override for galaxies. If None do not run.
    qso_config : dict of lephare.keyword or None
        Config values to override for QSO. If None do not run.
    """
    # ensure that all values in the keymap are keyword objects
    config = all_types_to_keymap(config)
    object_types = dict(STAR=star_config, GAL=gal_config, QSO=qso_config)

    # Run the filter command
    # load filters from config
    filter_lib = lp.FilterSvc.from_keymap(config)
    # Get location to store filter files
    filter_output = os.path.join(os.environ["LEPHAREWORK"], "filt", config["FILTER_FILE"].value)
    # Write filter files
    lp.write_output_filter(filter_output + ".dat", filter_output + ".doc", filter_lib)
    # # Write config to work directory
    write_yaml_config(config, f"{filter_output}_config.yaml")

    for object_type in object_types:
        # Skip object types that are not to be processed
        if object_types[object_type] is None:
            continue
        elif object_types[object_type] is _DEFAULT:
            updated_config = config
        else:
            # Overwrite config for this object type
            object_types[object_type] = all_types_to_keymap(object_types[object_type])
            updated_config = overwrite_config(config, object_types[object_type])

        # Write the updated config
        sed_out_name = f"{updated_config[f'{object_type}_LIB'].value}_{object_type.lower()}_config.yaml"
        sed_output = os.path.join(os.environ["LEPHAREWORK"], "lib_bin", sed_out_name)
        mag_out_name = f"{updated_config[f'{object_type}_LIB_OUT'].value}_{object_type.lower()}_config.yaml"
        mag_output = os.path.join(os.environ["LEPHAREWORK"], "lib_mag", mag_out_name)
        # Run sedtolib
        sedlib = lp.Sedtolib(config_keymap=updated_config)
        list_loc = updated_config[f"{object_type}_SED"].value
        # if find sed/ in the path, assume it is in lephare-data and not an absolute path
        if list_loc.find("sed/") >= 0:
            list_loc = os.path.join(lp.LEPHAREDIR, list_loc)
        sedtolib_kwargs = {f"{object_type.lower()}_sed": list_loc}
        sedlib.run(typ=object_type, **sedtolib_kwargs)
        write_yaml_config(updated_config, sed_output)
        # Run mag_gal
        maglib = lp.MagGal(config_keymap=updated_config)
        maglib.run(typ=object_type, verbose=False)
        write_yaml_config(updated_config, mag_output)


def overwrite_config(config1, config2):
    """Check that two config can be safely broadcast for a joint run"""
    if config2 is None:
        return config1
    config = config1
    # Redshift grid must be idenitical
    # assert True
    # filters must be indentical
    # assert True
    for k in config2:
        config[k] = config2[k]
    return config


def read_yaml_config(yaml_file_path):
    """Open a standard yaml file and render it as a dictionary of keywords

    Parameters
    ==========
    yaml_file_path : str
        Path to input yaml file.
    """
    with open(yaml_file_path, "r") as file:
        config = yaml.safe_load(file)
    for key in config:
        config[key] = lp.keyword(key, config[key])
    return config


def write_yaml_config(keymap, yaml_file_path):
    """Write a dictionary of keywords to a yaml file

    Parameters
    ==========
    keymap : dict of lephare.keyword
        The dictionary of keywords to be written to yaml.
    yaml_file_path : str
        Path to output yaml file.
    """
    keymap = lp.all_types_to_keymap(keymap)
    config_dict = {}
    for k in keymap:
        config_dict[keymap[k].name] = keymap[k].value
    with open(yaml_file_path, "w") as yaml_file:
        now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        yaml_file.write(f"# File written automatically at {now} by lephare version {lp.__version__}\n")
        yaml.dump(config_dict, yaml_file)


def write_para_config(keymap, para_file_path):
    """Write a dictionary of keywords to a para file

    Parameters
    ==========
    keymap : dict of lephare.keyword
        The dictionary of keywords to be written to yaml.
    para_file_path : str
        Path to output para file.
    """
    keymap = lp.all_types_to_keymap(keymap)
    now = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    para_contents = f"# File written automatically at {now} by lephare version {lp.__version__}\n"
    for k in keymap:
        para_contents += f"{k} {keymap[k].value}\n"
    with open(para_file_path, "w") as file_handle:
        file_handle.write(para_contents)
        file_handle.close()


def keymap_to_string_dict(keymap):
    """Convert a dictionary of keywords to a dictionary of strings"""
    config_dict = {}
    for k in keymap:
        config_dict[keymap[k].name] = keymap[k].value
    return config_dict


def string_dict_to_keymap(string_dict):
    """Convert a dictionary of strings to a dictionary of keywords"""
    keymap = {}
    for k in string_dict:
        keymap[k] = lp.keyword(k, str(string_dict[k]))
    return keymap


def all_types_to_keymap(input_config):
    """Convert all types of configs to keymap"""
    keymap = {}
    for k in input_config:
        if type(input_config[k]) == lp.keyword:
            keymap[k] = input_config[k]
        else:
            keymap[k] = lp.keyword(k, str(input_config[k]))
    return keymap

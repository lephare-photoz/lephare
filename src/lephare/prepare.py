import os

import yaml

import lephare as lp

__all__ = ["prepare", "overwrite_config", "read_yaml_config", "write_yaml_config"]


def prepare(config, star_config=None, gal_config=None, qso_config=None):
    """Run the prepare stages of LePHARE

    In order to run "zphota" we must create the filter files, run sedtolib to
    create the SED libraries, and finally run mag_gal, to create the
    magnitude libararies. We abract these tasks into a single prepare stage.

    We overide the config for each type if distinct config set. If no overide
    configs are set we use the same for each type.

    Parameters
    ==========
    config : dict of lephare.keyword
        The config base to run all tasks
    star_config : dict of lephare.keyword or None
        Config values to override for stars
    gal_config : dict of lephare.keyword or None
        Config values to override for galaxies
    qso_config : dict of lephare.keyword or None
        Config values to override for QSO.
    """
    object_types = {"STAR": star_config, "GAL": gal_config, "QSO": qso_config}
    # Run the filter command
    # load filters from config file
    filter_lib = lp.FilterSvc.from_config(config)
    # Get location to store filter files
    filter_output = os.path.join(os.environ["LEPHAREWORK"], "filt", config["FILTER_FILE"].value)
    # Write filter files
    lp.write_output_filter(filter_output + ".dat", filter_output + ".doc", filter_lib)
    # Write config to work directory
    write_yaml_config(config, f"{filter_output}_config.yaml")

    for object_type in object_types:
        updated_config = overwrite_config(config, object_types[object_type])
        # Run sedtolib
        sedlib = lp.Sedtolib(config_keymap=updated_config)
        sedlib.run(typ=object_type)
        # Write config
        sed_output = os.environ["LEPHAREWORK"] + "sedtolib"
        write_yaml_config(updated_config, f"{sed_output}_config.yaml")
        # Run mag_gal
        maglib = lp.MagGal(config_keymap=updated_config)
        maglib.run(typ=object_type)
        # Write config
        mag_output = os.environ["LEPHAREWORK"] + "mag_gal"
        write_yaml_config(updated_config, f"{mag_output}_config.yaml")


def overwrite_config(c1, c2):
    """Check that two config can be safely broadcast for a joint run"""
    if c2 is None:
        return c1
    config = c1
    # Redshift grid must be idenitical
    # assert True
    # filters must be indentical
    # assert True
    for k in c2:
        config[k] = config[k]
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
    config_dict = {}
    for k in keymap:
        config_dict[keymap[k].name] = keymap[k].value
    with open(yaml_file_path, "w") as yaml_file:
        yaml.dump(config_dict, yaml_file)

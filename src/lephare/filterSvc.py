# ruff: noqa: N999

import os
import xml.dom.minidom

import requests
import yaml

from lephare import LEPHAREDIR
from lephare._lephare import check_first_char, flt

__all__ = [
    "FilterSvc",
]

SVO_URL = "http://svo2.cab.inta-csic.es/theory/fps"


class FilterSvc:
    """Filter service for retrieving filters from various sources.

    This class defines a number of methods for loading and manipulating filters.
    """

    @classmethod
    def from_yaml(cls, yaml_file):
        """Load filter from a yaml file

        Parameters
        ----------
        yaml_file : `str`
            Path to yaml file

        Returns
        -------
        flt_array :  list of `lephare.flt`
            Array of wavelengths in Angstrom and filter transmissivity.
        """
        config = yaml.load(open(yaml_file), Loader=yaml.BaseLoader)["filters"]  # noqa: SIM115
        if "calib" in config:
            default_calib = config["calib"]
        if "trans" in config:
            default_trans = config["trans"]
        flt_array = []
        counter = 0
        for entry in config["list"]:
            name = entry["name"]
            # filters are ordered starting from one (FORTRAN legacy)
            counter += 1  # noqa: SIM113
            calib = entry.get("calib", default_calib)
            trans = entry.get("trans", default_trans)
            if name[:4] == "svo:":
                flt_obj = cls.from_svo(counter, name[4:], "AB", calib)
            else:
                flt_obj = cls.from_file(name, counter, trans, calib)
            flt_array.append(flt_obj)
        return flt_array

    @classmethod
    def from_config(cls, config_file):
        """Load filter from a .para config file

        Parameters
        ----------
        config_file : `str`
            Path to config file

        Returns
        -------
        flt_array : list of `lephare.flt`
            List of filters.
        """
        keywords = ["FILTER_REP", "FILTER_LIST", "TRANS_TYPE", "FILTER_CALIB", "FILTER_FILE"]
        keymap = {}
        with open(config_file) as fstream:
            count = 0
            while count < len(keywords):
                line = fstream.readline()
                for key in keywords:
                    if key in line and check_first_char(line):
                        keymap[key] = line.split()[1]
                        count += 1
        keymap["FILTER_REP"] = keymap["FILTER_REP"].replace("$LEPHAREDIR", LEPHAREDIR)

        filter_list = keymap["FILTER_LIST"].split(",")
        filter_calib = keymap["FILTER_CALIB"].split(",")
        filter_trans = keymap["TRANS_TYPE"].split(",")
        if len(filter_trans) == 1:
            filter_trans = len(filter_list) * filter_trans
        elif len(filter_trans) != len(filter_list):
            raise RuntimeError("FILTER_LIST and FILTER_TRANS do not have the same size")
        if len(filter_calib) == 1:
            filter_calib = len(filter_list) * filter_calib
        elif len(filter_calib) != len(filter_list):
            raise RuntimeError("FILTER_LIST and FILTER_CALIB do not have the same size")
        flt_array = []
        for i in range(len(filter_list)):
            name = os.path.join(keymap["FILTER_REP"], filter_list[i])
            oneflt = flt(i, name, int(filter_trans[i]), int(filter_calib[i]))
            flt_array.append(oneflt)
        return flt_array

    @classmethod
    def from_svo(cls, counter, filter_id, system="AB", calib=0):
        """Return filter from SVO

        Parameters
        ----------
        counter : `int`
            Filter number
        filter_id : `str`
            Id of filter in SVO format.
        system : `str`, optional
            Photometric system
        calib : `float`, optional
            Calibration value. Not currently used.

        Returns
        -------
        res : `lephare.flt`
            Filter in native lephare format.
        """
        res = FilterSvc.svo_request(counter, filter_id, system)
        return res

    @classmethod
    def from_file(cls, filename, counter=-1, trans=0, calib=0):
        """Return filter from SVO

        Parameters
        ----------
        filename : `str`
            Path to filter file.
        counter : `int`
            Filter number
        trans : `int`, optional
            Photometric system.
        calib : `int`, optional
            Calibration value.

        Returns
        -------
        f : `lephare.flt`
            Filter in native lephare format.
        """
        name = filename.replace("$LEPHAREDIR", LEPHAREDIR)
        f = flt(counter, name, trans, calib)
        return f

    @classmethod
    def svo_request(cls, counter, filter_id, system):
        """Retrieve a filter from the SVO

        Parameters
        ----------
        counter : `int`
            Filter number
        filter_id : `str`
            Id of filter in SVO format.
        system : `str`, optional
            Photometric system

        Returns
        -------
        res : `lephare.flt`
            Filter in native lephare format.
        """
        try:
            query = f"{SVO_URL}/fps.php?PhotCalID={filter_id}/{system}"
            r = requests.get(query, timeout=5)
        except ConnectionRefusedError:
            print(f"request {query} failed due to failure to connect to the server.")
            return None
        except requests.ConnectTimeout:
            print("Timeout on SVO server")
            return None
        try:
            dd = xml.dom.minidom.parseString(r.content)
        except xml.parsers.expat.ExpatError:
            print("SVO server down")
            return None
        # assert query ok
        for info in dd.getElementsByTagName("INFO"):
            if info.getAttribute("name") == "QUERY_STATUS" and info.getAttribute("value") != "OK":
                raise AssertionError(
                    f"QUERY_STATUS did not return OK; check the filter_id input: {SVO_URL}. "
                    f"For a list of valid filters, see {filter_id}"
                )
        # filter params: not used (yet?)
        params = dd.getElementsByTagName("PARAM")
        param_dict = {"Date": r.headers["Date"]}
        for param in params:
            param_dict[param.getAttribute("name")] = param.getAttribute("value")

        trans_type = int(param_dict["DetectorType"])
        name = filter_id.split("/")[1]

        # filter curve
        ##check units
        f1, f2 = dd.getElementsByTagName("FIELD")
        # these assert are meant to secure against non
        # standard entries in the SVO DB
        assert f1.getAttribute("unit") == "Angstrom"
        assert f1.getAttribute("datatype") == "double"
        assert f2.getAttribute("datatype") == "double"
        table = dd.getElementsByTagName("TABLEDATA")[0]
        data = table.getElementsByTagName("TR")
        with open("./" + name, "w") as stream:
            for _, d in enumerate(data):
                dd = d.getElementsByTagName("TD")
                l_val = float(dd[0].childNodes[0].data)
                t_val = float(dd[1].childNodes[0].data)
                stream.write(f"{l_val} {t_val}\n")

        # not super nice, but the only quick solution found to
        # avoid exposing trans and clean methods.
        flt_obj = flt(counter, "./" + name, trans_type, 0)
        # flt_obj.read("dummy")
        os.remove(name)
        # save the SVO additional info in an instance property
        flt_obj.svo_params = param_dict
        return flt_obj

    @classmethod
    def from_keymap(cls, keymap):
        """Load filter from a config keymap.

        Parameters
        ----------
        keymap : dict of lephare.keyword
            The config keymap

        Returns
        -------
        flt_array : list of `lephare.flt`
            List of filters.
        """
        keymap["FILTER_REP"].value = keymap["FILTER_REP"].value.replace("$LEPHAREDIR", LEPHAREDIR)

        filter_list = keymap["FILTER_LIST"].value.split(",")
        filter_calib = keymap["FILTER_CALIB"].value.split(",")
        filter_trans = keymap["TRANS_TYPE"].value.split(",")
        if len(filter_trans) == 1:
            filter_trans = len(filter_list) * filter_trans
        elif len(filter_trans) != len(filter_list):
            raise RuntimeError("FILTER_LIST and FILTER_TRANS do not have the same size")
        if len(filter_calib) == 1:
            filter_calib = len(filter_list) * filter_calib
        elif len(filter_calib) != len(filter_list):
            raise RuntimeError("FILTER_LIST and FILTER_CALIB do not have the same size")
        flt_array = []
        for i in range(len(filter_list)):
            name = os.path.join(keymap["FILTER_REP"].value, filter_list[i])
            oneflt = flt(i, name, int(filter_trans[i]), int(filter_calib[i]))
            flt_array.append(oneflt)
        return flt_array

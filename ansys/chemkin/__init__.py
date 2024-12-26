"""
Chemkin module for python
core
"""

# import numpy as np
from ctypes import c_int
import inspect
import os
import platform

import yaml

# import kinetics
from . import chemkin_wrapper as ck_wrapper
from .chemistry import *
from .color import Color
from .mixture import *
from .reactormodel import *

# try:
#    import importlib.metadata as importlib_metadata
# except ModuleNotFoundError:
#    import importlib_metadata

# __version__ = importlib_metadata.version(__name__.replace(".", "-"))

# show ansys (chemkin) version number
print(
    Color.YELLOW + f"chemkin version number = {ck_wrapper._ansys_ver:d}",
    end=Color.END,
)
# get ansys installation location
ansys_dir = str(ck_wrapper._ansys_dir)

_chemkin_platform = None
if platform.system() == "Windows":
    _chemkin_platform = "win64"
else:
    _chemkin_platform = "linuxx8664"

# get chemkin installation location
_chemkin_root = os.path.join(ansys_dir, "reaction", "chemkin" + _chemkin_platform)

# set default units to cgs
unit_code = c_int(1)
iError = ck_wrapper.chemkin.KINSetUnitSystem(unit_code)

# == Chemkin module global parameters -- DO NOT MODIFY without asking Chemkin development team members
frm = inspect.currentframe()
if frm is not None:
    _chemkinmodule_path = os.path.dirname(
        inspect.getfile(frm)
    )  # chemkin module home directory
_ChemkinHelpFile = (
    os.sep + "ChemkinKeywordTips.yaml"  # Chemkin keyword help data file in YAML format
)
_helploaded = False
# == end of global parameters
# print(_chemkinmodule_path + _ChemkinHelpFile)
# Chemkin keyword hints
if not _helploaded:
    # load Chemkin keyword dictionary from the YAML file
    with open(_chemkinmodule_path + _ChemkinHelpFile, "r") as hints:
        _CKdict = yaml.safe_load(hints)
        _helploaded = True


def keywordhints(mykey):
    """
    Get hints about the Chemkin keyword
    :param mykey: keyword (upper-case string)
    :return: None
    """
    # look up the keyword
    global _CKdict
    key = _CKdict.get(mykey.upper())
    if key is not None:
        # fetch the information about the keyword
        description, default, unit = key.values()
        # show the result
        print(Color.YELLOW + f"** tips about keyword '{mykey}'")
        print(f"     Description: {description}")
        print(f"     Default Value: {default}")
        print(f"     Units: {unit}", end=Color.END)
    else:
        print(Color.RED + f"** keyword '{mykey}' is not found", end=Color.END)


def phrasehints(phrase):
    """
    Get keyword hints by using key phrase in the description
    :param phrase: lower-case search string (string)
    :return: None
    """
    # initialization
    keys = []
    global _CKdict
    # search to find keyword descriptions that contain the phrase
    for s in _CKdict.values():
        if phrase.lower() in s.get("Description"):
            # get the dictionary index
            k = list(_CKdict.values()).index(s)
            # put the corresponding keywords into a candidate list
            keys.append(list(_CKdict.keys())[k])
    # show the hints for all candidate keywords
    if len(keys) > 0:
        for k in keys:
            keywordhints(k)
    else:
        print(
            Color.RED
            + f"** no keyword description containing the phrase '{phrase}' is found",
            end=Color.END,
        )

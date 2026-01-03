# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Chemkin help menu for keywords and key phrases."""

import importlib.resources
from pathlib import Path
from typing import Union
import webbrowser

import yaml

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.logger import logger

CKdict = {}  # chemkin hints
_help_loaded = False


def setup_hints():
    """Set up Chemkin keyword hints."""
    # Chemkin keyword help data file in YAML format
    _chemkin_resources_dir = importlib.resources.files("ansys.chemkin.core").joinpath(
        "data"
    )
    help_file = Path(_chemkin_resources_dir) / "ChemkinKeywordTips.yaml"
    global _help_loaded
    if not _help_loaded:
        global CKdict
        # load Chemkin keyword dictionary from the YAML file
        with Path.open(help_file, "r") as hints:
            CKdict = yaml.safe_load(hints)
            _help_loaded = True


def clear_hints():
    """Clear the Chemkin keyword data."""
    global _help_loaded
    global CKdict
    if _help_loaded:
        CKdict.clear()


def keyword_hints(mykey: str):
    """Get hints about the Chemkin keyword."""
    """Get hints about the Chemkin keyword.

    Parameters
    ----------
        mykey: string
            keyword phrase

    """
    # look up the keyword
    global CKdict
    key = CKdict.get(mykey.upper())
    if key is not None:
        # fetch the information about the keyword
        description, default, unit = key.values()
        # show the result
        print(f"** tips about keyword '{mykey}'")
        print(f"     Description: {description}")
        print(f"     Default Value: {default}")
        print(f"     Units: {unit}")
    else:
        msg = [Color.PURPLE, "keyword", mykey, "is not found.", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)


def phrase_hints(phrase: str):
    """Get keyword hints by using key phrase in the description."""
    """Get keyword hints by using key phrase in the description.

    Parameters
    ----------
    phrase: string
        search phrase

    """
    # initialization
    keys = []
    global CKdict
    # search to find keyword descriptions that contain the phrase
    for s in CKdict.values():
        if phrase.lower() in s.get("Description"):
            # get the dictionary index
            k = list(CKdict.values()).index(s)
            # put the corresponding keywords into a candidate list
            keys.append(list(CKdict.keys())[k])
    # show the hints for all candidate keywords
    if len(keys) > 0:
        for this_key in keys:
            keyword_hints(this_key)
    else:
        msg = [
            Color.PURPLE,
            "no keyword description containing the phrase",
            phrase,
            "can be found.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)


def help(topic: Union[str, None] = None):
    """Provide assistance on finding information about Chemkin keywords."""
    """Provide assistance on finding information about Chemkin keywords.

    Parameters
    ----------
        topic: string
            the keyword topic of which additional hints are requested

    """
    #
    if topic is None:
        # general information about getting help
        msg = [
            "For detailed information about all Chemkin keywords and",
            "reactor models,\n",
            Color.SPACEx6,
            "use \"ansys.chemkin.core.help('manual')\".",
        ]
        Color.ckprint("info", msg)
        msg = [
            "For usage of the real-gas cubic EOS",
            "in mixture thermodynamic property calculation,\n",
            Color.SPACEx6,
            "use \"ansys.chemkin.core.help('real-gas')\".",
        ]
        Color.ckprint("normal", msg)
        msg = [
            "For mixture equilibrium calculation options,\n",
            Color.SPACEx6,
            "use \"ansys.chemkin.core.help('equilibrium')\".",
        ]
        Color.ckprint("normal", msg)
        msg = [
            "For information about a Chemkin reactor model keyword,\n",
            Color.SPACEx6,
            "use \"ansys.chemkin.core.help('keyword')\".",
        ]
        Color.ckprint("normal", msg)
        msg = [
            "For batch reactors ignition delay time definitions,\n",
            Color.SPACEx6,
            "use \"ansys.chemkin.core.help('ignition')\".",
        ]
        Color.ckprint("normal", msg)
    elif topic.lower() in "manual manuals":
        # information about chemkin manuals
        msg = [
            "ansys.chemkin.core.manuals will open the Chemkin manuals page",
            "of the Ansys Help portal",
            "in a new tab of the default browser.\n",
            Color.SPACEx6,
            "* provide Ansys login credentials to access the manuals.\n",
            Color.SPACEx6,
            "* for Chemkin keywords: check out the Input manual.\n",
            Color.SPACEx6,
            "* for Chemkin reactor models: check out the Theory manual.\n",
        ]
        Color.ckprint("info", msg)
        manuals()
    elif topic.lower() in "keyword keywords":
        # information about getting information about specific reactor model keywords
        msg = [
            "For information about a Chemkin reactor model keyword,\n",
            Color.SPACEx6,
            "use \"ansys.chemkin.core.keyhints('<keyword>')\"\n",
            Color.SPACEx6,
            "ex: ansys.chemkin.core.keyhints('HTC')\n",
        ]
        Color.ckprint("normal", msg)
        msg = [
            "For information about keywords related to a phrase,\n",
            Color.SPACEx6,
            "use \"ansys.chemkin.core.phrase_hints('<phrase>')\"\n",
            Color.SPACEx6,
            "ex: ansys.chemkin.core.phrase_hints('tolerance')\n",
        ]
        Color.ckprint("normal", msg)
    elif topic.lower() in "ignition":
        # ignition definition options
        show_ignition_definitions()
    elif topic.lower() in "real gas real-gas":
        # real-gas model usage
        show_realgas_usage()
    elif topic.lower() in "equilibrium":
        # equilibrium calculation options
        show_equilibrium_options()
    else:
        msg = [
            Color.PURPLE,
            "cannot find help topic:",
            topic,
            "\n",
            Color.SPACE,
            "the valid topics are\n",
            Color.SPACEx6,
            "'manual', 'keyword', 'ignition', 'real gas', and 'equilibrium'.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)


def show_realgas_usage():
    """Show Chemkin real-gas model usage and options."""
    print(
        Color.YELLOW
        + "** the real-gas cubic equation of state is available "
        + "only when the mechanism contains the 'EOS_' data"
    )
    print(
        "   after the gas-phase mechanism is pre-processed, "
        + "the pre-processor will indicate "
        + "if the mechanism contains the necessary real-gas data"
    )
    print("   * for real-gas eligible mechanisms,")
    print("     > using the real-gas EOS with mixtures:")
    print("       to check the current activation status of the real-gas EOS, use")
    print("           status = ansys.chemkin.core.check_realgas_status()")
    print("              status = True means the real-gas EOS is active")
    print("                     = False means the ideal gas law is active")
    print(
        "       to activate the real-gas cubic EOS in mixture property calculation, use"
    )
    print("              <mixture_object>.use_realgas_cubicEOS()")
    print("       to select the mixing rule after the real-gas EOS is activated, use")
    print("              <mixture_object>.set_realgas_mixing_rule(rule)")
    print("              rule is the mixing rule option:")
    print("                   0: the Van der Waals mixing method (default)")
    print("                   1: the pseudo-critical property method")
    print(
        "       to deactivate the real-gas cubic EOS in mixture property "
        "calculation, use"
    )
    print("              <mixture_object>.use_idealgas_law()")
    print("     > using the real-gas EOS with reactor models:")
    print("       see reactor model keywords: 'RLGAS' and 'RLMIX'")
    print("              ex: ansys.chemkin.core.keyword_hints('RLGA'", end=Color.END)


def show_equilibrium_options():
    """Show the equilibrium calculation usage and options."""
    print(Color.YELLOW + "** equilibrium calculation usage: ")
    print("      EQ_mixture = ansys.chemkin.core.equilibrium(INIT_mixture, opt)")
    print("      INIT_mixture is the initial mixture (object)")
    print("      EQ_mixture is the final/equilibrium mixture (object)")
    print("      opt is the equilibrium calculation option: ")
    print("           1: SPECIFIED T AND P (default)")
    print("           2: SPECIFIED T AND V")
    print("           4: SPECIFIED P AND V")
    print("           5: SPECIFIED P AND H")
    print("           7: SPECIFIED V AND U")
    print("           8: SPECIFIED V AND H")
    print("** Chapman-Jouguet detonation calculation usage:")
    print("      speed_list, CJ_mixture = ansys.chemkin.core.detonation(INIT_mixture)")
    print("      INIT_mixture is the initial mixture (object)")
    print("      CJ_mixture is the C-J state mixture (object)")
    print("      speed_list is a list consists of two speed values at the C-J state: ")
    print(
        "           [sound_speed, detonation_wave_speed] in cm/sec",
        end=Color.END,
    )


def show_ignition_definitions():
    """Show the ignition definitions available in Chemkin."""
    # show ignition definition usage
    print(Color.YELLOW + "** ignition definition is assigned as a string, e.g., 'OH' ")
    print("   valid options are: ")
    print("    1: 'T_inflection'")
    print("    2: 'T_rise', <val>")
    print("    3: 'T_ignition', <val>")
    print("    4: 'Species_peak', '<target species>'", end=Color.END)


def manuals():
    """Access the Chemkin manuals page on the Ansys Help portal."""
    # Chemkin manual page
    chemkin_manual_url = (
        "https://ansyshelp.ansys.com/account/secured?returnurl="
        "/Views/Secured/prod_page.html?" + "pn=Chemkin&pid=ChemkinPro&lang=en"
    )
    # open the web page on a new tab of the default web browser
    webbrowser.open_new_tab(chemkin_manual_url)

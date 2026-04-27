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

"""Surface properties test for the Chemkin surface utilities."""

from pathlib import Path

import numpy as np

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set numpy printing option
np.set_printoptions(legacy="1.25")

# Create a chemistry set with surface chemistry
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
# find the floder where this py example file is located
script_dir = Path(__file__).resolve().parent
# use the relative path to find the location of the example mechanism files
example_mech_data = script_dir / ".." / "data"
print(f"example data = {str(example_mech_data)}")
# create a chemistry set based on the gasoline 14-components mechanism
MySurfMech = ck.Chemistry(label="multi_materials")
# set mechanism input files
# including the full file path is recommended
local_data_dir = "CH4_MB89"
MySurfMech.chemfile = str(example_mech_data / local_data_dir / "chem_ch4_MB89.inp")
MySurfMech.surffile = str(example_mech_data / local_data_dir / "surf_ch4_MB89.inp")
MySurfMech.thermfile = str(data_dir / "therm.dat")
# preprocess the mechanism files
ierror = MySurfMech.preprocess()
if ierror != 0:
    print("Error: Failed to preprocess the mechanism!")
    print(f"       Error code = {ierror}")
    exit()

# Access surface chemistry information in the Chemistry Set
# the number of surface materials in the Chemistry Set
print(f"Total number of surface materials = {MySurfMech.number_materials}")
# get all surface material names in the Chemistry Set with surface chemistry
m_symbols = MySurfMech.material_names
# use the material name to access the corresponding ``Material`` object
# in the Chemistry Set
# get the names of all phases (gas phase and surface (site and bulk) phases)
for s in m_symbols:
    print(
        f"Material '{s.rstrip()}' "
        f"contains phases: {MySurfMech.materials[s].phase_names}\n"
    )

# Display surface material phase and species information as well as
# the surface reactions of the material.

# list the commonly used properties of all surface materials of the Chemistry Set
for n, name in enumerate(m_symbols):
    # get the surface Material object by its name
    m = MySurfMech.materials[name]
    print("=" * 80)
    print(f"\nMaterial index = {n}")
    # display the general size information of the material
    m.information()
    print("=" * 40)
    # list the elemental compositions of all site species of the material
    print("elemental compositions of surface species on the material:")
    if m.num_site_species > 0:
        for i, p in enumerate(m.site_species_names):
            print(f"  site species {p.rstrip()}")
            k = m.site_species_map[p]
            for e, ele in enumerate(MySurfMech.element_symbols):
                print(f"    element {ele} = {m.material_species_composition(e, k)}")
        h_site = m.get_site_species_h(temp=400.0)
        for k, h in enumerate(h_site):
            # local index among all site species of this material
            print(
                f"  site species {k} {m.site_species_names[k]} "
                f"H = {h / ck.ERGS_PER_JOULE} [J/mole]"
            )
    # list the elemental compositions of all bulk species of the material
    if m.num_bulk_species > 0:
        for i, p in enumerate(m.bulk_species_names):
            print(f"  bulk species {p.rstrip()}")
            k = m.bulk_species_map[p]
            for e, ele in enumerate(MySurfMech.element_symbols):
                print(f"    element {ele} = {m.material_species_composition(e, k)}")
        # also the specific heat capacity of the bulk species
        cp_bulk = m.get_bulk_species_cp(temp=450.0)
        for k, c in enumerate(cp_bulk):
            # local index among all bulk species of this material
            print(
                f"  bulk species {k} {m.bulk_species_names[k]} "
                f"Cp = {c / ck.ERGS_PER_JOULE} [J/mole-K]"
            )
    # get all phase names of the Chemistry Set, the "Gas" phase is always included
    # as the first phase
    p_symbols = m.phase_names
    # list the surface phase types and the index of the first species of the
    # surface phase.
    # There are three phase types: "Gas", "site", and "bulk". The surface species
    # (site and bulk) have two indices: global (of all species including the gas
    # species) and local (of the same type of surface species).
    for p in p_symbols:
        if p == "Gas":
            continue
        phase_type = m.phases[p].phase_type
        #
        if phase_type == "site":
            print(
                f" surface site phase: {p} "
                f"site density = {m.phases[p].site_density} [mole/cm2]\n"
            )
            # the global index is 1-based in Chemkin
            print(
                f" 1-base species index on phase {p} = "
                f"{m.phases[p].first_species_index}"
            )
            # the global index of the material, the phase, and the species
            # is 0-base in PyChemkin
            print(
                f" 0-base species index on site phase {p} "
                f"= {m.first_site_species_index}"
            )
        else:
            # the global index is 1-based in Chemkin
            print(
                f" 1-base species index on phase {p} = "
                f"{m.phases[p].first_species_index}"
            )
            # the global index of the material, the phase, and the species
            # is 0-base in PyChemkin
            print(
                f" 0-base species index on bulk phase {p} "
                f"= {m.first_bulk_species_index}"
            )

    # list surface reaction information of the material, the reaction index
    # is 1-based in both Chemkin and PyChemkin.
    print("=" * 40)
    # the first surface reaction of the material
    rxn_id = 1
    # the first surface reaction of the material
    rxn_id = 1
    # reaction string of the first surface reaction of the material
    print(f" surface reaction # 1 = {m.get_surface_reaction_string(rxn_id)}")
    # get the Arrhenius rate parameters of all surface reactions of the material
    a_factor, beta, act_energy = m.get_surface_reaction_parameters()
    # display the Arrhenius parameters of the first surface reaction of this material
    print("    original rate parameters:")
    print(f"    A = {a_factor[0]}   B = {beta[0]}  Ea = {act_energy[0]} [K]")
    # save a copy of the original value
    act_energy_org = act_energy[0]
    # change the activation temperature of the first surface reaction to 3000 [K]
    m.set_surface_reaction_act_energy(rxn_id, 3000.0)
    # retrieve and display the modified Arrhenius parameters
    # of the first surface reaction
    a_factor, beta, act_energy = m.get_surface_reaction_parameters()
    print("    modified rate parameters:")
    print(f"    A = {a_factor[0]}   B = {beta[0]}  Ea = {act_energy[0]} [K]")
    # restore the activation temperature of the first surface reaction
    m.set_surface_reaction_act_energy(rxn_id, act_energy_org)
    # verify the change
    a_factor, beta, act_energy = m.get_surface_reaction_parameters()
    print("    restored rate parameters:")
    print(f"    A = {a_factor[0]}   B = {beta[0]}  Ea = {act_energy[0]} [K]")

resultfile = Path(current_dir) / "multiple_materials.result"
results = {}
results["state-site_H"] = h_site.tolist()
results["state-bulk_Cp"] = cp_bulk.tolist()
results["state-A_factor"] = a_factor.tolist()
results["state-activation_temperature"] = act_energy.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

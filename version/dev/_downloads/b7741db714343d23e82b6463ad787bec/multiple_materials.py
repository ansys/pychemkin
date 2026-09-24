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

r""".. _ref_multiple_surface_materials:

======================================================
Processing a surface mechanism with multiple materials
======================================================
When the surface mechanism input file is added to the *Chemistry Set* via the
``surffile`` parameter, the ``preprocess()`` method will process the
surface mechanism in the file. Typically, a surface mechanism consists of
one surface *material* block. The *material* block defines the surface phases,
the site and the bulk surface species of the surface phases, the thermodynamic
data of all surface (site and bulk) species, and the surface reactions describing
the chemical interactions between the gas-phase species and the surface species
and between the surface species of the surface material. The gas phase is always
the first phase of any surface material followed by all the site phases and then all
the bulk phases of the material. The species are sorted in two different ways:
the *global order* and the *local order*. The *global species index/order* combines
species from all phases and arranges them in the order of gas phase species,
species from all site phases of the material, followed by species from all bulk phases
of the material. The *local species index/order* lists the species from the same type
of surface phase. For example, local index of all site species of the surface material.
Here is a simple graph showing different species arrangements of a surface material:

.. figure:: global_species_mapping.png
  :scale: 60 %
  :alt: species index

When there are multiple surface materials defined in the surface mechanism input file,
all surface materials share the same gas phase (gas species and gas-phase reactions).
Each surface material consists of its own surface phases, surface species, and
surface reactions, and has independent *global* and *local* species orders/indices.
Subsequently, you can access species and reaction properties one surface material
at a time, and direct surface species interactions across different surface materials
are *not* supported.

For detailed information about the *Surface Chemkin*, you can check out the
*Chemkin theory* and the *Chemkin Input* manuals at the Ansys Product Help website.

When a ``Chemistry Set`` with surface chemistry is used to initiate a ``Mixture``
or a ``Stream`` object, Pychemkin will automatically set up the surface chemistry
sub-module of the ``Mixture`` or the ``Stream`` object. Similarly, the surface
chemistry and the surface governing equations will be ready when a ``ReactorModel``
is initiated with a surface chemistry enabled ``Mixture`` or ``Stream``. There is
no additional steps required to manually "activate" the surface chemistry in
Pychemkin objects other than the ``Chemistry Set``. The figure below shows the
initialization process of the "surface chemistry modules" in different Pychemkin
objects:

.. figure:: surface_chemistry_module.png
    :scale: 60 %
    :alt: surface chemistry module

This example demonstrates the procedures of preprocessing a ``Chemistry Set`` that
supports surface chemistry. It also shows how to access properties from two different
surface materials. Typically, you would want to get information about the symbols,
the index, the density, the molecular weights, and the elemental compositions of
the species. For site species, you can get the "site occupancy". For bulk species,
you can get the "enthalpy" and the "specific heat capacity".
"""

# sphinx_gallery_thumbnail_path = '_static/multiple_surface_materials.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from pathlib import Path

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

###############################################
# Create a chemistry set with surface chemistry
# =============================================
# This example uses mechanism files from the "examples\data" folder, and
# assumes that this Python example file is located in the
# "examples\surface_chemistry" folder.
# Please verify the location and the availability of these files and modify
# the directories below to point to the correct files before running
# this example.
#
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
# find the floder where this py example file is located
try:
    script_dir_obj = Path(__file__).parent.resolve()
except NameError:
    script_dir_obj = Path(current_dir)
script_dir = str(script_dir_obj)
# use the relative path to find the location of the example mechanism files
example_mech_data = script_dir_obj / ".." / "data"
print(f"example data = {str(example_mech_data)}")
# create a chemistry set based on the gasoline 14-components mechanism
MySurfMech = ck.Chemistry(label="multi_materials")
# set mechanism input files
# including the full file path is recommended
local_data_dir = "CH4_MB89"
MySurfMech.chemfile = str(example_mech_data / local_data_dir / "chem_ch4_MB89.inp")
MySurfMech.surffile = str(example_mech_data / local_data_dir / "surf_ch4_MB89.inp")
MySurfMech.thermfile = str(data_dir / "therm.dat")

#######################################
# Preprocess the surface chemistry set
# =====================================
# preprocess the mechanism files
_ = MySurfMech.preprocess()

###########################################################
# Access surface chemistry information in the Chemistry Set
# =========================================================
# Access the `Material` and the `SurfacePhase` objects in the `Chemistry Set`
# containing surface mechanism by their names given in the surface mechanism
# input file. For example, ``MySurfMech.materials['CLEANUP']`` for the
# ``Material`` object associated with the surface material 'CLEANUP' defined
# in the surface mechanism input file, or
# ``MySurfMech.materials['INITIATE'].phases['QUARTZ']`` for the
# ``SurfacePhase`` object associated with the 'QUARTZ' phase of the surface
# material 'INITIATE'.
#
# You can get the surface material names by using the
# ``MySurfMech.material_names`` method. The ``Material`` object provides
# information about the surface material such as the number of surface species
# and surface reactions, the surface phases, and the properties of the site
# species and the bulk species.
# Some of the available properties in the ``Material`` and the ``SurfacePhase``
# are shown below.

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

############################################################
# Information about the surface material, phase, and species
# ==========================================================
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
    # is 1-base in both Chemkin and PyChemkin.
    print("=" * 40)
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

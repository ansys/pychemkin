# Copyright (c) 2026 Synopsys, Inc. and ANSYS, Inc. and/or its affiliates.
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

"""Test for post-processing Chemkin XML solution from reactor models."""

from pathlib import Path

import matplotlib.pyplot as plt

from ansys.chemkin.core.chemkin_solution_utility import ChemkinSolutionImporter
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.utilities import (
    copy_file,
    delete_files_by_extension,
)

# check working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)

#######################################
# Provide the Chemkin XML solution file
# =====================================
# Find the data folder containing the XML solution file.

try:
    script_dir_obj = Path(__file__).parent.resolve()
except NameError:
    script_dir_obj = Path(current_dir)
script_dir = str(script_dir_obj)
# use the relative path to locate the XML solution data folder
local_data_folder = script_dir_obj / ".." / "data" / "1_D_solution"
print(f"solution data folder = {str(local_data_folder)}")
# set the xml solution data file
xml_solution_file = "XMLdata_opposed-flow_flame__gas_soot_radiation.zip"
xml_solution_path = local_data_folder / xml_solution_file

##########################################
# Initialize the Chemkin solution importer
# ========================================
# Initialize the Chemkin solution importer which will store all the
# raw solution data read from the XML solution file. Specify the working directory,
# clean up any existing solution files there before importing new solution data.
# Copy the XML solution file from the data directory to the working directory, and
# provide the full path of the XML file to the importer.
# The ``ChemkinSolutionImporter`` have a few options:
#
# 1. use ``get_help()`` method to get information about the Chemkin post-processor
#    ``GetSolution``
# 2. use ``generate_preference_file()`` method to create a preference file.
#    You can make changes to the preference file ``CKSolnList.txt`` to customize
#    the filters and the units before running the importer. Save a copy of
#    the preference file to a different location to avoid it being deleted during
#    the cleanup step.
# 3. Once the customized preference file is ready, use the
#    ``use_preference_file()`` method to apply it before reading the XML solution file.
#
# Then use the ``import_solution()`` and the ``process_solution_data()`` methods to
# import the XML solution file and to process the solution data for further processing.
#

# Instantiate the Chemkin solution importer object
ck_solution = ChemkinSolutionImporter()
# set the working directory
ck_solution.set_working_dir(current_dir)
# clean up the working directory before copying the new solution file
delete_files_by_extension(
    targetpath=current_dir, extensions=[".ckcsv", ".txt", ".xml", ".zip"]
)
# copy the solution file to the current working directory
copy_file(str(local_data_folder), current_dir, xml_solution_file)
# set the XML solution data file for post-processing
ck_solution.set_xml_data_file(str(Path(current_dir) / xml_solution_file))
# - optional step: generate the preference file from the solution data
# the preference file name is "CKSolnList.txt".
# ck_solution.generate_preference_file(str(Path(current_dir) / xml_solution_file))
#
# - optional step: edit the preference file before generating the solution data
# pause here if you want to edit the preference file before applying it
# Remember to save the edited preference file to a different location,
# so that it won't be deleted by the cleanup step.
#
# - optional step: apply the modified preference file
# copy it to the current working directory if it is stored elsewhere
# copy_file(str(local_data_folder), current_dir, "CKSolnList.txt")
# apply the preference file
# ck_solution.use_preference_file(str(Path(current_dir) / "CKSolnList.txt"))
#
# import solution data from the XML file
ck_solution.import_solution()
# convert the solution data into group objects for further analysis
ck_solution.process_solution_data()

#########################################
# Post-process the imported solution data
# =======================================
# Now you can get basic information about the solution data such as
# the number of solution groups, the number of solution points,
# the number of solution variables, etc.
# You use the methods provided by the ``ChemKinSolution`` object and/or the
# ``SolutionGroup`` object to manipulate the solution data for furrther analysis
# and plotting.

# get the list of solution group names from the imported solution data
solution_groups = ck_solution.get_solution_groups()
# print basic information about each solution group
for group in solution_groups:
    print(f"Solution Group: {group}")
    print(f"number of solution points: {ck_solution.get_number_of_points()}")
    print(f"Number of gas species: {ck_solution.get_number_of_gas_species()}")
    # print(f"variables: {ck_solution.get_variable_labels(group)}")

# print(f"species symbols: {ck_solution.get_gas_species_symbols()}")

# Instantiate a solution group data object from the first (and the only) solution group
# in the solution data.
grp_idx = 0
# get the name of the first solution group
this_group = solution_groups[grp_idx]
# get the solution group object for the first solution group
this_group_obj = ck_solution.get_group_object(this_group)
# define the independent variable
independent_var = "Distance"
# get the array of the independent variable and its unit from the solution group object
distance = this_group_obj.get_variable_array(independent_var)
distance_unit = this_group_obj.get_variable_unit(independent_var)
# define the solution variable to be plotted
soln_var_1 = "Particle_volume_fraction"
soln_profile_1 = this_group_obj.get_variable_array(soln_var_1)
if len(soln_profile_1) <= 0:
    print(f"Error...No data available for solution variable: {soln_var_1}")
    exit()
soln_unit_1 = this_group_obj.get_variable_unit(soln_var_1)
# soln_var_2 = "Temperature"
# soln_var_2 = "Axial_velocity"
soln_var_2 = "A4"
soln_profile_2 = this_group_obj.get_variable_array(soln_var_2)
if len(soln_profile_2) <= 0:
    print(f"Error...No data available for solution variable: {soln_var_2}")
    exit()
soln_unit_2 = this_group_obj.get_variable_unit(soln_var_2)

#####################################
# Plot the selected solution variable
# ===================================
# plot solution profile
line_colors = ["b", "r"]
fig, ax1 = plt.subplots()
# Plot the primary data on the left y-axis
color1 = line_colors[0]
ax1.plot(distance, soln_profile_1, color=color1, label=soln_var_1.replace("_", " "))
ax1.set_xlabel(f"{independent_var} {distance_unit}")
ax1.set_ylabel(
    f"{soln_var_1.replace('_', ' ')} {soln_unit_1.replace('_', ' ')}", color=color1
)
ax1.tick_params(axis="y", labelcolor=color1)

# Create the secondary axes sharing the x-axis
ax2 = ax1.twinx()

color2 = line_colors[1]
ax2.plot(distance, soln_profile_2, color=color2, label=soln_var_2.replace("_", " "))
ax2.set_xlabel(f"{independent_var} {distance_unit}")
ax2.set_ylabel(
    f"{soln_var_2.replace('_', ' ')} {soln_unit_2.replace('_', ' ')}", color=color2
)
ax2.tick_params(axis="y", labelcolor=color2)
# plot results
plt.savefig("plot_post_process_flame.png", bbox_inches="tight")

# return results for comparisons
resultfile = Path(current_dir) / "flame_post_process.results"
results = {}
results["state-distance"] = distance.to_numpy(dtype=float, copy=False)
results["species-volume_fraction"] = soln_profile_1.to_numpy(dtype=float, copy=False)
results["species-A4_mole_fraction"] = soln_profile_2.to_numpy(dtype=float, copy=False)
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

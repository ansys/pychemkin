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

r""".. _ref_post_process_cvd_flow:

==============================================
Post-process Chemkin shear-layer flow solution
==============================================

PyChemkin has a post-process module ``ChemkinSolutionImporter`` that can used
to import the XML solution data from selected Chemkin reactor models into Python.
This example shows how use some of objects and utilities from the
``ChemkinSolutionImporter``, the ``SolutionData``, and the ``SolutionGroup`` modules
to extract the solution variables from a planar shear-layer flow simulation.
In combination of other PyChemkin modules,
such as ``Chemistry`` and ``Mixture``, you can further analyze and visualize
the simulation results.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_post_process_cvd_flow.png'

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
xml_solution_file = "XMLdata_planar_shear_flow__tcs_cvd.zip"
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

# print(f"variables: {ck_solution.get_variable_labels(solution_groups[0])}")
# print(f"species symbols: {ck_solution.get_gas_species_symbols()}")

# Instantiate a solution group data objects from a list of solution group names
# (corresponding to different distances from the entrance) in the solution data.
group_names = [
    "slice_5.000000",
    "slice_10.000000",
    "slice_20.000000",
    "slice_29.800000",
]
# curve attributes
curvelist = ["g", "b--", "r:", "m-."]
#
soln_x = {}
soln_y = {}
soln_unit = {}
# loop through the selected solution slices (groups)
for this_group in group_names:
    # get the solution group object for the selected solution group
    this_group_obj = ck_solution.get_group_object(this_group)
    # define the independent variable
    independent_var = "Height"
    # get the array of the independent variable and its unit from
    # the solution group object
    distance_series = this_group_obj.get_variable_array(independent_var)
    distance = distance_series.to_numpy(dtype=float, copy=False)
    distance_unit = this_group_obj.get_variable_unit(independent_var)
    # define the solution variable to be plotted
    # soln_var = "Temperature"
    soln_var = "SI2CL6"
    # soln_var = "SICL4"
    soln_profile_series = this_group_obj.get_variable_array(soln_var)
    soln_profile = soln_profile_series.to_numpy(dtype=float, copy=False)
    if len(soln_profile) <= 0:
        print(f"Error...No data available for solution variable: {soln_var}")
        exit()
    soln_y[this_group] = distance
    soln_x[this_group] = soln_profile
    soln_unit[this_group] = this_group_obj.get_variable_unit(soln_var)

#####################################
# Plot the selected solution variable
# ===================================
# plot solution profiles at different distances from the entrance.
for i, grp in enumerate(group_names):
    plt.plot(
        soln_x[grp], soln_y[grp], curvelist[i], label=grp.replace("slice_", "x = ")
    )
plt.ylabel(f"{independent_var} {distance_unit}")
y_label = f"{soln_var} {soln_unit[group_names[0]].replace('_', ' ')}"
plt.xlabel(y_label)
plt.legend([grp.replace("slice_", "x = ") for grp in group_names], loc="upper right")
plt.title("Profile Evolution along Channel")
# plot results
plt.show()
# plt.savefig("plot_post_process_cvd_flow.png", bbox_inches="tight")

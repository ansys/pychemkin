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

r""".. _ref_run_chemkin_applications:

==================================
Run chemkin applications in Python
==================================

PyChemkin does not support all chemkin functionalities and reactor models. However,
it is still possible to run simulations of unsupported applications in Python.

This example describes the procedures to run any chemkin reactor model simulation in
Python environment.
By doing so, it opens up the possibility for interactions between chemkin applications
and other applications under the Python environment. Thus, this example is primarily
for Chemkin users
(1) who know how to generate the application-specific input file and
(2) who are familiar with the application keywords.

The entire process largely follows the same work flow for running any chemkin application
from the command line. The first step is to identify the os platform and the correct
location of the Ansys Chemkin installation.
The second step is to set up and preprocess the ``Chemistry Set``.
The third step is to specify the application text input file, text output file, and
the XML solution file. Then, create the batch script to run the simulation and
post-process the solutions.
The last step is to run the batch script created with the ``subprocess.run()`` method.
The "final" solution data should be in MS EXCEL csv format stored in a number of
"CKSoln_solution_xxx.csv" files. The original XML solution file should also be available.

The advantages of direct execution of the chemkin application
(without using the PyChemkin methods) are

1. all chemkin applications are available;
2. access to all application options;
3. can run chemkin applications with custom-built user subroutines;
4. complete set of solution (including the A-factor sensitivities) is created and
    can be further post-processed by various utilities.

On the other hand, the disadvantages are

1. problem setup is less interactive, you need to know how to manipulate the contents of
    a text file to change options/parameters in Python;
2. additional steps are required to post-process the solution files.

This example employs the chemkin rotating-disk reactor model for silicon deposition
from silane. You can easily modify this example for running other chemkin applications.
For details of the chemkin applications and their inputs/options, please check out the
"Theory" and the "Input" manuals at the Ansys Product Help website.
Instructions of running chemkin applications from the command line and the "get solution"
utilities can be found in the "Getting Started" and the "Visualization" manuals.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_cvd_profiles.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import pandas as pd
from pathlib import Path
import platform
import subprocess
import time

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.chemistry import chemkin_bin_dir
from ansys.chemkin.core import Color
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.utilities import copy_file, delete_files_by_extension
import matplotlib.pyplot as plt


# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True

###########################################
# Establish the chemkin runtime environment
# =========================================
# The batch script to establish the chemkin runtime environment is platform
# dependent. The first task is to find out the platform and the shell environment.
# Then you will run the environment setup script accordingly.

# check platform
under_win = False
plat = platform.system()
if plat == "Windows":
    # Windows
    under_win = True
    # set the chemkin runtime environment batch script
    env_file = "run_chemkin_env_setup.bat"
    # define the chemkin application
    ck_app = "CKReactorRotatingDiskCVD.exe"
    # post-processer
    get_soln = "GetSolution.exe"
    soln_transpose = "CKSolnTranspose.exe"
elif plat == "Linux":
    # Linux
    # set the chemkin runtime environment for the /bash shell script
    env_file = "chemkin_setup.ksh"
    # define the chemkin application
    ck_app = "CKReactorRotatingDiskCVD"
    # post-processer
    get_soln = "GetSolution"
    soln_transpose = "CKSolnTranspose"
else:
    msg = [
        Color.RED,
        "unsupported platform",
        str(platform.system()),
        "\n",
        "PyChemkin does not support the current os.",
        Color.END,
    ]
    this_msg = Color.SPACE.join(msg)
    logger.critical(this_msg)
    exit()

# run the chemkin steup script to set up the chemkin runtime environment
# define the chemkin bin folder that holds all the executables
ckbin = chemkin_bin_dir()
# set the full path to the chemkin environment script
run_env_path = str(Path(ckbin) / env_file)

#############################################################
# Set up all directories required to run chemkin applications
# ===========================================================
# Define the working directory ``current_dir`` used to run this Python script.
# Use the location of this Python source file ``script_dir`` to locate the data
# folder ``example_mech_sata`` containing the mechanism files and the application
# input files. And the data folder of the Ansys chemkin installation
# ``data_dir`` that stores common thermodynamic and the transport data files.
#
# .. note::
#       This example uses mechanism files and input files from the 'examples\data'
#       folder, and assumes that this Python example file is located in the
#       'examples\advanced' folder. Please verify the location and the availability
#       of these files and modify the directories below to point to the correct files
#       before running this example.
#

# get the working directory
current_dir_obj = Path.cwd()
current_dir = str(current_dir_obj)
logger.debug("working directory: " + current_dir)
# get the current py file location
try:
    script_dir_obj = Path(__file__).parent.resolve()
except NameError:
    script_dir_obj = current_dir_obj / ".." / "data"
script_dir = str(script_dir_obj)
logger.debug("script directory: " + script_dir)
# use the relative path to find the location of the example mechanism files
example_mech_data = script_dir_obj / ".." / "data"
logger.debug("example input directory: " + str(example_mech_data))
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
logger.debug("database directory: " + str(data_dir))

################################
# Instantiate the Chemistry Set
# ==============================
# Define the full path to the mechanism input files. This CVD example requires
# the gas mechanism, the surface mechanism, the thermodynamic data, and
# the transport data files. These files can be in different directories so make sure
# the correct directory path is used. For example, the mechanism files are in the
# 'examples\data' folder, but the data files are in the 'data' folder under the Ansys
# Chemkin installation.

# set up all the mechanism input files for preprocessing
local_data_dir = "SICH4_CVD"
chem_input = str(example_mech_data / local_data_dir / "sih4_cvd_chem.inp")
surf_input = str(example_mech_data / local_data_dir / "sih4_cvd_surf.inp")
therm_input = str(data_dir / "therm.dat")
tran_input = str(data_dir / "tran.dat")
# create the Chemistry Set
MyCVDMech = ck.Chemistry(
    chem=chem_input,
    surf=surf_input,
    therm=therm_input,
    tran=tran_input,
    label="SiH4_CVD",
)

#######################################
# Preprocess the surface chemistry set
# =====================================
# Once the ``Chemistry Set`` is created, run the preprocessor. The *asc files will
# have their default names, for instance, 'chem.asc', and they will sit in the
# current working folder.

# preprocess the mechanism files
_ = MyCVDMech.preprocess()

######################################
# Prepare the application input files
# ====================================
# Clean up the working folder by removing all application solution file.
# Copy the required input files from the 'example/data' folder to the working folder.

# delete the existing solution file *.zip
delete_files_by_extension(targetpath=current_dir, extensions=[".zip", ".dtd", ".bat"])
# copy the chemkindata.dtd from the data directory
copy_file(
    sourcepath=str(example_mech_data),
    targetpath=current_dir,
    filename="chemkindata.dtd",
)

# define full path to the chemkin application
run_file = Path(ckbin) / ck_app
# check the application
if not run_file.is_file():
    print(f"Error...the Chemkin application {str(run_file)} does not exist.")
    print("   please make sure the path is correct...")
    exit()

# define the application text input and output files
app_input = "rotating_disk__sih4_cvd.inp"
app_output = "rotating_disk__sih4_cvd.out"
# define the full path to the text files
# the input file must be prepared before running the application
# the input file for this example is already prepared and is stored in the
# 'examples\data' folder
in_file_path = example_mech_data / local_data_dir / app_input
# check the input file
if not in_file_path.is_file():
    print(f"Error...the application input file {str(in_file_path)} cannot be found.")
    print("   please make sure the input file exists and the path is correct...")
    exit()
# the text output file will be saved in the current working folder during the simulation
out_file_path = str(current_dir_obj / app_output)

# XML solution file
xml_file = str(current_dir_obj / "XMLdata.zip")

############################################################
# Create the run batch script to run the chemkin application
# ==========================================================
# The Chemkin applications and the ``GetSolution`` utilities rely on the
# runtime environment established by the ``run_env_path`` script. Therefore,
# all the simulation runs and the post-processing task involving the ``GetSolution``
# utilities must be included in the same batch script created below. In other words,
# you must establish the chemkin runtime environment in a batch script before running
# any chemkin applications. You can perform other non-chemkin related tasks after the
# ``subprocess.run()`` step.

if under_win:
    # Windows
    # define the run script name at the current working floder
    bat_file = current_dir_obj / "RUN_CKjob.bat"
    # create the batch script on the fly
    logger.debug("creating the batch script...")
    with bat_file.open(mode="w") as file:
        file.write("@REM Generated Batch Job Command File\n")
        file.write("REM get the runtime environment defined\n")
        file.write(f"@CALL \"{run_env_path}\"\n")
        file.write("@SET MKL_NUM_THREADS=1\n")
        file.write("@SET OMP_NUM_THREADS=1\n")
        file.write("@SET XERCES_DISABLE_DTD=1\n")
        file.write("REM RUN THE JOB AND CAPTURE EXIT STATUS\n")
        file.write(f"{ck_app} -i {in_file_path} -o {out_file_path} -x {xml_file}\n")
        file.write(r"@ECHO %errorlevel%" + "\n")
        file.write("@REM subprocess will terminate the script\n")
        file.write("@REM when it detects an error during the simulation\n")
        file.write("@ECHO running 'GetSolution'...\n")
        file.write(f"{get_soln} {xml_file}\n")
        file.write(f"{soln_transpose}\n")
        file.write("exit\n")
        file.close()
    # use subprocess.run() to run the script just created
    # set the start wall time
    start_time = time.time()
    logger.debug("simulation started...")
    try:
        result = subprocess.run([bat_file], capture_output=True, text=True, check=True)
        logger.debug("Command output:")
        logger.debug(result.stdout)
    except subprocess.CalledProcessError as e:
        logger.debug(f"Error executing command: {e.cmd} "
                     f"with non-zero exit status {e.returncode}.")
        logger.debug(f"Stderr: {e.stderr}")
        exit()
    except subprocess.TimeoutExpired as e:
        logger.debug(f"Command timed out: {e.timeout}")
        exit()
    except FileNotFoundError as e:
        logger.debug(f"Command not found: {e}")
        exit()
    except OSError as e:
        logger.debug(f"OS error occurred: {e}")
        exit()
else:
    # Linux
    # define the run script name at the current working floder
    bat_file = current_dir_obj / "RUN_CKjob.sh"
    # create the bash script on the fly
    logger.debug("creating the batch script...")
    with bat_file.open(mode="w") as file:
        file.write("#!/bin/bash -v\n")
        file.write("### Generated Batch Job Script File\n")
        file.write("### get the runtime environment defined\n")
        file.write(f". {run_env_path}\n")
        file.write("export MKL_NUM_THREADS=1\n")
        file.write("export OMP_NUM_THREADS=1\n")
        file.write("export XERCES_DISABLE_DTD=1\n")
        file.write("### RUN THE JOB AND CAPTURE EXIT STATUS\n")
        file.write(f"{ck_app} -i {in_file_path} -o {out_file_path} -x {xml_file}\n")
        file.write("### subprocess will terminate the script\n")
        file.write("### when it detects an error during the simulation\n")
        file.write(r"@echo $?" + "\n")
        file.write("@echo running 'GetSolution'...\n")
        file.write(f"{get_soln} {xml_file}\n")
        file.write(f"{soln_transpose}\n")
        file.write("exit\n")
        file.close()
    # use subprocess.run() to run the script just created
    # set the start wall time
    start_time = time.time()
    logger.debug("simulation started...")
    try:
        result = subprocess.run([bat_file], capture_output=True, text=True, check=True)
        logger.debug("Command output:")
        logger.debug(result.stdout)
    except subprocess.CalledProcessError as e:
        logger.debug(f"Error executing command: {e.cmd} "
                     f"with non-zero exit status {e.returncode}.")
        logger.debug(f"Stderr: {e.stderr}")
        exit()
    except subprocess.TimeoutExpired as e:
        logger.debug(f"Command timed out: {e.timeout}")
        exit()
    except FileNotFoundError as e:
        logger.debug(f"Command not found: {e}")
        exit()
    except OSError as e:
        logger.debug(f"OS error occurred: {e}")
        exit()

# compute the total runtime
runtime = time.time() - start_time
logger.debug("simulation done...")
print(f"Total simulation duration: {runtime} [sec]\n")

#############################
# Post process the solution
# ===========================
# Once the simulation is completed successfully, you have two options to post-process
# the solution:
# 
# 1. you can use open the chemkin GUI and use the ``visualizer``
#    from the top banner. In this case, you can stop at this point. You only need the
#    XML solution file located in the working folder.
# 2. you can use the chemkin ``GetSolution`` utilities to process the raw solution. The
#    steps will be described below.
#

if interactive:
    user_input = input(
        ">>> Press [Enter] to continue, or type 'x' and "
        "press [Enter] to terminate the example: "
    )
    if user_input.lower() == "x":
        logger.debug("example terminated by user...")
        exit() 

#############################
# Load and parse the solution
# ===========================
# Load the solution from the csv files generated by the chemkin utilities
# ``GetSolution`` and ``CKSolnTranspose``. You can find the raw solution
# file ``xml_file`` and processed solution files *.csv in the working folder.
# Use the ``pandas.read_csv()`` method to load the csv file in interest and
# extract the solution variables by their header names. Of course, if you are
# running this example in Windows, you can simply load the csv solution files
# into MS EXCEL.
#
# .. note::
#   The ``pandas`` package is required to run this part of the example.
#

logger.debug("loading csv solution file...")
# read the solutions from the csv file for post-processing
csv_solution = "CKSoln_solution_no_1.csv"
# Read the CSV file into a DataFrame
soln = pd.read_csv(csv_solution)
logger.debug("solution loading done...")
# get the solution variable labels from the loaded data
soln_variables = list(soln.keys())
# get solution profiles along the centerline of the rotating CVD reactor
# get the solution index of the variables to be plotted
dist_id = -1
vel_id = -1
temp_id = -1
sih4_id = -1
si2h6_id = -1
for i, p in enumerate(soln_variables):
    if "Distance" in p:
        # distance from the substrate/disk [cm]
        dist_id = i
    elif "Axial_velocity" in p:
        # axial velocity [cm/sec]
        vel_id = i
    elif "Temperature" in p:
        # gas temperature [K]
        temp_id = i
    elif "fraction_SIH4_" in p:
        # SiH4 mole fraction
        sih4_id = i
    elif "fraction_H2_" in p:
        # H2 mole fraction
        h2_id = i
# convert solution profiles into 1-D arrays
# distance from the surface along the reactor centerline [cm]
dist = soln[soln_variables[dist_id]]
# axial gas velocity [cm/sec]
ax_vel = soln[soln_variables[vel_id]]
# gas temperature [K]
temp = soln[soln_variables[temp_id]]
# SiH4 mole fraction [-]
mole_sih4 = soln[soln_variables[sih4_id]]
# H2 mole fraction [-]
mole_h2 = soln[soln_variables[h2_id]]

############################
# Plot the solution profiles
# ==========================
# Plot the species, the gas temperature, and the gas axial velocity profiles along the
# centerline of the CVD reactor.
#
#
plt.subplots(2, 2, sharey="row", figsize=(10, 8))
plt.suptitle("Rotating Disk CVD Reactor", fontsize=16)
plt.subplot(221)
plt.plot(temp, dist, "b-")
plt.xlabel("Gas Temperature [K]")
plt.ylabel("Distance From Surface [cm]")
plt.subplot(222)
plt.plot(-ax_vel, dist, "m-")
plt.xlabel("Gas Axial Velocity Towards Surface [cm/sec]")
plt.subplot(223)
plt.plot(mole_sih4, dist, "r-")
plt.xlabel("SiH4 Mole Fraction")
plt.ylabel("Distance From Surface [cm]")
plt.subplot(224)
plt.plot(mole_h2, dist, "g-")
plt.xlabel("H2 Mole Fraction")

# clean up
ck.done()

# plot results
if interactive:
    plt.show()
else:
    plt.savefig("plot_cvd_profiles.png", bbox_inches="tight")

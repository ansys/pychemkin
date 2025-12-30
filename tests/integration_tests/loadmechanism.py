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

"""Test for defining and preprocessing a chemistry set."""

#################################################
# Import PyChemkin package and start the logger
# ===============================================

from pathlib import Path

# import PyChemkin packages
import ansys.chemkin.core as ck
from ansys.chemkin.core import Color
from ansys.chemkin.core.logger import logger

# check the working directory
current_dir = str(Path.cwd())
logger.debug("working directory: " + current_dir)

# set PyChemkin verbose mode
ck.set_verbose(True)


#####################################
# Create a chemistry set
# ===================================
# The first mechanism to load is the GRI 3.0 mechanism for methane combustion.
# This mechanism and its associated data files come with the standard Ansys Chemkin
# installation under the subdirectory *"\reaction\data"*.

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir

# set mechanism input files
# including the full file path is recommended
# the gas-phase reaction mechanism file (GRI 3.0)
chemfile = str(mechanism_dir / "grimech30_chem.inp")
# the thermodynamic data file
thermfile = str(mechanism_dir / "grimech30_thermo.dat")
# the transport data file
tranfile = str(mechanism_dir / "grimech30_transport.dat")

# create a chemistry set instance based on the GRI 3.0 methane combustion mechanism
MyGasMech = ck.Chemistry(
    chem=chemfile,
    therm=thermfile,
    tran=tranfile,
    label="GRI 3.0",
)

###################################
# Pre-Process the ``Chemistry Set``
# =================================

# preprocess the mechanism files
ierror = MyGasMech.preprocess()

# display the pre-process status
print()
if ierror != 0:
    # When a non-zero value is returned from the process, check the text output files
    # chem.out, tran.out, or summary.out for potential error messages about
    # the mechanism data.
    print(f"Preprocessing error encountered. Code = {ierror:d}.")
    print(f"see the summary file {MyGasMech.summaryfile} for details")
    exit()
else:
    Color.ckprint("OK", ["PreProcess success", "!!!"])
    print("mechanism information:")
    print(f"number of elements = {MyGasMech.number_elements:d}")
    print(f"number of gas species = {MyGasMech.number_species:d}")
    print(f"number of gas reactions = {MyGasMech.number_gas_reactions:d}")

#####################################
# Display basic mechanism information
# ===================================

print(f"\nelement and species information of mechanism {MyGasMech.label}")

print("=" * 50)

# extract element symbols as a list
elelist = MyGasMech.element_symbols
# get atomic masses as numpy 1D double array
awt = MyGasMech.atomic_weight
# print element information
for k in range(len(elelist)):
    print(f"element # {k + 1:3d}: {elelist[k]:16} mass = {awt[k]:f}")

print("=" * 50)

# extract gas species symbols as a list
specieslist = MyGasMech.species_symbols
# get species molecular masses as numpy 1D double array
wt = MyGasMech.species_molar_weight
# print gas species information
for k in range(len(specieslist)):
    print(f"species # {k + 1:3d}: {specieslist[k]:16} mass = {wt[k]:f}")
print("=" * 50)


##################################################################
# Create the second ``Chemistry Set`` instance in the same project
# ================================================================
# The second mechanism to be loaded into the project is the C2-NOx mechanism
# for the combustion of C1-C2 hydrocarbons. This mechanism differs from
# the GRI mechanism in the sense that it is self-contained, that is,
# the thermodynamic and the transport data of all species are included
# in the mechanism input file C2_NOx_SRK.inp which sits in the same
# reaction data directory as the GRI mechanism input files.
#
# The same steps to instantiate the first ``Chemistry Set`` object can be
# applied to process any set of reaction mechanism files. Here,
# slightly different procedures are employed to instantiate the second
# ``Chemistry Set`` object. The object is first created, and
# the reaction mechanism files are specified one by one afterwards.
# The reaction mechanism file in this case contains all the necessary
# thermodynamic and transport data, therefore no need to specify
# the therm and the tran data files. Note that, an additional step is
# required to instruct the pre-processor to include the transport data.


# set the 2nd mechanism directory (the default Chemkin mechanism data directory)
mechanism_dir = data_dir

# create a chemistry set based on C2_NOx using an alternative method
My2ndMech = ck.Chemistry(label="C2 NOx")

# set mechanism input files individually
# this mechanism file contains all the necessary thermodynamic and transport data
# therefore no need to specify the therm and the tran data files
My2ndMech.chemfile = str(mechanism_dir / "C2_NOx_SRK.inp")

# instruct the preprocessor to include the transport properties
# only when the mechanism file contains all the transport data
My2ndMech.preprocess_transportdata()


##########################################
# Pre-Process the second ``Chemistry Set``
# ========================================

# preprocess the 2nd mechanism files
ierror = My2ndMech.preprocess()

# display the pre-process status
print()
if ierror != 0:
    # When a non-zero value is returned from the process, check the text output files
    # chem.out, tran.out, or summary.out for potential error messages about
    # the mechanism data.
    print(f"Preprocessing error encountered. Code = {ierror:d}.")
    print(f"see the summary file {My2ndMech.summaryfile} for details")
    exit()
else:
    Color.ckprint("OK", ["PreProcess success", "!!!"])
    print("mechanism information:")
    print(f"number of elements = {My2ndMech.mm:d}")
    print(f"number of gas species = {My2ndMech.kk:d}")
    print(f"number of gas reactions = {My2ndMech.ii_gas:d}")


#####################################
# Display basic mechanism information
# ===================================

print(f"\nelement and species information of mechanism {My2ndMech.label}")
print("=" * 50)

# extract element symbols as a list
elelist = My2ndMech.element_symbols
# get atomic masses as numpy 1D double array
awt = My2ndMech.awt
# print element information
for k in range(len(elelist)):
    print(f"element # {k + 1:3d}: {elelist[k]:16} mass = {awt[k]:f}")

print("=" * 50)

# extract gas species symbols as a list
specieslist = My2ndMech.species_symbols
# get species molecular masses as numpy 1D double array
wt = My2ndMech.wt
# print gas species information
for k in range(len(specieslist)):
    print(f"species # {k + 1:3d}: {specieslist[k]:16} mass = {wt[k]:f}")

print("=" * 50)
# return results for comparisons
resultfile = Path(current_dir) / "loadmechanism.result"
results = {}
results["state-AWT"] = awt.tolist()
results["state-WT"] = wt.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

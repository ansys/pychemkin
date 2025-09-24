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
import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin.mixture import Mixture, species_diffusion_velocity
from ansys.chemkin.logger import logger

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# including the full file path is recommended
chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

# preprocess the mechanism files
iError = MyGasMech.preprocess()
if iError != 0:
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {iError}")
    exit()

# get species sybols
specieslist = MyGasMech.species_symbols

# define grid points
xj = 0.1
xjp1 = 0.2
dX = xjp1 - xj

# create mixture at grid point J
mixtureJ = Mixture(MyGasMech)
mixtureJ.X = ck.Air.X()
mixtureJ.temperature = 900.0
mixtureJ.pressure = ck.Patm
mixtureJ.volume = 1.0

# create mixture at grid point J + 1
mixtureJP1 = Mixture(MyGasMech)
mixtureJP1.X = [("H2", 0.2), ("CO2", 0.4), ("CH4", 0.4)]
mixtureJP1.temperature = 400.0
mixtureJP1.pressure = ck.Patm
mixtureJP1.volume = 1.0

# compute species diffusion velocities from J to J + 1
YV = species_diffusion_velocity(mixtureJP1, mixtureJ, mode="mix", tdiff=False)
print("\nspecies mass diffusion velocity without thermal diffusivity:\n")
k = 0
for v in YV:
    if mixtureJ.X[k] + mixtureJP1.X[k] > 1.0e-8:
        print(f"species index = {k+1} {specieslist[k]}: diffusion velocity = {-v/dX} [cm/sec]")
    k += 1

# compute species diffusion velocities with thermal diffusivity from J to J + 1
print("\nspecies mass diffusion velocity with thermal diffusivity:\n")
YV_withTDIFF = species_diffusion_velocity(mixtureJP1, mixtureJ, mode="mix", tdiff=True)
k = 0
for v in YV_withTDIFF:
    if mixtureJ.X[k] + mixtureJP1.X[k] > 1.0e-8:
        print(f"species index = {k+1} {specieslist[k]}: diffusion velocity = {-v/dX} [cm/sec]")
    k += 1

# return results for comparisons
resultfile = os.path.join(current_dir, "diffusionvelocity.result")
results = {}
results["state-diffusion_velocity"] = YV.tolist()
results["state-diffusion_velocity"] = YV_withTDIFF.tolist()
#
r = open(resultfile, "w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

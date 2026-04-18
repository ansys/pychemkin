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

"""Transport utility test for computing species diffusion velocities."""

from pathlib import Path

import numpy as np

import ansys.chemkin.core as ck  # Chemkin
from ansys.chemkin.core.logger import logger
from ansys.chemkin.core.mixture import Mixture, species_diffusion_velocity

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

# set mechanism directory (the default Chemkin mechanism data directory)
data_dir = Path(ck.ansys_dir) / "reaction" / "data"
mechanism_dir = data_dir
# including the full file path is recommended
chemfile = str(mechanism_dir / "grimech30_chem.inp")
thermfile = str(mechanism_dir / "grimech30_thermo.dat")
tranfile = str(mechanism_dir / "grimech30_transport.dat")
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0")

# preprocess the mechanism files
ierror = MyGasMech.preprocess()
if ierror != 0:
    print("Error: failed to preprocess the mechanism!")
    print(f"       error code = {ierror}")
    exit()

# get species symbols
specieslist = MyGasMech.species_symbols

# define grid points
xj = 0.1
xjp1 = 0.2
dx = xjp1 - xj

# create mixture at grid point J
mixture_j = Mixture(MyGasMech)
mixture_j.x = ck.Air.x()
mixture_j.temperature = 900.0
mixture_j.pressure = ck.P_ATM
mixture_j.volume = 1.0

# create mixture at grid point J + 1
mixture_jp1 = Mixture(MyGasMech)
mixture_jp1.x = [("H2", 0.2), ("CO2", 0.4), ("CH4", 0.4)]
mixture_jp1.temperature = 400.0
mixture_jp1.pressure = ck.P_ATM
mixture_jp1.volume = 1.0

# compute species diffusion velocities from J to J + 1
v_y = species_diffusion_velocity(mixture_jp1, mixture_j, mode="mix", tdiff=False)
print("\nspecies mass diffusion velocity without thermal diffusivity:\n")
k = 0
for v in v_y:
    if mixture_j.x[k] + mixture_jp1.x[k] > 1.0e-8:
        print(
            f"species index = {k + 1} {specieslist[k]}: "
            f"diffusion velocity = {-v / dx} [cm/sec]"
        )
    k += 1

# compute species diffusion velocities with thermal diffusivity from J to J + 1
print("\nspecies mass diffusion velocity with thermal diffusivity:\n")
v_y_with_tdiff = species_diffusion_velocity(
    mixture_jp1, mixture_j, mode="mix", tdiff=True
)
k = 0
for v in v_y_with_tdiff:
    if mixture_j.x[k] + mixture_jp1.x[k] > 1.0e-8:
        print(
            f"species index = {k + 1} {specieslist[k]}: "
            f"diffusion velocity = {-v / dx} [cm/sec]"
        )
    k += 1

# return results for comparisons
resultfile = Path(current_dir) / "diffusionvelocity.result"
results = {}
results["state-diffusion_velocity"] = v_y.tolist()
results["state-diffusion_velocity"] = v_y_with_tdiff.tolist()
#
r = resultfile.open(mode="w")
r.write("{\n")
for k, v in results.items():
    r.write(f'"{k}": {v},\n')
r.write("}\n")
r.close()

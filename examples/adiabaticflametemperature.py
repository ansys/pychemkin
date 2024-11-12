import os

import chemkin as ck  # Chemkin
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set verbose mode
ck.setverbose(True)

# This is a pychemkin equivalent of equil_test07

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the diesel 14 components mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")

iError = MyGasMech.preprocess()

oxid = ck.Mixture(MyGasMech)
# set mass fraction
oxid.X = [("O2", 1.0)]
oxid.temperature = 295.15
oxid.pressure = ck.Patm  # 1 atm
# mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric?

fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.X = [("CH4", 1.0)]
fuel.temperature = oxid.temperature
fuel.pressure = oxid.pressure

mixture = ck.Mixture(MyGasMech)
mixture.pressure = oxid.pressure
mixture.temperature = oxid.temperature
products = ["CO2", "H2O"]

points = 12
deq = 0.1
equiv_ini = 0.5

T = np.zeros(points, dtype=np.double)
equiv = np.zeros_like(T, dtype=np.double)

add_frac = np.zeros(MyGasMech.KK, dtype=np.double)
iError = 0
for i in range(points):
    equiv_current = equiv_ini
    iError = mixture.XbyEquivalenceRatio(
        MyGasMech, fuel.X, oxid.X, add_frac, products, equivalenceratio=equiv_current
    )
    if iError != 0:
        raise RuntimeError
    result = ck.equilibrium(mixture, opt=5)
    T[i] = result.temperature
    equiv[i] = equiv_current
    equiv_ini = equiv_ini + deq
print(T)

plt.plot(equiv, T, "bs--")
plt.xlabel("Equivalence ratio")
plt.ylabel("Temperature [K]")
plt.show()

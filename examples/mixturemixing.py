import os

import chemkin as ck  # Chemkin

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set verbose mode
ck.setverbose(True)
# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on GRI 3.0
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# transport data not needed
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# create the fuel mixture
# note: mixture pressures are not specified because pressure is not required for the calculations here
# the mixing process is assumed to take place at fixed pressure; i.e., the mixtures are at the same pressure
fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.X = [("CH4", 1.0)]
fuel.temperature = 300.0
# create the air mixture
air = ck.Mixture(MyGasMech)
# set mole fraction
air.X = [("O2", 0.21), ("N2", 0.79)]
air.temperature = 300.0
# mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric?)
mixture_recipe = [(fuel, 1.0), (air, 17.19)]
# create the new mixture (the air-fuel ratio is by mass)
premixed = ck.isothermalmixing(
    recipe=mixture_recipe, mode="mass", finaltemperature=300.0
)
# list molar composition
premixed.listcomposition(mode="mole")
print()
# now create an argon mixture
ar = ck.Mixture(MyGasMech)
# species composition
ar.X = [("AR", 1.0)]
# mixture temperature
ar.temperature = 600.0
# dilute the premixed mixture adiabatically with the ar mixture by 30% by volume
dilute_recipe = [(premixed, 0.7), (ar, 0.3)]
# create the diluted mixture
diluted = ck.adiabaticmixing(recipe=dilute_recipe, mode="mole")
# list molar composition
diluted.listcomposition(mode="mole")
# show the mixture temperatures
print(f"the diluted mixture temperature is  {diluted.temperature:f} [K]")
print(f"the ar mixture temperature is       {ar.temperature:f} [K]")
print(f"the premixed mixture temperature is {premixed.temperature:f} [K]")

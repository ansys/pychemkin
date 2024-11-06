import os

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin batch reactor models (transient)
from chemkin.batchreactors.batchreactor import (
    GivenPressureBatchReactor_EnergyConservation,
)

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set verbose mode
ck.setverbose(False)

# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the diesel 14 components mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
# pre-process
iError = MyGasMech.preprocess()
if iError == 0:
    print("mechanism information:")
    print(f"number of gas species = {MyGasMech.KK:d}")
    print(f"number of gas reactions = {MyGasMech.IIGas:d}")
else:
    exit()
#
# create air-fuel mixture with equivalence ratio = 1.1
oxid = ck.Mixture(MyGasMech)
# set mole fraction
oxid.X = [("O2", 1.0), ("N2", 3.76)]
oxid.temperature = 900
oxid.pressure = ck.Patm  # 1 atm
fuel = ck.Mixture(MyGasMech)
# set mole fraction
fuel.X = [("C3H8", 0.1), ("CH4", 0.8), ("H2", 0.1)]
fuel.temperature = oxid.temperature
fuel.pressure = oxid.pressure
mixture = ck.Mixture(MyGasMech)
mixture.pressure = oxid.pressure
mixture.temperature = oxid.temperature
products = ["CO2", "H2O", "N2"]
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)
# create the air-fuel mixture by using the equivalence ratio method
iError = mixture.XbyEquivalenceRatio(
    MyGasMech, fuel.X, oxid.X, add_frac, products, equivalenceratio=1.1
)
if iError != 0:
    raise RuntimeError
if ck.verbose():
    mixture.listcomposition(mode="mole")
# get the original rate parameters
Afactors = np.zeros(MyGasMech.IIGas, dtype=np.double)
Beta = np.zeros_like(Afactors, dtype=np.double)
ActiveEnergy = np.zeros_like(Afactors, dtype=np.double)
Afactors, Beta, ActiveEnergy = MyGasMech.getreactionparameters()
if ck.verbose():
    for i in range(MyGasMech.IIGas):
        print(f"reaction: {i + 1}")
        print(f"A = {Afactors[i]}")
        print(f"Ea = {ActiveEnergy[i]}\n")
# compute the ignition delay time
#
# create a constant pressure batch reactor (with energy equation)
#
MyCONP = GivenPressureBatchReactor_EnergyConservation(mixture, label="CONP")
# show initial gas composition inside the reactor
MyCONP.listcomposition(mode="mole")
# set other reactor parameters
# reactor volume [cm3]
MyCONP.volume = 10.0
# simulation end time [sec]
MyCONP.time = 2.0
# turn ON adaptive solution saving
MyCONP.adaptivesolutionsaving(mode=True, steps=20)
# set tolerance
MyCONP.settolerances(absolute_tolerance=1.0e-10, relative_tolerance=1.0e-8)
# set ignition delay
# ck.showignitiondefinition()
MyCONP.setignitiondelay(method="T_inflection")
# run the nominal case
runstatus = MyCONP.run()
#
if runstatus == 0:
    # get ignition delay time
    delaytime_org = MyCONP.getignitiondelay()
    print(f"ignition delay time = {delaytime_org} [msec]")
else:
    # if get this, most likely the END time is too short
    print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    print("failed to find the ignition delay time of the nominal case")
    exit()
# brute force A factor sensitivities of ignition delay time
IGsen = np.zeros_like(Afactors, dtype=np.double)
perturb = 0.001  # increase by 0.1%
# loop over all reactions
for i in range(MyGasMech.IIGas):
    Anew = Afactors[i] * (1.0 + perturb)
    # update the A factor
    ireac = i + 1
    MyGasMech.setreactionAFactor(ireac, Anew)
    # run the reactor model
    runstatus = MyCONP.run()
    #
    if runstatus == 0:
        # get ignition delay time
        delaytime = MyCONP.getignitiondelay()
        print(f"ignition delay time = {delaytime} [msec]")
        # compute the normalized sensitivity coefficient
        IGsen[i] = (delaytime - delaytime_org) / perturb
        # restore the A factor
        MyGasMech.setreactionAFactor(ireac, Afactors[i])
    else:
        # if get this, most likely the END time is too short
        print(f"trouble finding ignition delay time for raection {ireac}")
        print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
        exit()

# print top ignition delay time sensitivity coefficients
top = 5
# rank the positive coefficients
posindex = np.argpartition(IGsen, -top)[-top:]
poscoeffs = IGsen[posindex]

# rank the negative coefficients
NegIGsen = -1.0 * IGsen
negindex = np.argpartition(NegIGsen, -top)[-top:]
negcoeffs = IGsen[negindex]
# print the top sensitivity coefficients
if ck.verbose():
    print("positive sensitivity coefficients")
    for i in range(top):
        print(f"reaction {posindex[i] + 1}: coefficient = {poscoeffs[i]}")
    print()
    print("negative sensitivity coefficients")
    for i in range(top):
        print(f"reaction {negindex[i] + 1}: coefficient = {negcoeffs[i]}")
#
# create a rate plot
plt.rcParams.update({"figure.autolayout": True})
plt.subplots(2, 1, sharex="col", figsize=(10, 5))
# convert reaction # from integers to strings
rxnstring = []
for i in range(len(posindex)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.getgasreactionstring(posindex[i] + 1))
# use horizontal bar chart
plt.subplot(211)
plt.barh(rxnstring, poscoeffs, color="orange", height=0.4)
#
# convert reaction # from integers to strings
rxnstring.clear()
fnegindex = np.flip(negindex)
fnegcoeffs = np.flip(negcoeffs)
for i in range(len(negindex)):
    # the array index starting from 0 so the actual reaction # = index + 1
    rxnstring.append(MyGasMech.getgasreactionstring(fnegindex[i] + 1))
plt.subplot(212)
plt.barh(rxnstring, fnegcoeffs, color="orange", height=0.4)
plt.xlabel("Sensitivity Coefficients")
#
plt.show()

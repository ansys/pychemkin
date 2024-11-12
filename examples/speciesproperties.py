import os

import chemkin as ck  # Chemkin
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

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
MyGasMech.tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# extract element symbols as a list
elelist = MyGasMech.elementsymbols
# extract gas species symbols as a list
specieslist = MyGasMech.speciessymbols
# list of gas species interested
plotspeclist = ["CH4", "O2", "N2"]
# find elemental compositions of selected species
print(" ")
for s in plotspeclist:
    speciesID = MyGasMech.getspecindex(s)
    print("species " + specieslist[speciesID])
    print("elemental composition")
    for elemID in range(MyGasMech.MM):
        num_elem = MyGasMech.SpeciesComposition(elemID, speciesID)
        print(f"    {elelist[elemID]:>4}: {num_elem:2d}")
    print("=" * 10)
print()
#
# plot Cv value at different temperatures for selected gas species
#
plt.figure(figsize=(12, 6))
# temperature increment
dTemp = 20.0
# number of property data points
points = 100
# curve attributes
curvelist = ["g", "b--", "r:"]
# create arrays
# species specific heat capacity at constant volume data
Cv = np.zeros(points, dtype=np.double)
# temperature data
T = np.zeros(points, dtype=np.double)
# start of the plotting loop #1
k = 0
# loop over the selected gas species
for s in plotspeclist:
    # starting temperature at 300K
    Temp = 300.0
    # loop over temperature data points
    for i in range(points):
        HeatCapacity = MyGasMech.SpeciesCv(Temp)
        ID = MyGasMech.getspecindex(s)
        T[i] = Temp
        # convert ergs to joules
        Cv[i] = HeatCapacity[ID] / ck.ergsperjoule
        Temp += dTemp
    plt.subplot(121)
    plt.plot(T, Cv, curvelist[k])
    k += 1
# plot Cv versus temperature
plt.xlabel("Temperature [K]")
plt.ylabel("Cv [J/mol-K]")
plt.legend(plotspeclist, loc="upper left")
# create arrays
# species conductivity
kappa = np.zeros(points, dtype=np.double)
# start of the plotting loop #2
k = 0
# loop over the selected gas species
for s in plotspeclist:
    # starting temperature at 300K
    Temp = 300.0
    # loop over temperature data points
    for i in range(points):
        conductivity = MyGasMech.SpeciesCond(Temp)
        ID = MyGasMech.getspecindex(s)
        T[i] = Temp
        # convert ergs to joules
        kappa[i] = conductivity[ID] / ck.ergsperjoule
        Temp += dTemp
    plt.subplot(122)
    plt.plot(T, kappa, curvelist[k])
    k += 1
# plot conductivity versus temperature
plt.xlabel("Temperature [K]")
plt.ylabel("Conductivity [J/cm-K-sec]")
plt.legend(plotspeclist, loc="upper left")
# calculate species binary diffusion coefficients
# at 2 atm and 500K
diffcoef = MyGasMech.SpeciesDiffusionCoeffs(2.0 * ck.Patm, 500.0)
ID1 = MyGasMech.getspecindex(plotspeclist[0])
ID2 = MyGasMech.getspecindex(plotspeclist[1])
c = diffcoef[ID1][ID2]
print(
    f"diffusion coefficient for {plotspeclist[0]} against {plotspeclist[1]} is {c:e} cm2/sec"
)
# display the plots
plt.show()

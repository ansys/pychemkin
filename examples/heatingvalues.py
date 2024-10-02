import os

import numpy as np  # number crunching

import chemkin as ck

# create the mechanism file with fuel species and complete combustion products only
# no reaction
# check working directory
current_dir = os.getcwd()
# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
# create a chemistry set object
MyGasMech = ck.Chemistry(label="EQ")
#
# create a new mechanism input file
#
mymechfile = os.path.join(current_dir, "fuels_chem.inp")
m = open(mymechfile, "w")
# the mechanism contains only the necessary species (fuel, oxygen, and major combustion products)
# decalre elements
m.write("ELEMENT c h o END\n")
# declare species (h2o(l) is used to compute High heating values)
m.write("SPECIES\n")
m.write("ch4 c4h10 nc5h12 nc7h16 ic8h18 nc9h20 nc10h22\n")
m.write("hmn c6h5ch3 c6h5c2h5\n")
m.write("chx mch decalin etbe mtbe\n")
m.write("o2 co2 h2o\n")
m.write("END\n")
# no reaction is needed for equilibrium calculation
m.write("REACTION\n")
m.write("END\n")
# close the mechnaism file
m.close()
#
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = mymechfile
MyGasMech.thermfile = os.path.join(
    data_dir,
    "ModelFuelLibrary",
    "Full",
    "Gasoline-Diesel-Biodiesel_PAH_NOx_therm_MFL2023.dat",
)
# pre-process
iError = MyGasMech.preprocess()
if iError == 0:
    print(ck.Color.GREEN + ">>> preprocess OK", end="\n" + ck.Color.END)
else:
    print(ck.Color.RED + ">>> preprocess failed!", end="\n" + ck.Color.END)
    exit()
#
# set pressure & temperature condition
thispressure = ck.Patm
thistemperature = 298.15
# create the unburned fuel-oxygen mixture
unburned = ck.Mixture(MyGasMech)
unburned.pressure = thispressure
unburned.temperature = thistemperature
# find the index for water vapor
watervaporID = MyGasMech.getspecindex("h2o")
# water latent heat [erg/g-water] at 298.15 [K]
latentheat = 2444.421181749129 * ck.ergsperjoule
# prepare the fuel mixture
fuel = ck.Mixture(MyGasMech)
fuel.pressure = thispressure
fuel.temperature = thistemperature
# list of fuel compositions of which the heating values will be computed
fuels = [[("ch4", 1.0)], [("c4h10", 1.0)], [("nc7h16", 0.2), ("ic8h18", 0.8)]]
# specify oxidizers = pure oxygen
oxid = ck.Mixture(MyGasMech)
oxid.X = [("o2", 1.0)]
oxid.pressure = thispressure
oxid.temperature = thistemperature
# get o2 index
oxyID = MyGasMech.getspecindex("o2")
# specify the complete combustion product species
products = ["co2", "h2o"]
# no added species
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)
#
# compute fuel heating values
#
LHV = np.zeros(len(fuels), dtype=np.double)
HHV = np.zeros_like(LHV, dtype=np.double)
fuelcount = 0
for f in fuels:
    fuel.X = f
    iError = unburned.XbyEquivalenceRatio(
        MyGasMech, fuel.X, oxid.X, add_frac, products, equivalenceratio=1.0
    )
    # get the mixture enthalpy of the initial mixture [erg/g]
    Hunburned = unburned.HML() / unburned.WTM
    # compute the complete combustion state (fixed temperature and pressure)
    # burned = ck.equilibrium(unburned, opt=1)
    burned = unburned.FindEquilibrium()
    # get the mixture enthalpy of the final mixture [erg/g]
    Hburned = burned.HML() / burned.WTM
    # get total fuel mass fraction
    fmass = 0.0e0
    bmassfrac = unburned.Y
    for i in range(MyGasMech.KK):
        if i != oxyID:
            fmass += bmassfrac[i]
    # water vapor mass fraction in the urned mixture
    wmass = burned.Y[watervaporID]
    #
    if np.isclose(fmass, 0.0, atol=1.0e-10):
        print(f">>> error finding fuel species {f} in the unburned mixture")
        exit()
    # compute the heating values [erg/g-fuel]
    LHV[fuelcount] = -(Hburned - Hunburned) / fmass
    HHV[fuelcount] = -(Hburned - (Hunburned + latentheat * wmass)) / fmass
    fuelcount += 1

# display results
print(f"Fuel Heating Values at {thistemperature} [K] and {thispressure*1.0e-6} [bar]\n")
for i in range(len(fuels)):
    print(f"fuel composition:  {fuels[i]}")
    print(f" LHV [kJ/g-fuel]:  {LHV[i] / ck.ergsperjoule / 1.0e3}")
    print(f" HHV [kJ/g-fuel]:  {HHV[i] / ck.ergsperjoule / 1.0e3}\n")

del HHV, LHV
del fuel, oxid, unburned, burned
# delete the local mechanism file just created
os.remove(mymechfile)

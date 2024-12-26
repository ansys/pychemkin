import os

import numpy as np  # number crunching

import chemkin as ck
from chemkin import Color

# check working directory
current_dir = os.getcwd()
# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
# set pressure & temperature condition
thispressure = ck.Patm
thistemperature = 298.15


#
def getwaterheatofvaporization(temp):
    """
    Compute water heat of vaporization [erg/g-water] at the given temperature
    Use the enthalpy difference between water vapor and liquid water at the temperature
    Enthalpy data depend on temperature only
    There are empirical formulas for heat of vaporization, for example, DIPPR EQ
    :param temp: temperature [K] (double scalar)
    """
    # compute water heat of vaporization
    # create a chemistry set object
    WaterMech = ck.Chemistry(label="Water Only")
    #
    # create a new mechanism input file
    #
    waterfile = os.path.join(current_dir, "water_chem.inp")
    w = open(waterfile, "w")
    # the water mechanism contains two species:
    # water vapor (H2O) and liquid water (H2O(L))
    # decalre elements
    w.write("ELEMENT H O END\n")
    w.write("SPECIES\n")
    w.write("H2O   H2O(L)\n")
    w.write("END\n")
    # no reaction is needed for thermodynamic property calculation
    w.write("REACTION\n")
    w.write("END\n")
    # close the mechnaism file
    w.close()
    # set mechanism input files
    # inclusion of the full file path is recommended
    WaterMech.chemfile = waterfile
    WaterMech.thermfile = os.path.join(
        data_dir,
        "therm.dat",
    )
    # pre-process
    iError = WaterMech.preprocess()
    if iError != 0:
        return 0.0
    # get species enthalpies [erg/mol] at 298.15 [K]
    waterenthalpies = WaterMech.SpeciesH(thistemperature)
    # compute heat of vaporization of water [erg/g-water]
    # H_water_vapor - H_liquid_water
    heatvaporization = (waterenthalpies[0] - waterenthalpies[1]) / WaterMech.WT[0]
    # remove the temporary water mechanism file
    os.remove(waterfile)
    return heatvaporization


#
# create the mechanism file with fuel species and complete combustion products only
# no reaction
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
# declare species
# ch4: Methane                  c4h10: n-Butane
# nc5h12: n-Pentane             nc7h16: n-Heptane
# ic8h18: iso-Octane            nc9h20: n-Nonane
# nc10h22: n-Decane             hmn: Heptamethylnonane
# c6h5ch3: Toluene              c6h5c2h5: Ethylbenzene
# chx: Cyclohexane              mch: Methylcyclohexane
# decalin: Decalin              ch3oh: Methanol
# c2h5oh: Alcohol               ch3och3: Dimethyl ether (DME)
# nc4h9oh: n-Butanol            mb: Mythyl butanoate
# md: Methyl decanoate          mhd: Methyl stearate
m.write("SPECIES\n")
m.write("ch4 c4h10 nc5h12 nc7h16 ic8h18 nc9h20 nc10h22\n")
m.write("hmn c6h5ch3 c6h5c2h5 chx mch decalin\n")
m.write("ch3oh c2h5oh ch3och3 nc4h9oh\n")
m.write("mb md mhd\n")
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
    print(Color.GREEN + ">>> preprocess OK", end=Color.END)
else:
    print(Color.RED + ">>> preprocess failed!", end=Color.END)
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
# water heat of vaporization [erg/g-water] at 298.15 [K]
# either call this method before creating the current Chemistry Set
# or use the activate method to switch back to the current Chemistry Set after the call
heatvaporization = getwaterheatofvaporization(thistemperature)
# switch back to MyGasMech
MyGasMech.activate()
# prepare the fuel mixture
fuel = ck.Mixture(MyGasMech)
fuel.pressure = thispressure
fuel.temperature = thistemperature
# list of fuel compositions (mole/volume fractions) of which the heating values will be computed
# [Methane, n-Butane, PRF RON 80, biodiesel]
fuels = [
    [("ch4", 1.0)],
    [("c4h10", 1.0)],
    [("nc7h16", 0.2), ("ic8h18", 0.8)],
    [("mhd", 0.9), ("ch3oh", 0.1)],
]
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
    # re-set the fuel composition (mole/volume fractions)
    fuel.X = f
    # create a soichiometric fuel-oxygen mixture
    iError = unburned.XbyEquivalenceRatio(
        MyGasMech, fuel.X, oxid.X, add_frac, products, equivalenceratio=1.0
    )
    # get the mixture enthalpy of the initial mixture [erg/g]
    Hunburned = unburned.HML() / unburned.WTM
    # compute the complete combustion state (fixed temperature and pressure)
    # this step mimics the complete burning of the initial fuel-oxygen mixture at constant pressure
    # and the subsequent cooling of the combustion prodcts back to the original temperature
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
        # no fuel species exists in the unburned mixture
        print(f">>> error no fuel species in the unburned mixture {f}")
        exit()
    # compute the heating values [erg/g-fuel]
    LHV[fuelcount] = -(Hburned - Hunburned) / fmass
    HHV[fuelcount] = -(Hburned - (Hunburned + heatvaporization * wmass)) / fmass
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

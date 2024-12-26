import os

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin homonegeous charge compression ignition (HCCI) engine model (transient)
from chemkin.engines.HCCI import HCCIengine

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set verbose mode
ck.setverbose(True)
# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the GRI mechanism
MyGasMech = ck.Chemistry(label="GRI 3.0")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
MyGasMech.thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
MyGasMech.tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
# preprocess the mechanism files
iError = MyGasMech.preprocess()
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.X = [("CH4", 0.9), ("C3H8", 0.05), ("C2H6", 0.05)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = 1.5 * ck.Patm
fuelmixture.temperature = 400.0
# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.X = [("O2", 0.21), ("N2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = 1.5 * ck.Patm
air.temperature = 400.0
# create the unburned fuel-air mixture
fresh = ck.Mixture(MyGasMech)
# products from the complete combustion of the fuel mixture and air
products = ["CO2", "H2O", "N2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 0.8
iError = fresh.XbyEquivalenceRatio(
    MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=equiv
)
if iError != 0:
    raise RuntimeError
# list the composition of the unburned fuel-air mixture
fresh.listcomposition(mode="mole")
# set mixture temperature and pressure (equivalent to setting the initial temperature and pressure of the reactor)
fresh.temperature = 447.0
fresh.pressure = 1.065 * ck.Patm
# set exhaust gas recirculation (EGR) ratio with volume fraction
EGRratio = 0.3
# compute the EGR stream composition in mole fractions
add_frac = fresh.getEGRmolefraction(EGRratio, threshold=1.0e-8)
# recreate the initial mixture with EGR
iError = fresh.XbyEquivalenceRatio(
    MyGasMech,
    fuelmixture.X,
    air.X,
    add_frac,
    products,
    equivalenceratio=equiv,
    threshold=1.0e-8,
)
# list the composition of the fuel+air+EGR mixture
fresh.listcomposition(mode="mole", bound=1.0e-8)
# HCCI engine
# create a single-zone HCCI engine object
MyEngine = HCCIengine(reactor_condition=fresh, nzones=1)
# show initial gas composition inside the reactor
MyEngine.listcomposition(mode="mole", bound=1.0e-8)
#
# set engine parameters
# cylinder bore diameter [cm]
MyEngine.bore = 12.065
# engine stroke [cm]
MyEngine.stroke = 14.005
# connecting rod length [cm]
MyEngine.connectingrod = 26.0093
# compression ratio [-]
MyEngine.compressionratio = 16.5
# engine speed [RPM]
MyEngine.RPM = 1000
# set piston pin offset distance [cm]
MyEngine.setpistonpinoffset(offset=-0.5)
# set other parameters
# simulation start CA [degree]
MyEngine.startingCA = -142.0
# simulation end CA [degree]
MyEngine.endingCA = 116.0
# list the engine parameters
MyEngine.listengineparameters()
print(f"engine displacement volume {MyEngine.getdisplacementvolume()} [cm3]")
print(f"engine clearance volume {MyEngine.getclearancevolume()} [cm3]")
print(f"number of zone(s) = {MyEngine.getnumberofzones()}")
# wall heat transfer model
# set model parameters
# "dimensionless": [<a> <b> <c> <Twall>]
# "dimensional": [<a> <b> <c> <Twall>]
# "hohenburg": [<a> <b> <c> <d> <e> <Twall>]
heattransferparameters = [0.035, 0.71, 0.0]
# set cylinder wall temperature [K]
Twall = 400.0
MyEngine.setwallheatransfer("dimensionless", heattransferparameters, Twall)
# incylinder gas velocity correlation parameter (Woschni)
# [<C11> <C12> <C2> <swirl ratio>]
GVparameters = [2.28, 0.308, 3.24, 0.0]
MyEngine.setgasvelocitycorrelation(GVparameters)
# set piston head top surface area [cm2]
MyEngine.setpistonheadarea(area=124.75)
# set cylinder clearance surface area [cm2]
MyEngine.setcylinderheadarea(area=123.5)
# output controls
# set the number of crank angles between saving solution
MyEngine.CAstepforsavingsolution = 0.5
# set the number of crank angles between printing solution
MyEngine.CAstepforprintingsolution = 10.0
# turn ON adaptive solution saving
MyEngine.adaptivesolutionsaving(mode=True, steps=20)
# turn OFF adaptive solution saving
# MyEngine.adaptivesolutionsaving(mode=False)
# set tolerance
MyEngine.settolerances(absolute_tolerance=1.0e-12, relative_tolerance=1.0e-10)
# get solver parameters
ATOL, RTOL = MyEngine.tolerances
print(f"default absolute tolerance = {ATOL}")
print(f"default relative tolerance = {RTOL}")
# turn on the force non-negative solutions option in the solver
MyEngine.forcenonnegative = True
# specify the ignition definitions
# ck.showignitiondefinition()
MyEngine.setignitiondelay(method="T_inflection")
# stop the simulation when ignition is detected
# MyEngine.stopafterignition()
# show solver option
# show the number of crank angles between printng solution
print(f"crank angles between solution printing: {MyEngine.CAstepforprintingsolution}")
# show other transient solver setup
print(f"forced non-negative solution values: {MyEngine.forcenonnegative}")
# show the additional keywords given by user
MyEngine.showkeywordinputlines()
# run the single-zone HCCI engine model
runstatus = MyEngine.run()
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> RUN COMPLETED <<<", end=Color.END)
# get ignition delay "time"
delayCA = MyEngine.getignitiondelay()
print(f"ignition delay CA = {delayCA} [degree]")
# get heat release information
HR10, HR50, HR90 = MyEngine.getengineheatrelease()
print("Engine Heat Release Information")
print(f"10% heat release CA = {HR10} [degree]")
print(f"50% heat release CA = {HR50} [degree]")
print(f"90% heat release CA = {HR90} [degree]\n")
# post-process the solutions
MyEngine.processenginesolution()
# get the number of solution time points
solutionpoints = MyEngine.getnumbersolutionpoints()
print(f"number of solution points = {solutionpoints}")
# get the time profile
timeprofile = MyEngine.getsolutionvariableprofile("time")
# convert time to crank angle
CAprofile = np.zeros_like(timeprofile, dtype=np.double)
count = 0
for t in timeprofile:
    CAprofile[count] = MyEngine.toCA(timeprofile[count])
    count += 1
# get the cylinder pressure profile
presprofile = MyEngine.getsolutionvariableprofile("pressure")
presprofile *= 1.0e-6
# get the volume profile
volprofile = MyEngine.getsolutionvariableprofile("volume")
# create arrays for mixture density, NO mole fraction, and mixture specific heat capacity
denprofile = np.zeros_like(timeprofile, dtype=np.double)
Cpprofile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the mixture at the time point
    solutionmixture = MyEngine.getsolutionmixtureatindex(solution_index=i)
    # get gas density [g/cm3]
    denprofile[i] = solutionmixture.RHO
    # get mixture specific heat capacity profile [erg/mole-K]
    Cpprofile[i] = solutionmixture.CPBL() / ck.ergsperjoule * 1.0e-3
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.subplot(221)
plt.plot(CAprofile, presprofile, "r-")
plt.ylabel("Pressure [bar]")
plt.subplot(222)
plt.plot(CAprofile, volprofile, "b-")
plt.ylabel("Volume [cm3]")
plt.subplot(223)
plt.plot(CAprofile, denprofile, "g-")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Density [g/cm3]")
plt.subplot(224)
plt.plot(CAprofile, Cpprofile, "m-")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Cp [kJ/mole]")
# display the plots
plt.show()

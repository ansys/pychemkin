import os
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin spark ignition (SI) engine model (transient)
from chemkin.engines.SI import SIengine

# check working directory
current_dir = os.getcwd()
print("current working directory: " + current_dir)
# set verbose mode
ck.setverbose(True)
# set mechanism directory (the default chemkin mechanism data directory)
data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
mechanism_dir = data_dir
# create a chemistry set based on the gasoline 14 components mechanism
MyGasMech = ck.Chemistry(label="Gasoline")
# set mechanism input files
# inclusion of the full file path is recommended
MyGasMech.chemfile = os.path.join(mechanism_dir, "gasoline_14comp_WBencrypt.inp")
# preprocess the mechanism files
iError = MyGasMech.preprocess()
print("mechanism information:")
print(f"number of gas species = {MyGasMech.KK:d}")
print(f"number of gas reactions = {MyGasMech.IIGas:d}")
# create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
# create the fuel mixture
fuelmixture = ck.Mixture(MyGasMech)
# set fuel composition
fuelmixture.X = [("ic8h18", 0.9), ("nc7h16", 0.1)]
# setting pressure and temperature is not required in this case
fuelmixture.pressure = ck.Patm
fuelmixture.temperature = 353.0
# create the oxidizer mixture: air
air = ck.Mixture(MyGasMech)
air.X = [("o2", 0.21), ("n2", 0.79)]
# setting pressure and temperature is not required in this case
air.pressure = ck.Patm
air.temperature = 353.0
# create the unburned fuel-air mixture
fresh = ck.Mixture(MyGasMech)
# products from the complete combustion of the fuel mixture and air
products = ["co2", "h2o", "n2"]
# species mole fractions of added/inert mixture. can also create an additives mixture here
add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros
# mean equivalence ratio
equiv = 1.0
iError = fresh.XbyEquivalenceRatio(
    MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=equiv
)
if iError != 0:
    raise RuntimeError
# list the composition of the unburned fuel-air mixture
fresh.listcomposition(mode="mole")
# set mixture temperature and pressure (equivalent to setting the initial temperature and pressure of the reactor)
fresh.temperature = fuelmixture.temperature
fresh.pressure = fuelmixture.pressure
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
# SI engine
# create an SI engine object
MyEngine = SIengine(reactor_condition=fresh)
# show initial gas composition inside the reactor
MyEngine.listcomposition(mode="mole", bound=1.0e-8)
#
# set engine parameters
# cylinder bore diameter [cm]
MyEngine.bore = 8.5
# engine stroke [cm]
MyEngine.stroke = 10.82
# connecting rod length [cm]
MyEngine.connectingrod = 17.853
# compression ratio [-]
MyEngine.compressionratio = 12
# engine speed [RPM]
MyEngine.RPM = 600
# set other parameters
# simulation start CA [degree]
MyEngine.startingCA = -120.2
# simulation end CA [degree]
MyEngine.endingCA = 139.8
# list the engine parameters
MyEngine.listengineparameters()
print(f"engine displacement volume {MyEngine.getdisplacementvolume()} [cm3]")
print(f"engine clearance volume {MyEngine.getclearancevolume()} [cm3]")
# set mass burned fraction profile for the SI engine combustion
# >>> option 1
# use Wiebe function
# start of combustion crank angle
MyEngine.setburntiming(SOC=-14.5, duration=45.6)
MyEngine.wiebeparameters(n=4.0, b=7.0)
# >>> option 2
# use anchor points
# MyEngine.setanchorpoints(CA10=5.2, CA50=14.22, CA90=22.01)
# >>> option 3
# use burned mass fraction profile data
# start of combustion crank angle
# MyEngine.setburntiming(SOC=-14.5, duration=45.6)
# nMBFdata = 9
# MBangles = np.zeros(nMBFdata, dtype=np.double)
# MBFrac = np.zeros_like(MBangles, dtype=np.double)
# normalized crank angles (CA - SOC) / duration
# MBangles = [0.0, 0.26985, 0.3371, 0.3742, 0.4322, 0.63, 0.704, 0.8, 1.0]
# mass burned fraction
# MBFrac =   [0.0, 0.01, 0.03, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0]
# MyEngine.setmassburnedprofile(crankangles=MBangles, fractions=MBFrac)
# <<<
# set minimum zonal mass [g]
MyEngine.setminimumzonemass(minmass=1.0e-5)
# wall heat transfer model
# set model parameters
# "dimensionless": [<a> <b> <c> <Twall>]
# "dimensional": [<a> <b> <c> <Twall>]
# "hohenburg": [<a> <b> <c> <d> <e> <Twall>]
heattransferparameters = [0.1, 0.8, 0.0]
# set cylinder wall temperature [K]
Twall = 434.0
MyEngine.setwallheatransfer("dimensionless", heattransferparameters, Twall)
# incylinder gas velocity correlation parameter (Woschni)
# [<C11> <C12> <C2> <swirl ratio>]
GVparameters = [2.28, 0.318, 0.324, 0.0]
MyEngine.setgasvelocitycorrelation(GVparameters)
# set piston head top surface area [cm2]
MyEngine.setpistonheadarea(area=56.75)
# set cylinder clearance surface area [cm2]
MyEngine.setcylinderheadarea(area=56.75)
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
MyEngine.settolerances(absolute_tolerance=1.0e-15, relative_tolerance=1.0e-6)
# get solver parameters
ATOL, RTOL = MyEngine.tolerances
print(f"default absolute tolerance = {ATOL}")
print(f"default relative tolerance = {RTOL}")
# turn on the force non-negative solutions option in the solver
MyEngine.forcenonnegative = True
# show solver option
# show the number of crank angles between printng solution
print(f"crank angles between solution printing: {MyEngine.CAstepforprintingsolution}")
# show other transient solver setup
print(f"forced non-negative solution values: {MyEngine.forcenonnegative}")
# show the additional keywords given by user
MyEngine.showkeywordinputlines()
# set the start wall time
start_time = time.time()
# run the single-zone HCCI engine model
runstatus = MyEngine.run()
# compute the total runtime
runtime = time.time() - start_time
# check run status
if runstatus != 0:
    # run failed!
    print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
    exit()
# run success!
print(Color.GREEN + ">>> RUN COMPLETED <<<", end=Color.END)
print(f"total simulation duration: {runtime} [sec]")
#
# post-process the solution profiles in selected zone
unburnedzone = 1
burnedzone = 2
zonestrings = ["Unburned Zone", "Burned Zone"]
thiszone = burnedzone
MyEngine.processenginesolution(zoneID=thiszone)
plottitle = zonestrings[thiszone - 1] + " Solution"
# post-process cylinder-averged solution profiles
# MyMZEngine.processaverageenginesolution()
# plottitle = "Cylinder Averaged Solution"
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
# get zonal temperature profile [K]
tempprofile = MyEngine.getsolutionvariableprofile("temperature")
# get the zonal volume profile
volprofile = MyEngine.getsolutionvariableprofile("volume")
# create arrays for zonal CO mole fraction
# find CO species index
COindex = MyGasMech.getspecindex("co")
COprofile = np.zeros_like(timeprofile, dtype=np.double)
# loop over all solution time points
for i in range(solutionpoints):
    # get the zonal mixture at the time point
    solutionmixture = MyEngine.getsolutionmixtureatindex(solution_index=i)
    # get zonal CO mole fraction
    COprofile[i] = solutionmixture.X[COindex]
#
# post-process cylinder-averged solution
MyEngine.processaverageenginesolution()
# get the cylinder volume profile
cylindervolprofile = MyEngine.getsolutionvariableprofile("volume")
# create arrays for cylinder-averaged mixture temperature
cylindertempprofile = MyEngine.getsolutionvariableprofile("temperature")
#
# plot the profiles
plt.subplots(2, 2, sharex="col", figsize=(12, 6))
plt.suptitle(plottitle, fontsize=16)
plt.subplot(221)
plt.plot(CAprofile, presprofile, "r-")
plt.ylabel("Pressure [bar]")
plt.subplot(222)
plt.plot(CAprofile, volprofile, "b-")
plt.plot(CAprofile, cylindervolprofile, "b--")
plt.ylabel("Cylinder Volume [cm3]")
plt.subplot(223)
plt.plot(CAprofile, tempprofile, "g-")
plt.plot(CAprofile, cylindertempprofile, "g--")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("Mixture Temperature [K]")
plt.subplot(224)
plt.plot(CAprofile, COprofile, "m-")
plt.xlabel("Crank Angle [degree]")
plt.ylabel("CO Mole Fraction")
# display the plots
plt.show()

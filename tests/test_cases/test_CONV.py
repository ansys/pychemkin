import os

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin batch reactor models (transient)
from chemkin.batchreactors.batchreactor import (
    GivenVolumeBatchReactor_EnergyConservation,
)

# import pytest


# @pytest.mark.skip(reason="Causes segfault")
def test_CONV():
    # check working directory
    current_dir = os.getcwd()
    print("current working directory: " + current_dir)
    # set verbose mode
    ck.setverbose(True)
    # set mechanism directory (the default chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # create a chemistry set based on the diesel 14 components mechanism
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
    fuelmixture.X = [("CH4", 1.0)]
    # setting pressure and temperature is not required in this case
    fuelmixture.pressure = 5.0 * ck.Patm
    fuelmixture.temperature = 1500.0
    # create the oxidizer mixture: air
    air = ck.Mixture(MyGasMech)
    air.X = [("O2", 0.21), ("N2", 0.79)]
    # setting pressure and temperature is not required in this case
    air.pressure = 5.0 * ck.Patm
    air.temperature = 1500.0
    # create the premixed mixture to be defined
    premixed = ck.Mixture(MyGasMech)
    # products from the complete combustion of the fuel mixture and air
    products = ["CO2", "H2O", "N2"]
    # species mole fractions of added/inert mixture. can also create an additives mixture here
    add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros
    iError = premixed.XbyEquivalenceRatio(
        MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=0.7
    )
    if iError != 0:
        raise RuntimeError
    # list the composition of the premixed mixture
    premixed.listcomposition(mode="mole")
    # set mixture temperature and pressure (equivalent to setting the initial temperature and pressure of the reactor)
    premixed.temperature = 800.0
    premixed.pressure = 3.0 * ck.Patm
    # Rapid Compression Machine
    # create a constant volume batch reactor (with energy equation)
    #
    MyCONV = GivenVolumeBatchReactor_EnergyConservation(premixed, label="RCM")
    # set the initial reactor temperature (see the warning message in the run output)
    # MyCONV.temperature = 800.0  # K
    # show initial gas composition inside the reactor
    MyCONV.listcomposition(mode="mole")
    # set other reactor properties
    MyCONV.volume = 10.0  # cm3
    MyCONV.time = 0.1  # sec
    # output controls
    # set timestep between saving solution
    MyCONV.timestepforsavingsolution = 0.001
    # turn ON saving to XML solution file (default)
    MyCONV.XML_Output = True
    # turn ON adaptive solution saving
    MyCONV.adaptivesolutionsaving(mode=True, value_change=100, target="TEMPERATURE")
    # turn OFF adaptive solution saving
    # MyCONV.adaptivesolutionsaving(mode=False)
    # set tolerance
    MyCONV.settolerances(absolute_tolerance=1.0e-10, relative_tolerance=1.0e-8)
    # get solver parameters
    ATOL, RTOL = MyCONV.tolerances
    print(f"default absolute tolerance = {ATOL}")
    print(f"default relative tolerance = {RTOL}")
    # turn on the force non-negative solutions option in the solver
    MyCONV.forcenonnegative = True
    # set RCM volume profile
    # number of profile data points
    npoints = 3
    # position array of the profile data
    x = np.zeros(npoints, dtype=np.double)
    # value array of the profile data
    volprofile = np.zeros_like(x, dtype=np.double)
    # set reactor volume data points
    x = [0.0, 0.01, 2.0]  # [sec]
    volprofile = [10.0, 4.0, 4.0]  # [cm3]
    # set the volume profile
    MyCONV.setvolumeprofile(x, volprofile)
    # change timestep between saving solution
    MyCONV.timestepforsavingsolution = 0.01
    # turn OFF adaptive solution saving
    # MyCONV.adaptivesolutionsaving(mode=False)
    # set ignition delay
    # get ignition definitions
    # ck.showignitiondefinition()
    MyCONV.setignitiondelay(method="T_inflection")
    # stop the simulation when ignition is detected
    # MyCONV.stopafterignition()
    # show solver option
    print(f"timestep between solution printing: {MyCONV.timestepforprintingsolution}")
    # show timestep between printing solution
    print(f"forced non-negative solution values: {MyCONV.forcenonnegative}")
    # show the additional keywords given by user
    MyCONV.showkeywordinputlines()
    # run the CONV reactor model
    runstatus = MyCONV.run()
    # check run status
    if runstatus != 0:
        # run failed!
        print(Color.RED + ">>> RUN FAILED <<<", end=Color.END)
        exit()
    # run success!
    print(Color.GREEN + ">>> RUN COMPLETED <<<", end=Color.END)
    # get ignition delay time
    delaytime = MyCONV.getignitiondelay()
    print(f"ignition delay time = {delaytime} [msec]")
    # post-process the solutions
    MyCONV.processsolution()
    # get the number of solution time points
    solutionpoints = MyCONV.getnumbersolutionpoints()
    print(f"number of solution points = {solutionpoints}")
    # mix1 = MyCONV.getsolutionmixture(0.03)
    # mix2 = MyCONV.getsolutionmixture(0.05)
    # get the time profile
    timeprofile = MyCONV.getsolutionvariableprofile("time")
    # get the temperature profile
    tempprofile = MyCONV.getsolutionvariableprofile("temperature")
    # get CH4 mass fraction profile
    # CH4profile = MyCONV.getsolutionvariableprofile('CH4')
    CH4profile = np.zeros_like(timeprofile, dtype=np.double)
    #
    # more involving post-processing by using Mixtures
    #
    # create arrays for CH4 ROP and mixture viscosity
    CH4ROPprofile = np.zeros_like(timeprofile, dtype=np.double)
    viscprofile = np.zeros_like(timeprofile, dtype=np.double)
    # CurrentROP = np.zeros(MyGasMech.KK, dtype=np.double)
    # find co species index
    CH4_index = MyGasMech.getspecindex("CH4")
    # loop over all solution time points
    for i in range(solutionpoints):
        # get the mixture at the time point
        solutionmixture = MyCONV.getsolutionmixtureatindex(solution_index=i)
        # get CH4 mole fraction profile
        CH4profile[i] = solutionmixture.X[CH4_index]
        # get CH4 ROP profile
        currentROP = solutionmixture.ROP()
        CH4ROPprofile[i] = currentROP[CH4_index]
        # get mixture vicosity profile
        viscprofile[i] = solutionmixture.mixtureviscosity()
    # plot the profiles
    plt.subplots(2, 2, sharex="col", figsize=(12, 6))
    plt.subplot(221)
    plt.plot(timeprofile, tempprofile, "r-")
    plt.ylabel("Temperature [K]")
    plt.subplot(222)
    plt.plot(timeprofile, CH4profile, "b-")
    plt.ylabel("CH4 Mole Fraction")
    plt.subplot(223)
    plt.plot(timeprofile, CH4ROPprofile, "g-")
    plt.xlabel("time [sec]")
    plt.ylabel("CH4 Production Rate [mol/cm3-sec]")
    plt.subplot(224)
    plt.plot(timeprofile, viscprofile, "m-")
    plt.xlabel("time [sec]")
    plt.ylabel("Mixture Viscosity [g/cm-sec]")
    # display the plots
    plt.savefig("CONV_solution.png", bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    test_CONV()

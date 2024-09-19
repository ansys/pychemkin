import os
import time

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

import chemkin as ck  # Chemkin
from chemkin import Color

# chemkin batch reactor models (transient)
from chemkin.batchreactor import GivenPressureBatchReactor_EnergyConservation

# import pytest


# @pytest.mark.skip(reason="Temporarily disabled for demonstration purposes")
def test_ignitiondelay():
    # check working directory
    current_dir = os.getcwd()
    print("current working directory: " + current_dir)
    # set verbose mode
    ck.setverbose(True)
    # set mechanism directory (the default chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # create a chemistry set based on GRI 3.0
    diesel = ck.Chemistry(label="diesel 14comp")
    # set mechanism input files
    # inclusion of the full file path is recommended
    diesel.chemfile = os.path.join(mechanism_dir, "gasoline_14comp_WBencrypt.inp")
    # preprocess the mechanism files
    iError = diesel.preprocess()
    # create a premixed fuel-oxidizer mixture by assigning the equivalence ratio
    # create the fuel mixture
    fuelmixture = ck.Mixture(diesel)
    # set fuel = composition PRF 60
    fuelmixture.X = [("ic8h18", 0.6), ("nc7h16", 0.4)]
    # setting pressure and temperature is not required in this case
    fuelmixture.pressure = 5.0 * ck.Patm
    fuelmixture.temperature = 1500.0
    # create the oxidizer mixture: air
    air = ck.Mixture(diesel)
    air.X = [("o2", 0.21), ("n2", 0.79)]
    # setting pressure and temperature is not required in this case
    air.pressure = 5.0 * ck.Patm
    air.temperature = 1500.0
    # create the premixed mixture to be defined
    premixed = ck.Mixture(diesel)
    # products from the complete combustion of the fuel mixture and air
    products = ["co2", "h2o", "n2"]
    # species mole fractions of added/inert mixture. can also create an additives mixture here
    add_frac = np.zeros(diesel.KK, dtype=np.double)  # no additives: all zeros
    iError = premixed.XbyEquivalenceRatio(
        diesel, fuelmixture.X, air.X, add_frac, products, equivalenceratio=1.0
    )
    if iError != 0:
        raise RuntimeError
    # list the composition of the premixed mixture
    premixed.listcomposition(mode="mole")
    # set mixture temperature and pressure (equivalent to setting the initial temperature and pressure of the reactor)
    premixed.pressure = 40.0 * ck.Patm
    premixed.temperature = 700.0
    #
    # create a constant pressure batch reactor (with energy equation)
    #
    MyCONP = GivenPressureBatchReactor_EnergyConservation(premixed, label="CONP")
    # set the initial reactor temperature (see the warning message in the run output)
    MyCONP.temperature = 950.0  # K
    # show initial gas composition inside the reactor
    MyCONP.listcomposition(mode="mole")
    # set other reactor parameters
    MyCONP.volume = 10.0  # cm3
    MyCONP.time = 1.0  # sec
    # output controls
    # set timestep between saving solution
    MyCONP.timestepforsavingsolution = 0.001
    # turn ON saving to XML solution file (default)
    MyCONP.XML_Output = True
    # turn ON adaptive solution saving
    MyCONP.adaptivesolutionsaving(mode=True, value_change=100, target="TEMPERATURE")
    # set tolerance
    MyCONP.settolerances(absolute_tolerance=1.0e-10, relative_tolerance=1.0e-8)
    # change timestep between saving solution
    MyCONP.timestepforsavingsolution = 0.01
    # set ignition delay
    # ck.showignitiondefinition()
    MyCONP.setignitiondelay(method="T_inflection")
    # stop after ignition is detected (not recommended for ignition delay time calculations)
    # MyCONP.stopafterignition()
    # show solver option
    print(f"timestep between solution printing: {MyCONP.timestepforprintingsolution}")
    # show timestep between printing solution
    print(f"forced non-negative solution values: {MyCONP.forcenonnegative}")
    #
    # loop over initial reactor temperature to create an ignition delay time plot
    #
    npoints = 20
    delta_temp = 20.0
    init_temp = premixed.temperature
    delaytime = np.zeros(npoints, dtype=np.double)
    temp_inv = np.zeros_like(delaytime, dtype=np.double)
    # set the start wall time
    start_time = time.time()
    for i in range(npoints):
        # update the initial reactor temperature
        MyCONP.temperature = init_temp  # K
        # show the additional keywords given by user
        # MyCONP.showkeywordinputlines()
        # run the reactor model
        runstatus = MyCONP.run()
        #
        if runstatus == 0:
            # plot 1/T instead of T
            temp_inv[i] = 1.0e0 / init_temp
            # get ignition delay time
            delaytime[i] = MyCONP.getignitiondelay()
            print(f"ignition delay time = {delaytime[i]} [msec]")
        else:
            # if get this, most likely the END time is too short
            print(Color.RED + ">>> RUN FAILED <<<", end="\n" + Color.END)
        init_temp += delta_temp
    # compute the total runtime
    runtime = time.time() - start_time
    print(f"total simulation duration: {runtime} [sec] for {npoints} cases")
    # create an ignition delay versus 1/T plot for the PRF fuel (should exhibit the NTC region)
    plt.rcParams.update({"figure.autolayout": True})
    plt.semilogy(temp_inv, delaytime, "bs--")
    plt.xlabel("1/T [1/K]")
    plt.ylabel("Ignition delay time [msec]")
    plt.savefig("ignitiondelays.png", bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    test_ignitiondelay()

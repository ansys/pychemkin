import os

import chemkin as ck  # Chemkin
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching
import pytest


@pytest.mark.skip(reason="Causes segfault")
def test_equilibriumcomposition():
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
    assert iError == 0
    # create the fuel mixture
    fuel = ck.Mixture(MyGasMech)
    # set mole fraction
    fuel.X = [("CH4", 0.8), ("H2", 0.2)]
    fuel.temperature = 300.0
    fuel.pressure = ck.Patm  # 1 atm
    # create the air mixture
    air = ck.Mixture(MyGasMech)
    # set mass fraction
    air.Y = [("O2", 0.23), ("N2", 0.77)]
    air.temperature = 300.0
    air.pressure = ck.Patm  # 1 atm
    # mix the fuel and the air with an air-fuel ratio of 17.19 (almost stoichiometric?)
    mixture_recipe = [(fuel, 1.0), (air, 17.19)]
    # create the new mixture (the air-fuel ratio is by mass)
    premixed = ck.isothermalmixing(
        recipe=mixture_recipe, mode="mass", finaltemperature=300.0
    )
    # find the equilibrium composition at different temperature
    # and create a NO mole fraction versus temperature plot
    # NO species index
    NO_index = MyGasMech.getspecindex("NO")
    # set up plotting temperatures
    Temp = 500.0
    dTemp = 20.0
    points = 100
    # curve
    T = np.zeros(points, dtype=np.double)
    NO = np.zeros_like(T, dtype=np.double)
    # start the temperature loop
    for k in range(points):
        # reset mixture temperature
        premixed.temperature = Temp
        # find the equilibrium state mixture at the given mixture temperature and pressure
        eqstate = ck.equilibrium(premixed, opt=1)
        #
        NO[k] = eqstate.X[NO_index] * 1.0e6  # convert to ppm
        T[k] = Temp
        Temp += dTemp
    # create plot
    plt.plot(T, NO, "bs--", markersize=3, markevery=4)
    plt.xlabel("Temperature [K]")
    plt.ylabel("NO [ppm]")
    plt.savefig("NO_equilibrium.png", bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    test_equilibriumcomposition()

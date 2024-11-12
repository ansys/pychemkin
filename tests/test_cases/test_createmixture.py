import copy
import os

import chemkin as ck  # Chemkin
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# import pytest


# @pytest.mark.skip(reason="Temporarily disabled for demonstration purposes")
def test_createmixture():
    # check working directory
    current_dir = os.getcwd()
    print("current working directory: " + current_dir)
    # set verbose mode
    ck.setverbose(True)
    # set mechanism directory (the default chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # create a chemistry set based on the C2 NOx mechanism
    MyGasMech = ck.Chemistry(label="C2 NOx")
    # set mechanism input files
    # inclusion of the full file path is recommended
    MyGasMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")
    # this mechanism file contains all the necessary thermodynamic and transport data
    # therefore no need to specify the therm and the tran data files
    # instruct the preprocessor to include the transport properties (when the tran data file is not provided)
    MyGasMech.preprocesstransportdata()
    # preprocess the mechanism files
    iError = MyGasMech.preprocess()
    assert iError == 0
    # get species molecular masses as numpy 1D double array
    WT = MyGasMech.WT
    print(WT)
    # create a mixture called premixed based on the MyGasMech chemistry set
    premixed = ck.Mixture(MyGasMech)
    # set mixture pressure [dynes/cm2]
    mixpressure = 2.0  # given in atm
    # convert to dynes/cm2
    premixed.pressure = mixpressure * ck.Patm
    # set mixture temperature [K]
    premixed.temperature = 500.0
    # molar composition of the mixture
    mixture_recipe = [("CH4", 0.08), ("N2", 0.6), ("O2", 0.2), ("H2O", 0.12)]
    # set mixture mole fractions
    premixed.X = mixture_recipe
    # find the mixture mean molecular mass
    print(f"mean molecular mass = {premixed.WTM:f} gm/mole")
    print("=" * 40)
    # convert automatically to mass fractions
    print("mixture mass fractions (raw data):")
    print(str(premixed.Y))
    # switch back to mole fractions
    print("\nmixture mole fractions (raw data):")
    print(str(premixed.X))
    print("\nformatted mixture composition output:")
    print("=" * 40)
    premixed.listcomposition(mode="mole")
    print("=" * 40)
    # create a 'hard' copy of the premixed mixture ('soft' copying, i.e., anotherpremixed = premixed, works, too)
    anotherpremixed = copy.deepcopy(premixed)
    print("\nformatted mixture composition of the copied mixture:")
    anotherpremixed.listcomposition(mode="mole")
    print("=" * 40)
    #
    # test mixture properties
    #
    plt.subplots(2, 2, sharex="col", figsize=(9, 9))
    dTemp = 20.0
    points = 100
    # curve attributes
    curvelist = ["g", "b--", "r:"]
    # list of pressure values for the plot
    press = [1.0, 5.0, 10.0]  # given in atm
    # create arrays for the plot
    # mixture density
    rho = np.zeros(points, dtype=np.double)
    # mixture enthalpy
    enthalpy = np.zeros_like(rho, dtype=np.double)
    # mixture viscosity
    visc = np.zeros_like(rho, dtype=np.double)
    # mixture averaged diffusion coefficient of CH4
    diff_CH4 = np.zeros_like(rho, dtype=np.double)
    # species index of CH4
    CH4_index = MyGasMech.getspecindex("CH4")
    # temperature data
    T = np.zeros_like(rho, dtype=np.double)
    # start of the plotting loop #1
    k = 0
    # loop over the pressure values
    for j in range(len(press)):
        # starting temperature at 300K
        temp = 300.0
        # loop over temperature data points
        for i in range(points):
            # set mixture pressure [dynes/cm2]
            premixed.pressure = press[j] * ck.Patm
            # set mixture temperature [K]
            premixed.temperature = temp
            # get mixture density [gm/cm3]
            rho[i] = premixed.RHO
            # get mixture enthalpy [ergs/mol] and convert it to [kJ/mol]
            enthalpy[i] = premixed.HML() * 1.0e-3 / ck.ergsperjoule
            # get mixture viscosity [gm/cm-sec]
            visc[i] = premixed.mixtureviscosity()
            # get mixture-averaged diffusion coefficient of CH4 [cm2/sec]
            diffcoeffs = premixed.mixturediffusioncoeffs()
            diff_CH4[i] = diffcoeffs[CH4_index]
            T[i] = temp
            temp += dTemp
        # create sub plots
        # plot mixture density versus temperature
        plt.subplot(221)
        plt.plot(T, rho, curvelist[k])
        plt.ylabel("Density [g/cm3]")
        plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper right")
        # plot mixture enthalpy versus temperature
        plt.subplot(222)
        plt.plot(T, enthalpy, curvelist[k])
        plt.ylabel("Enthalpy [kJ/mole]")
        plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper left")
        # plot mixture viscosity versus temperature
        plt.subplot(223)
        plt.plot(T, visc, curvelist[k])
        plt.xlabel("Temperature [K]")
        plt.ylabel("Viscosity [g/cm-sec]")
        plt.legend(("1 atm", "5 atm", "10 atm"), loc="lower right")
        # plot mixture averaged CH4 diffusion coefficient versus temperature
        plt.subplot(224)
        plt.plot(T, diff_CH4, curvelist[k])
        plt.xlabel("Temperature [K]")
        plt.ylabel(r"$D_{CH_4}$ [cm2/sec]")
        plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper left")
        k += 1
    # plot legends
    plt.legend(("1 atm", "5 atm", "10 atm"), loc="upper left")
    # display the plots
    plt.savefig("mixtureproperties.png", bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    test_createmixture()

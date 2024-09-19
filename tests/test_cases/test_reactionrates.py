import os

import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching
#import pytest

import chemkin as ck  # Chemkin


# @pytest.mark.skip(reason="Temporarily disabled for demonstration purposes")
def test_reactionrates():
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
        MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=1.0
    )
    if iError != 0:
        raise RuntimeError
    # list the composition of the premixed mixture
    premixed.listcomposition(mode="mole")
    #
    # compute reaction rates
    #
    # temperature and pressure are requires to compute the reaction rates
    premixed.pressure = 5.0 * ck.Patm
    premixed.temperature = 1600.0
    # get the net species molar rates of production [mole/cm3-sec]
    rop = premixed.ROP()
    # list the nonzero rates in descending order
    print()
    specrate_order, species_rates = premixed.listROP()
    # get the forward and the reverse rates of each reaction
    kf, kr = premixed.RxnRates()
    print()
    print(f"reverse reaction rates: (raw values of all {MyGasMech.IIGas:d} reactions)")
    print(str(kr))
    print("=" * 40)
    # list the nonzero net reaction rates
    rxn_order, net_rxn_rates = premixed.listreactionrates()
    # create a rate plot
    plt.rcParams.update({"figure.autolayout": True})
    plt.subplots(2, 1, sharex="col", figsize=(10, 5))
    # convert reaction # from integers to strings
    rxnstring = []
    for i in range(len(rxn_order)):
        # the array index starting from 0 so the actual reaction # = index + 1
        rxnstring.append(MyGasMech.getgasreactionstring(rxn_order[i] + 1))
    # use horizontal bar chart
    plt.subplot(211)
    plt.barh(rxnstring, net_rxn_rates, color="blue", height=0.4)
    # use log scale on x axis
    plt.xscale("symlog")
    # plt.ylabel('reaction')
    plt.text(-3.0e-4, 0.5, "T = 1600K", fontsize=10)
    # change the mixture temperature
    premixed.temperature = 1800.0
    # get the list the nonzero net reaction rates at the new temperature
    rxn_order, net_rxn_rates = premixed.listreactionrates()
    plt.subplot(212)
    # convert reaction # from integers to strings
    rxnstring.clear()
    for i in range(len(rxn_order)):
        # the array index starting from 0 so the actual reaction # = index + 1
        rxnstring.append(MyGasMech.getgasreactionstring(rxn_order[i] + 1))
    plt.barh(rxnstring, net_rxn_rates, color="orange", height=0.4)
    plt.xlabel("reaction rate [mole/cm3-sec]")
    # plt.ylabel('reaction')
    plt.text(-3.0e-4, 0.5, "T = 1800K", fontsize=10)
    # use log scale on x axis
    plt.xscale("symlog")
    plt.savefig("reactionrates.png", bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    test_reactionrates()

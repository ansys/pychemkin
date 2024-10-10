import os

import pytest

import chemkin as ck  # Chemkin
from chemkin import Color


@pytest.mark.skip(reason="Temporarily disabled for demonstration purposes")
def test_loadmechanism():
    # check working directory
    current_dir = os.getcwd()
    print("current working directory: " + current_dir)
    # set verbose mode
    ck.setverbose(True)
    # set mechanism directory (the default chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # set mechanism input files
    # inclusion of the full file path is recommended
    chemfile = os.path.join(mechanism_dir, "grimech30_chem.inp")
    thermfile = os.path.join(mechanism_dir, "grimech30_thermo.dat")
    tranfile = os.path.join(mechanism_dir, "grimech30_transport.dat")
    # create a chemistry set based on GRI 3.0
    MyGasMech = ck.Chemistry(
        chem=chemfile, therm=thermfile, tran=tranfile, label="GRI 3.0"
    )
    # preprocess the mechanism files
    iError = MyGasMech.preprocess()
    print()
    if iError != 0:
        print(f"PreProcess: error encountered...code = {iError:d}")
        print(f"see the summary file {MyGasMech.summaryfile} for details")
        exit()
    else:
        print(Color.GREEN + "PreProcess success!!", end=Color.END)
        print("mechanism information:")
        print(f"number of elements = {MyGasMech.MM:d}")
        print(f"number of gas species = {MyGasMech.KK:d}")
        print(f"number of gas reactions = {MyGasMech.IIGas:d}")

    print(f"\nelement and species information of mechanism {MyGasMech.label}")
    print("=" * 50)
    # extract element symbols as a list
    elelist = MyGasMech.elementsymbols
    # get atomic masses as numpy 1D double array
    AWT = MyGasMech.AWT
    # print element information
    for k in range(len(elelist)):
        print(f"element # {k+1:3d}: {elelist[k]:16} mass = {AWT[k]:f}")

    print("=" * 50)
    # extract gas species symbols as a list
    specieslist = MyGasMech.speciessymbols
    # get species molecular masses as numpy 1D double array
    WT = MyGasMech.WT
    # print gas species information
    for k in range(len(specieslist)):
        print(f"species # {k+1:3d}: {specieslist[k]:16} mass = {WT[k]:f}")
    print("=" * 50)
    #
    # load the second mechanism
    #
    # set the 2nd mechanism directory (the default chemkin mechanism data directory)
    mechanism_dir = data_dir
    # create a chemistry set based on C2_NOx using an alternative method
    My2ndMech = ck.Chemistry(label="C2 NOx")
    # set mechanism input files individually
    # this mechanism file contains all the necessary thermodynamic and transport data
    # therefore no need to specify the therm and the tran data files
    My2ndMech.chemfile = os.path.join(mechanism_dir, "C2_NOx_SRK.inp")
    # instruct the preprocessor to include the transport properties
    # only when the mechanism file contains all the transport data
    My2ndMech.preprocesstransportdata()
    # preprocess the 2nd mechanism files
    iError = My2ndMech.preprocess()
    print()
    if iError != 0:
        print(f"PreProcess: error encountered...code = {iError:d}")
        print(f"see the summary file {My2ndMech.summaryfile} for details")
        exit()
    else:
        print(Color.GREEN + "PreProcess success!!", end=Color.END)
        print("mechanism information:")
        print(f"number of elements = {My2ndMech.MM:d}")
        print(f"number of gas species = {My2ndMech.KK:d}")
        print(f"number of gas reactions = {My2ndMech.IIGas:d}")

    print(f"\nelement and species information of mechanism {My2ndMech.label}")
    print("=" * 50)
    # extract element symbols as a list
    elelist = My2ndMech.elementsymbols
    # get atomic masses as numpy 1D double array
    AWT = My2ndMech.AWT
    # print element information
    for k in range(len(elelist)):
        print(f"element # {k+1:3d}: {elelist[k]:16} mass = {AWT[k]:f}")

    print("=" * 50)
    # extract gas species symbols as a list
    specieslist = My2ndMech.speciessymbols
    # get species molecular masses as numpy 1D double array
    WT = My2ndMech.WT
    # print gas species information
    for k in range(len(specieslist)):
        print(f"species # {k+1:3d}: {specieslist[k]:16} mass = {WT[k]:f}")
    print("=" * 50)


if __name__ == "__main__":
    test_loadmechanism()

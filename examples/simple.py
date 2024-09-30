import os

import chemkin  # import PyChemkin

# create a Chemistry Set for GRI 3.0 mechanism in the data directory
mechanism_dir = os.path.join(chemkin.ansys_dir, "reaction", "data")
# set up mechanism file names
mech_file = os.path.join(mechanism_dir, "grimech30_chem.inp")
therm_file = os.path.join(mechanism_dir, "grimech30_thermo.dat")
tran_file = os.path.join(mechanism_dir, "grimech30_transport.dat")
# instantiate Chenistry Set 'GasMech'
GasMech = chemkin.Chemistry(
    chem=mech_file, therm=therm_file, tran=tran_file, label="GRI 3.0"
)
# pre-process the Chemistry Set
status = GasMech.preprocess()
# check preprocess status
if status != 0:
    # failed
    print(f"PreProcess: error encountered...code = {status:d}")
    print(f"see the summary file {GasMech.summaryfile} for details")
    exit()
# Create Mixture 'air' based on 'GasMech'
air = chemkin.Mixture(GasMech)
# set 'air' condition
# mixture pressure in [dynes/cm2]
air.pressure = 1.0 * chemkin.Patm
# mixture temperature in [K]
air.temperature = 300.0
# mixture composition in mole fractions
air.X = [("O2", 0.21), ("N2", 0.79)]
#
print(f"pressure    = {air.pressure/chemkin.Patm} [atm]")
print(f"temperature = {air.temperature} [K]")
# print the 'air' composition in mass fractions
air.listcomposition(mode="mass")
# get 'air' mixture density [g/cm3]
print(f"the mixture density   = {air.RHO} [g/cm3]")
# get 'air' mixture viscosity [g/cm-sec] or [poise]
print(f"the mixture viscosity = {air.mixtureviscosity()*100.0} [cP]")

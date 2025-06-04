Mixture
=======

*Mixture* is the "core" token in **PyChemkin**. It represents a gas-phase mixture by storing
its pressure, temperature, and species composition. *Stream* is an extended "Class" of *Mixture*
with the additional attribute of "mass/volumetric flow rate". *Mixture* and *Stream* include utilities
that can be used to evaluate the mixture properties (density, mixture enthalpy, mixture viscosity, ...)
and the reaction rates. There are also methods to manipulate the *Mixture*/*Stream* such as merging two mixtures
adiabatically.

The examples demonstrate how to instantiate the *Mixture*/*Stream* objects, how to extract information and properties
of a *Mixture*/*Stream*, and how to manipulate *Mixture*/*Stream* objects.

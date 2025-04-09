Getting Started
===============

Prerequisites
-------------
- `Ansys Chemkin`_ release 2025R2 or newer and a valid license
- `Python`_ 3.9 or newer
- `Numpy`_ 1.14.0 or newer
- `PyYAML`_ 6.0 or newer
- `matplotlib`_ for running examples and tutorials

.. _Ansys Chemkin: https://www.ansys.com/products/fluids/ansys-chemkin-pro

.. _Python: https://www.python.org/downloads/windows/

.. _Numpy: https://numpy.org/install

.. _PyYAML: https://pypi.org/project/PyYAML/

.. _matplotlib: https://matplotlib.org/stable/install/index.html

Installation
------------
- Installing the ``ansys-chemkin`` package

   **PyChemkin** is available for download as the ``ansys-chemkin`` package from the **PyAnsys** GitHub repository.
   Download the source codes and build the wheelhouse on local machine with the ``flit_core`` tools

   ::

      python -m build


   ``pip`` is the preferred installation method

   ::

      pip install dist\ansys_chemkin-*.whl


- Verifying the installation

  Invoke the Python interpreter interface from Windows' command prompt
  and try to import the ``ansys-chemkin`` package as

   ::

      >>> import ansys.chemkin


  - **PyChemkin** is correctly installed if python returns

   ::

      Chemkin version number = xxx


  - It is likely there is no Ansys product installed locally if python does not return anything

   ::

      >>>


  - The local Ansys Chemkin version needs to be updated to at least version *2025 Release 2* if Python returns

   ::

      ** PyChemkin does not support Chemkin versions older than 2025 R2



.. note:: A valid Ansys license is still required to run **PyChemkin** after the installation.

Example
-------
- Creating **PyChemkin** project
   The structure of a **PyChemkin** project follows the basic workflow of setting up and running a reactor model in Ansys Chemkin GUI:
      1. Create the *Chemistry Set* which includes the mechanism, the thermodynamic data, and/or the transport data.
      2. Pre-process the *Chemistry Set*.
      3. Create a *Mixture* based on the *Chemistry Set*.
      4. Specify the *Mixture* properties such as temperature, pressure, volume (optional), and species composition.
      5. Instantiate the *Reactor* by using the *Mixture* created.
      6. Set up the simulation:
         - reactor properties that are not provided by the initial *Mixture*.
         - solver parameters such as tolerances and solver timestep size.
         - solution saving controls such solution saving interval and adaptive solution saving.
      7. Run the simulation.
      8. Process the solution.

   Here is a **PyChemkin** project to compute the density of mixture ``air``

   .. code-block:: python
      :caption: A simple **PyChemkin** example
      :linenos:

         import os

         # import PyChemkin
         import ansys.chemkin as chemkin

         # create a Chemistry Set for GRI 3.0 mechanism in the data directory
         mech_dir = os.path.join(chemkin.ansys_dir, "reaction", "data")
         # set up mechanism file names
         mech_file = os.path.join(mech_dir, "grimech30_chem.inp")
         therm_file = os.path.join(mech_dir, "grimech30_thermo.dat")
         tran_file = os.path.join(mech_dir, "grimech30_transport.dat")
         # instantiate Chemistry Set 'GasMech'
         GasMech = chemkin.Chemistry(chem=mech_file, therm=therm_file,  tran=tran_file,  label='GRI 3.0')
         # pre-process the Chemistry Set
         status = GasMech.preprocess()
         # check preprocess status
         if status != 0:
             # failed
             print(f'PreProcess: error encountered...code = {status:d}')
             print(f'see the summary file {GasMech.summaryfile} for details')
             exit()
         # Create Mixture 'air' based on 'GasMech'
         air = chemkin.Mixture(GasMech)
         # set 'air' condition
         # mixture pressure in [dynes/cm2]
         air.pressure = 1.0 * chemkin.Patm
         # mixture temperature in [K]
         air.temperature = 300.0
         # mixture composition in mole fractions
         air.X = [('O2', 0.21), ('N2', 0.79)]
         #
         print(f"pressure    = {air.pressure/chemkin.Patm} [atm]")
         print(f"temperature = {air.temperature} [K]")
         # print the 'air' composition in mass fractions
         air.list_composition(mode='mass')
         # get 'air' mixture density [g/cm3]
         print(f"the mixture density = {air.RHO} [g/cm3]")

Support
-------
- Getting Help
   - help within **PyChemkin**
      - use ``ansys.chemkin.help()`` to get general information about a topic, for instance, 'ignition'.
      - use ``ansys.chemkin.keywordhints()`` to get description and syntax of a *reactor* keyword.
      - use ``ansys.chemkin.phrasehints()`` to get a list of *reactor* keywords related to a phrase.

   - online help
      `Ansys Learning Hub`_

   - technical support
      Contact Ansys account representatives


- Documents
   - `Chemkin Documents`_


.. _Ansys Learning Hub: https://learninghub.ansys.com/pages/17/home-page

.. _Chemkin Documents: https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?pn=Chemkin&pid=ChemkinPro&lang=en/


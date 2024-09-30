Getting Started
===============

- Installing the ``ansys-chemkin`` package

   **PyChemkin** is available for download as the ``ansys-chemkin`` package on `Test PyPI`_.

   ``pip`` is the preferred installation method

   ::

      python -m pip install -i https://test.pypi.org/simple/ ansys-chemkin


.. note:: Please refer to the `Prerequisites`_ for all required Python extensions to install/run **PyChemkin**.


- Verifying the installation

  Invoke the Python interpreter interface from Windows' command prompt
  and try to import the ``ansys-chemkin`` package as

   ::

      >>> import chemkin


  - **PyChemkin** is correctly installed if Python returns

   ::

      chemkin version number = xxx


  - It is likely there is no Ansys product installed locally if Python does not return anything

   ::

      >>>


  - The local Ansys Chemkin version needs to be updated to at least version *2025 Release 1* if Python returns

   ::

      ** PyChemkin does not support Chemkin versions older than 2025R1



.. note:: A valid Ansys license is still required to run **PyChemkin** after the installation.


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

         # import PyChemkin
         import chemkin
         # create a Chemistry Set for GRI 3.0 mechanism in the data directory
         mech_dir = chemkin.ansys_dir + '\\reaction\\data'
         # set up mechanism file names
         mech_file = mech_dir+'\\grimech30_chem.inp'
         therm_file = mech_dir+'\\grimech30_thermo.dat'
         tran_file = mech_dir+'\\grimech30_transport.dat'
         # instantiate Chenistry Set 'GasMech'
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
         air.listcomposition(mode='mass')
         # get 'air' mixture density [g/cm3]
         print(f"the mixture density = {air.RHO} [g/cm3]")


- Getting Help
   - help within **PyChemkin**
      - use ``chemkin.help()`` to get general information about a topic, for instance, 'ignition'.
      - use ``chemkin.keywordhints()`` to get description and syntax of a *reactor* keyword.
      - use ``chemkin.phrasehints()`` to get a list of *reactor* keywords related to a phrase.

   - online help
      `Ansys Learning Hub`_

   - technical support
      Contact Ansys account representatives


- Documents
   - `Chemkin Documents`_


.. _Test PyPI: https://test.pypi.org/project/ansys-chemkin/

.. _Prerequisites: Prerequisites.rst

.. _Ansys Learning Hub: https://learninghub.ansys.com/pages/17/home-page

.. _Chemkin Documents: https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?pn=Chemkin&pid=ChemkinPro&lang=en/


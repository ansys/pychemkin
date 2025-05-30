User guide
==========

Introduction
------------

PyChemkin inherits all Ansys Chemkin functionalities, providing Pythonic interfaces to Chemkin utilities and reactor models. In addition, PyChemkin introduces a hierarchy of four key objects to enhance user experiences within the Python framework:

- **Chemistry set:** A collection of utilities that handle chemistry and species properties.
- **Mixture:** A basic element representing a gas mixture that can be manipulated and transformed.
- **Stream:** A mixture object with an additional property of mass flow rate.
- **Reactor:** An instance of a Chemkin reactor model that transforms the initial mixture into another mixture.

The following figure shows PyChemkin's mixture-centric concept. It illustrates the basic types of operations applicable to a mixture, which inherits all properties of a chemistry set. A mixture can be combined with a constraint, equilibrium (with constraints), or a process (by a reactor model).

.. image:: _static/mixture_concept.png
   :width: 600
   :alt: Mixture-centric concept in PyChemkin

Workflow for a basic reactor model
----------------------------------

The workflow for setting up and running a basic reactor model in PyChemkin is the same as in the Ansys Chemkin GUI:

#. Create a chemistry set, which includes the mechanism, thermodynamic data, and/or transport data.
#. Preprocess the chemistry set.
#. Create a mixture or a stream based on the chemistry set.
#. Specify mixture or stream properties, such as temperature, pressure, volume (optional), and species composition.
#. Instantiate the reactor using the created mixture.
#. Set up the simulation:

   - Specify reactor properties not provided by the initial mixture or stream, such as heat loss rate to the surroundings and end time.
   - Specify solver parameters, such as tolerances and solver timestep size.
   - Specify saving controls, such as the solution saving interval and adaptive solution saving.

#. Run the simulation.
#. Process the solution.
#. Customize the output if necessary and generate plots.

Code example
------------

This code example computes the density of a mixture named ``air``:

.. code-block:: python
   :linenos:

      import os

      # import PyChemkin
      import ansys.chemkin as chemkin

      # create a chemistry set for the GRI 3.0 mechanism in the data directory
      mech_dir = os.path.join(chemkin.ansys_dir, "reaction", "data")
      # set up mechanism file names
      mech_file = os.path.join(mech_dir, "grimech30_chem.inp")
      therm_file = os.path.join(mech_dir, "grimech30_thermo.dat")
      tran_file = os.path.join(mech_dir, "grimech30_transport.dat")
      # instantiate a chemistry set named 'GasMech'
      GasMech = chemkin.Chemistry(chem=mech_file, therm=therm_file,  tran=tran_file,  label='GRI 3.0')
      # preprocess the chemistry set
      status = GasMech.preprocess()
      # check preprocess status
      if status != 0:
          # failed
          print(f'Preprocessing: Error encountered. Code = {status:d}.')
          print(f'See the summary file {GasMech.summaryfile} for details.')
          exit()
      # create a mixture named 'air' based on the 'GasMech' chemistry set
      air = chemkin.Mixture(GasMech)
      # set 'air' condition
      # mixture pressure in [dynes/cm2]
      air.pressure = 1.0 * chemkin.Patm
      # mixture temperature in [K]
      air.temperature = 300.0
      # mixture composition in mole fractions
      air.X = [('O2', 0.21), ('N2', 0.79)]
      #
      print(f"Pressure    = {air.pressure/chemkin.Patm} [atm].")
      print(f"Temperature = {air.temperature} [K].")
      # print the 'air' composition in mass fractions
      air.list_composition(mode='mass')
      # get 'air' mixture density [g/cm3]
      print(f"Mixture density = {air.RHO} [g/cm3].")

For more examples, see :ref:`examples`.

Support
-------

Pythonic methods
~~~~~~~~~~~~~~~~~

To get help within PyChemkin, use these Pythonic methods:

- ``ansys.chemkin.help()``: Get general information about a topic, such as ``ignition``.
- ``ansys.chemkin.keywordhints()``: Get the description and syntax of a reactor keyword.
- ``ansys.chemkin.phrasehints()``: Get a list of reactor keywords related to a phrase.

Product training and documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For product training and documentation, see these Ansys resources:

- `Ansys Learning Hub`_
- `Chemkin product documentation`_


.. _Ansys Learning Hub: https://learninghub.ansys.com/pages/17/home-page
.. _Chemkin product documentation: https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?pn=Chemkin&pid=ChemkinPro&lang=en/
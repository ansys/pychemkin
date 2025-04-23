Introduction
============

What is **PyChemkin**? **PyChemkin**, in its core, inherits all Chemkin functionalities. However, in addition to providing a collection of pythonic interfaces to Chemkin utilities and reactor models, **PyChemkin** adopts a few new concepts to deliver enhanced user experiences within the Python framework.

**PyChemkin** introduces a hierarchy of three key objects:
   1. the *Chemistry Set* object: collection of utilities that handle chemistry and species properties.
   2. the *Mixture* object: the *core* concept in **PyChemkin** representing a gas mixture and the basic element that can be manipulated and transformed.
   3. the *Stream* object: a *Mixture* object with the additional property of mass flow rate.
   4. the *Reactor* object: instance of a Chemkin reactor model that transforms the initial *Mixture* to another *Mixture*.

   The figure below is a schematic of **PyChemkin**'s *Mixture*-centric concept showing the basic types of operations applicable to the *Mixture* object in **PyChemkin**: create (by inheriting all properties of a *Chemistry Set*), combine/mix (with constraint), equilibrium (with constraints), and process (by a reactor model)

   .. image:: _static/mixture_concept.png
      :width: 600
      :alt: Mixture-centric Concept in PyChemkin

The structure of a **PyChemkin** project follows the basic workflow of setting up and running a reactor model in Ansys Chemkin GUI:
   1. Create the *Chemistry Set* which includes the mechanism, the thermodynamic data, and/or the transport data.
   2. Pre-process the *Chemistry Set*.
   3. Create a *Mixture* or a *Stream* based on the *Chemistry Set*.
   4. Specify the *Mixture/Stream* properties such as temperature, pressure, volume (optional), and species composition.
   5. Instantiate the *Reactor* by using the *Mixture* created.
   6. Set up the simulation.

        - reactor properties that are not provided by the initial *Mixture/Stream*. for example, heat loss rate to the surroundings and end time.
        - solver parameters such as tolerances and solver timestep size.
        - solution saving controls such as solution saving interval and adaptive solution saving.
   7. Run the simulation.
   8. Process the solution.
   9. Customize Output and/or generate plots.

The tutorials introduce the essential operations in **PyChemkin** with step by step instructions. The tutorials are arranged from basic task such as loading and preprocessing the *Chemistry Set* to more complex simulation project such as running an ignition delay parameter study. It is recommended to go through them in the order listed.


.. note::
   Look up `Chemkin Documents`_  for detailed information about Chemkin.


.. _Chemkin Documents: https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?pn=Chemkin&pid=ChemkinPro&lang=en/

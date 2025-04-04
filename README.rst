PyChemkin
=========

|

|pyansys| |python| |MIT|

.. |pyansys| image:: https://img.shields.io/badge/Py-Ansys-ffc107.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAIAAACQkWg2AAABDklEQVQ4jWNgoDfg5mD8vE7q/3bpVyskbW0sMRUwofHD7Dh5OBkZGBgW7/3W2tZpa2tLQEOyOzeEsfumlK2tbVpaGj4N6jIs1lpsDAwMJ278sveMY2BgCA0NFRISwqkhyQ1q/Nyd3zg4OBgYGNjZ2ePi4rB5loGBhZnhxTLJ/9ulv26Q4uVk1NXV/f///////69du4Zdg78lx//t0v+3S88rFISInD59GqIH2esIJ8G9O2/XVwhjzpw5EAam1xkkBJn/bJX+v1365hxxuCAfH9+3b9/+////48cPuNehNsS7cDEzMTAwMMzb+Q2u4dOnT2vWrMHu9ZtzxP9vl/69RVpCkBlZ3N7enoDXBwEAAA+YYitOilMVAAAAAElFTkSuQmCC
   :target: https://docs.pyansys.com/
   :alt: PyAnsys

.. |python| image:: https://img.shields.io/pypi/pyversions/ansys-geometry-core?logo=pypi
   :target: https://pypi.org/project/ansys-geometry-core/
   :alt: Python

.. |pypi| image:: https://img.shields.io/pypi/v/ansys-geometry-core.svg?logo=python&logoColor=white&label=PyPI
   :target: https://pypi.org/project/ansys-geometry-core
   :alt: PyPI

.. |MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/blog/license/mit
   :alt: MIT

.. contents::

Overview
--------

PyChemkin (Ansys-chemkin package) provides pythonic access to Ansys Chemkin. It facilitates programmatic customization of Chemkin simulation workflow within the Python ecosystem and permits access to Chemkin property and rate utilities as well as selected reactor models:

* Process Chemkin-compatible gas-phase mechanisms
* Evaluate species and mixture thermodynamic and transport properties
* Compute reaction rate of progress and species rate of production (ROP)
* Combine gas mixtures isothermally or adiabatically
* Find the equilibrium state of a gas mixture
* Run gas-phase batch reactor models

For more information on Chemkin, see the Ansys Chemkin page on the Ansys website.

Installation
^^^^^^^^^^^^

 ``pip`` is the preferred installation method. You can use `pip <https://pypi.org/project/pip/>`_ to install PyChemkin.

.. code:: bash

    pip install ansys-chemkin


.. note:: Please refer to the `Prerequisites`_ for all required Python extensions to install/run **PyChemkin**.

Verifying the installation
""""""""""""""""""""""""""

  Invoke the Python interpreter interface from Windows' command prompt
  and try to import the ``ansys-chemkin`` package as

.. code:: python

    import ansys.chemkin



* **PyChemkin** is correctly installed if python returns

.. code:: python

    Chemkin version number = xxx


* It is likely there is no Ansys product installed locally if python does not return anything.

* The local Ansys Chemkin version needs to be updated to at least version *2025 Release 2* if Python returns

.. code:: python

    ** PyChemkin does not support Chemkin versions older than 2025 R2

.. note:: A valid Ansys license is still required to run **PyChemkin** after the installation.

For more information, see `Getting Started`_.

Basic usage
^^^^^^^^^^^

Here is a **PyChemkin** project to compute the density of mixture ``air``. This code shows how to import
PyChemkin (Ansys Chemkin) and use some basic capabilities:

.. code:: python

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

For comprehensive usage information, see the Tutorials in the PyChemkin documentation.

Documentation and issues
^^^^^^^^^^^^^^^^^^^^^^^^
Documentation for the latest stable release of PyChemkin is hosted at PyChemkin documentation.

In the upper right corner of the documentation's title bar, there is an option for switching from
viewing the documentation for the latest stable release to viewing the documentation for the
development version or previously released versions.

User manuals and tutorials for Chemkin can be found at `Chemkin Documents`_.

On the PyAnsys Chemkin Issues page,
you can create issues to report bugs and request new features. On the PyChemkin Discussions
page or the `Discussions <https://discuss.ansys.com/>`_
page on the Ansys Developer portal, you can post questions, share ideas, and get community feedback.

To reach the project support team, email `pyansys.core@ansys.com <mailto:pyansys.core@ansys.com>`_.

.. LINKS AND REFERENCES
.. _Prerequisites: ./doc/source/Prerequisites.rst
.. _Getting Started: ./doc/source/Get_Start.rst
.. _Chemkin Documents: https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?pn=Chemkin&pid=ChemkinPro&lang=en/

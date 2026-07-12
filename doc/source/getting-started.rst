Getting started
===============

Install prerequisites
---------------------

- `Ansys Chemkin`_ 2025 R2 or later with a valid license
- `Python`_ 3.9 or later
- `NumPy`_ 1.14.0 or later
- `PyYAML`_ 6.0 or later
- `Matplotlib`_ to run examples

.. note:: Using the latest Ansys Chemkin version is highly recommended.

.. _Python: https://www.python.org/downloads/windows/
.. _NumPy: https://numpy.org/install
.. _PyYAML: https://pypi.org/project/PyYAML/
.. _Matplotlib: https://matplotlib.org/stable/install/index.html
.. _Flit: https://flit.pypa.io/en/stable/
.. _pip: https://pip.pypa.io/en/stable/index.html

Install PyChemkin
-----------------

1. Install PyChemkin.

   Download the ``ansys-chemkin-core`` package from the PyAnsys GitHub repository. Build the wheelhouse locally using `Flit`_:

   ::

      python -m build

   Install the package using `pip`_:

   ::

      pip install dist\ansys_chemkin-*.whl

   Alternatively, you can install the package directly from the PyPI repository:

   ::

      pip install ansys-chemkin-core

2. Verify the installation.

   Open the Python interpreter from the Windows command prompt and import the
   ``ansys-chemkin-core`` package:

   ::

      >>> import ansys.chemkin.core

   If PyChemkin is installed correctly, Python displays a statement like this:

   ::

      Chemkin version number = xxx
      PyChemkin version number = x.x.x

   PyChemkin is probably not installed locally if Python displays nothing:

   ::

      >>>

   If Python displays the following statement, update the local Ansys Chemkin installation
   to 2025 R2 or later:

   ::

      PyChemkin does not support Chemkin versions older than 2025 R2

.. note::
   1. You must have a valid Ansys license to run PyChemkin after installation.

   2. PyChemkin will search for the latest version of the local Ansys Chemkin installation
      and use it automatically. Use the ``PYCK_CHEMKIN_VER`` environment variable to specify
      the desired local Ansys Chemkin installation. For example, set ``PYCK_CHEMKIN_VER="261"``
      to force PyChemkin to use Ansys Chemkin 2026 R1.

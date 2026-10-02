Getting started
===============

Prerequisites
-------------

Before installing PyChemkin, install the following:

* `Ansys Chemkin`_ 2025 R2 or later, with a valid license.
* `Python`_ 3.10 or later, and earlier than 4.0.

PyChemkin declares its Python dependencies in ``pyproject.toml``. They are
installed automatically when PyChemkin is installed with ``pip``.

.. note:: Using the latest Ansys Chemkin version is recommended.

.. _Ansys Chemkin: https://www.ansys.com/products/fluids/ansys-chemkin
.. _Python: https://www.python.org/downloads/

Configure Ansys Chemkin
-----------------------

PyChemkin selects the newest supported local Chemkin installation. To select a
specific installation, define its ``ANSYSxxx_DIR`` environment variable, where
``xxx`` is the Chemkin release number. The value must be the installation's
``ANSYS`` directory. For example, on Windows:

.. code-block:: powershell

   $env:ANSYS261_DIR = "C:\Program Files\ANSYS Inc\v261\ANSYS"

When there are multiple ``ANSYSxxx_DIR`` environment variables defined, PyChemkin selects
the newest one by default. Use the ``PYCK_CHEMKIN_VER`` environment variable to specify
the desired local Ansys Chemkin installation. For example, set ``PYCK_CHEMKIN_VER="261"``
to force PyChemkin to use Ansys Chemkin 2026 R1.

On Linux, define the environment variable in the shell before starting Python.
You also need to source the ``chemkin_setup.ksh`` or ``chemkin_setup.csh`` script
provided in the Chemkin ``bin`` directory when using Linux.

Install PyChemkin
-----------------

Install the published package from PyPI with:

.. code-block:: console

   python -m pip install ansys-chemkin-core

To install the current source tree for development, run these commands from
the repository root:

.. code-block:: console

   python -m pip install --upgrade pip
   python -m pip install -e .

Build and install a wheel
--------------------------

To build a wheel from the source tree, install the build frontend and run it
from the repository root:

.. code-block:: console

   python -m pip install build
   python -m build

The wheel and source distribution are written to ``dist``. Install the wheel
with:

.. code-block:: console

   python -m pip install dist\ansys_chemkin_core-*.whl

Verify the installation
-----------------------

Start Python and import the package:

.. code-block:: pycon

   >>> import ansys.chemkin.core

The import initializes the native Chemkin library and reports the Chemkin and
PyChemkin versions. If initialization fails, verify that the selected Chemkin
installation exists, the ``ANSYSxxx_DIR`` variable is correct, and a valid
license is available.

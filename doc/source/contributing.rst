Contribute
==========

Thank you for your interest in contributing to PyChemkin. Contributions for making the project
better can include fixing bugs, adding new features, and improving the documentation.

Overall guidance on contributing to a PyAnsys library appears in the
`Contributing`_ topic in the *PyAnsys developer's guide*. Ensure that you
are thoroughly familiar with this guide before attempting to contribute to
PyChemkin.

The following contribution information is specific to PyChemkin.

- Clone the repository.
    To clone and install the latest PyChemkin release in development mode, run
    these commands:

    ::

        git clone https://github.com/ansys/pychemkin/
        cd pychemkin
        python -m pip install --upgrade pip
        pip install -e .

- Follow the code style.
    PyChemkin follows the PEP 8 standard as described in `PEP 8`_ in the *PyAnsys developer’s guide*
    and uses ``pre-commit`` for style checking.

    To ensure your code meets minimum code styling standards, run these commands:

    ::

        pip install pre-commit
        pre-commit run --all-files

    You can also install this as a pre-commit hook by running this command:

    ::

        pre-commit install

- Run the tests.
    Prior to running the tests, you must run this command to install the test dependencies:

    ::

        pip install - .tests

    To run the tests, navigate to the root directory of the repository and run this command:

    ::

        pytest


- Build the documentation.
    Prior to building the documentation, you must run this command to install the documentation dependencies:

    ::

        pip install -e .doc

    To build the documentation, run the following commands:

    ::

        cd doc

    - On Linux:

      ::

          make html

    - On Windows:

      ::

          ./make.bat html


The documentation is built in the ``docs/_build/html`` directory.

.. _Contributing: https://dev.docs.pyansys.com/how-to/contributing.html
.. _PEP 8: https://dev.docs.pyansys.com/coding-style/pep8.html#pep-8

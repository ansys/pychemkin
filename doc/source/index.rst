.. title:: PyChemkin

.. figure:: _static/logo/PyChemkin.png
    :align: center
    :width: 640px

PyChemkin (the Ansys-chemkin package) provides pythonic access to Ansys Chemkin-CFD-API. It facilitates programmatic customization
of Chemkin simulation workflow within the Python ecosystem and permits access to Chemkin property and rate utilities as well as
selected reactor models:

* Process Chemkin-compatible gas-phase mechanisms
* Evaluate species and mixture thermodynamic and transport properties
* Compute reaction rate of progress and species rate of production (ROP)
* Combine gas mixtures isothermally or adiabatically
* Find the equilibrium state of a gas mixture
* Run gas-phase batch, plug-flow, and perfectly-stirred reactor models
* Calculate the laminar flame speed of a combustible mixture
* Create and solve steady-state reactor network


.. grid::   1 2 3 3
    :gutter: 1 2 3 3
    :padding: 1 2 3 3

    .. grid-item-card:: Getting started :fa:`person-running`
        :link: getting_started
        :link-type: doc

        Learn how to install and verify the PyChemkin package.

    .. grid-item-card:: Introduction :fa:`book-open-reader`
        :link: introduction
        :link-type: doc

        Short introduction of the unique feaures of PyChemkin.

    .. grid-item-card:: User guide :fa:`book-open-reader`
        :link: https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/prod_page.html?pn=Chemkin&pid=ChemkinPro&lang=en/
        :link-type: url

        Access Chemkin manuals for details about the parameters, the options,
        and the theory related to the properties and the reactor models.

    .. grid-item-card:: API reference :fa:`book-bookmark`
        :link: autoapi/index
        :link-type: doc

        Understand PyChemkin API endpoints, their capabilities,
        and how to interact with them programmatically.

    .. grid-item-card:: Examples :fa:`scroll`
        :link: auto_examples/index
        :link-type: doc

        Explore examples that show how to use PyChemkin utilities and
        reactor models to perform different types of simulations.

    .. grid-item-card:: Contribute :fa:`people-group`
        :link: contributing
        :link-type: doc

        Learn how to contribute to the PyChemkin codebase
        or documentation.

.. jinja:: main_toctree

    .. toctree::
       :hidden:
       :maxdepth: 3

       getting_started
       introduction
       Tutorials <auto_examples/index.rst>
       contributing


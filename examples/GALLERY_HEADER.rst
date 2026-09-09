.. _examples:

Examples
========

This section shows how to use PyChemkin utilities and reactor models to perform different types of simulations.
Early examples demonstrate how to perform essential operations, such as loading and preprocessing the chemistry set.
Later examples explain how to tackle more complex simulations, such as running an ignition delay parameter study.
Thus, you should go through the examples in the order listed.

.. note::
    Some examples require mechanism and data files that are not included through the
    download links provided with their descriptions. These supplementary files
    are available in the ``examples/data`` folder of the *PyChemkin* repository.
    The corresponding Python source files are available in the ``examples`` folder.
    Make sure the full file paths to the mechanism and the data files are correct
    in the example scripts before running them.

.. note::
    When you use the Jupyter Notebook version of the example projects (``*.ipynb``), please
    remember to *close the project after you finish viewing and/or running the project*. When the project
    is open, it might hold on to your Ansys license. If you have too many Jupyter Notebook projects
    active, you might run out of your Ansys licenses. Closing the Jupyter Notebook project will
    release the license. If the *Jupyter Notebook web interface* is used to load the PyChemkin
    example projects, you can navigate to the top menu bar and select ``File > Close and Shut Down Notebook``
    to close the project and release the license. If the *VS Code* is used, you can simply close the tab
    associated with the example.

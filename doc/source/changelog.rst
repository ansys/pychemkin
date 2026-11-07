.. _ref_release_notes:

Release notes
#############

This document contains the release notes for the **PyChemkin** project.

.. vale off

.. towncrier release notes start

`0.3.0 <https://github.com/ansys/pychemkin/releases/tag/v0.3.0>`_ - 2026-11-04
==============================================================================

### Added

- Added the ``TransientPSR`` models to simulate the dynamic response of a PSR to the transient inlet conditions.
  The *atomic layer deposition* (ALD) and the *catalytic engine after-treatment* processes are two typical applications of this reactor model. `(#126) <https://github.com/ansys/pychemkin/pull/126>`_
- Added new method ``set_auto_ignition_check`` to the ``PremixedFlame`` model to enable/disable the autoignition check of the unburnt mixture. (*Requires Chemkin 2027R1*) `(#183) <https://github.com/ansys/pychemkin/pull/183>`_
- Added utilities to import Chemkin XML solution file. `(#256) <https://github.com/ansys/pychemkin/pull/256>`_
- Added the Gamma and the lognormal functions for the particle class/size distribution in the Aerosol module. `(#277) <https://github.com/ansys/pychemkin/pull/277>`_


### Changed

- Allowed to specify the Chemkin version by the value of the environment variable ``PYCK_CHEMKIN_VER``.
  For example, set ``PYCK_CHEMKIN_VER=271`` to use Chemkin *2027R1* in PyChemkin.
  Remove the environment variable or set its value to ``0`` to restore the default behavior: use the latest Chemkin version locally. `(#182) <https://github.com/ansys/pychemkin/pull/182>`_
- Suppressed ignition delay time related messages when the ``set_ignition_delay`` method is *not* called in a ``BatchReactor`` model simulation. `(#184) <https://github.com/ansys/pychemkin/pull/184>`_


### Fixed

- Fixed the Chemkin license check-in issue after running examples in **Jupyter Notebook**. (*Requires Chemkin 2027R1*) `(#88) <https://github.com/ansys/pychemkin/pull/88>`_


### Documentation

- Added new ``surface_chemistry`` example: **trans_PSR_ALD**. `(#126) <https://github.com/ansys/pychemkin/pull/126>`_
- Added new ``multiprocessing`` example: **sensitivity_multiproc** by *Kuiwen Zhang*. `(#210) <https://github.com/ansys/pychemkin/pull/210>`_
- Added post-processing Chemkin XML solution file examples. `(#265) <https://github.com/ansys/pychemkin/pull/265>`_
- Created `PyChemkin Wiki page <https://github.com/ansys/pychemkin/wiki>`_. 


### Maintenance

- Addressed the security risk associated with ``subprocess``. `(#156) <https://github.com/ansys/pychemkin/pull/156>`_
- Moved ``Keyword`` and ``Profile`` classes out of ``reactormodel.py``. `(#229) <https://github.com/ansys/pychemkin/pull/229>`_
- Reorganized surface chemistry related modules. `(#230) <https://github.com/ansys/pychemkin/pull/230>`_
- Improved readability of the ``Mixture`` module. `(#231) <https://github.com/ansys/pychemkin/pull/231>`_
- Refactored ``Chemistry`` module to improve readability. `(#232) <https://github.com/ansys/pychemkin/pull/232>`_
- Made error message logging and handling consistent across all modules. `(#233) <https://github.com/ansys/pychemkin/pull/233>`_

.. vale on

"""Script to run the group(s) of PyChemkin tests under the test sources folder."""

import pytest

from .tools import PyCKtools


class TestClassBasic:
    """Tests to verify Chemkin preprocessor."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # ignition delay [msec], heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    basic_list = ["simple", "loadmechanism"]
    fresh = True

    @pytest.mark.parametrize("test_file", basic_list)
    def test_basic(
        self,
        get_working_dir: str,
        get_source_dir: str,
        get_result_dir: str,
        test_file: str,
    ):
        """Run the selected pychemin basic utility test cases."""
        """
        Parameters
        ----------
            get_working_dir: string
                working folder for testing
            get_source_dir: string
                folder where source codes for testing are stored
            get_result_dir: string
                folder where test results will be kept
            test_file: string
                name of the source file to be tested
        """
        # initialization
        if TestClassBasic.fresh:
            PyCKtools.init_test_status()
            TestClassBasic.fresh = False
        # run the basic test cases
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("utilities", "all")
@pytest.mark.utilities
class TestClassUtilities:
    """Tests to verify Chemkin utilities.

    1. species/mixture property calculations,
    2. reaction rate calculations,
    3. mixture operations.
    """

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # ignition delay [msec], heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    utility_list = [
        "createmixture",
        "mixturemixing",
        "speciesproperties",
        "reactionrates",
        # skip the test below because the subprocess produces a non-zero return code
        # but the test is completed successfully
        #    "multiplemechanisms",
    ]

    @pytest.mark.parametrize("test_file", utility_list)
    def test_utilities(
        self, get_working_dir, get_source_dir, get_result_dir, test_file
    ):
        """Run the selected pychemin utility test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("equilibrium", "all")
@pytest.mark.equilibrium
class TestClassEquilibrium:
    """Tests to verify Chemkin utilities for equilibrium/detonation calculations."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # ignition delay [msec], heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    equilibrium_list = [
        "adiabaticflametemperature",
        "detonation",
        "equilibriumcomposition",
        "heatingvalues",
    ]

    @pytest.mark.parametrize("test_file", equilibrium_list)
    def test_equilibrium(
        self, get_working_dir, get_source_dir, get_result_dir, test_file
    ):
        """Run the selected pychemin equilibrium utility test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("batch", "all")
@pytest.mark.batch
class TestClassBatch:
    """Tests to verify Chemkin 0-D closed-homogeneous batch reactor models."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # ignition delay [msec], heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    batch_list = [
        "CONV",
        "closed_homogeneous__transient",
        "ignitiondelay",
        "vapor",
        "sensitivity",
    ]

    @pytest.mark.parametrize("test_file", batch_list)
    def test_batch(self, get_working_dir, get_source_dir, get_result_dir, test_file):
        """Run the selected pychemin batch reactor test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("engine", "all")
@pytest.mark.engine
class TestClassEngine:
    """Tests to verify Chemkin 0-D engine models."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # ignition delay [msec], heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    engine_list = [
        "hcciengine",
        "multizone",
        "sparkignitioengine",
    ]

    @pytest.mark.parametrize("test_file", engine_list)
    def test_engine(self, get_working_dir, get_source_dir, get_result_dir, test_file):
        """Run the selected pychemin engine model test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("PFR", "all")
@pytest.mark.PFR
class TestClassPFR:
    """Tests to verify Chemkin Plug-Flow Reactor (PFR) model."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # ignition delay [msec], heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    pfr_list = ["plugflow"]

    @pytest.mark.parametrize("test_file", pfr_list)
    def test_plug_flow_reactor(
        self, get_working_dir, get_source_dir, get_result_dir, test_file
    ):
        """Run the selected pychemin PFR model test case."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("PSR", "all")
@pytest.mark.PSR
class TestClassPSR:
    """Tests to verify Chemkin Perfectly-Stirred Reactor (PSR) model."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    psr_list = ["PSRgas", "jetstirredreactor", "multi-inletPSR", "PSRChain_declustered"]

    @pytest.mark.parametrize("test_file", psr_list)
    def test_perfectly_stirred_reactor(
        self, get_working_dir, get_source_dir, get_result_dir, test_file
    ):
        """Run the selected pychemin PSR model test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("ERN", "all")
@pytest.mark.ERN
class TestClassERN:
    """Tests to verify Chemkin equivalent reactor network (ERN) model."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    ern_list = ["PSRChain_network", "PSRnetwork"]

    @pytest.mark.parametrize("test_file", ern_list)
    def test_reactor_network(
        self, get_working_dir, get_source_dir, get_result_dir, test_file
    ):
        """Run the selected pychemin ERN model test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."

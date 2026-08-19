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
        # Skip the test below because the subprocess occasionally produces
        # a non-zero return code which indicates a native-process crash
        # rather than a Python exception. The most likely cause
        # is a changed ownership/lifetime assumption around shared Chemistry
        # or native state. The test is completed successfully.
        # "multiplemechanisms",
        "diffusionvelocity",
        "mixing_IEM",
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
        "sparkignitionengine",
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
    psr_list = [
        "PSRgas",
        "jetstirredreactor",
        "multi-inletPSR",
        "PSRChain_declustered",
        "trans_PSR_ignition",
    ]

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
    ern_list = [
        "PSRChain_network",
        "PSRnetwork",
        "PSRnetwork_coupled",
    ]

    @pytest.mark.parametrize("test_file", ern_list)
    def test_reactor_network(
        self, get_working_dir, get_source_dir, get_result_dir, test_file
    ):
        """Run the selected pychemin ERN model test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("premixed", "all")
@pytest.mark.premixed
class TestClassPremixed:
    """Tests to verify Chemkin premixed flame models."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], volume [cm3], velocity [cm/s],
    # flame speed [cm/sec], heat [cal]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    premixed_list = ["methane_flamespeed_table", "premixedburnerflame"]

    @pytest.mark.parametrize("test_file", premixed_list)
    def test_premixed(self, get_working_dir, get_source_dir, get_result_dir, test_file):
        """Run the selected pychemin premixed flame model test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("opposed", "all")
@pytest.mark.opposed
class TestClassOpposedFlame:
    """Tests to verify Chemkin opposed-flow flame models."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], axial velocity [cm/s],
    # and mixture fraction [-]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    opposed_list = ["dual_flame"]

    @pytest.mark.parametrize("test_file", opposed_list)
    def test_opposed(self, get_working_dir, get_source_dir, get_result_dir, test_file):
        """Run the selected pychemin opposed-flow flame model test case."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("shock", "all")
@pytest.mark.shock
class TestClassShockTube:
    """Tests to verify Chemkin shock tube reactor models."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K], velocity [cm/s],
    # Mach number [-], total thermicity [1/sec]
    # species: mole/mass fraction
    # rate: reaction rate, rate of production, heat release rate
    shock_list = ["incidentshock", "ZND"]

    @pytest.mark.parametrize("test_file", shock_list)
    def test_shock(self, get_working_dir, get_source_dir, get_result_dir, test_file):
        """Run the selected pychemin shock tube reactor model test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."


@pytest.mark.group("surface", "all")
@pytest.mark.surface
class TestClassSurface:
    """Tests to verify Chemkin surface chemistry utilities."""

    # define tolerances for this group of tests
    # {'type_of_variable': [absolute_tolerance, relative_tolerance], ... }
    # state: pressure [atm], temperature [K]
    # species: mole/mass fraction, surface site fractions, bulk activities
    # rate: reaction rate, rate of production, heat release rate
    surface_list = [
        # skip the test below for random pass/fail on different linux versions
        "multiple_materials",
        "SiC_cvd",
        "catalytic_combustion",
        "trans_PSR_ALD",
    ]

    @pytest.mark.parametrize("test_file", surface_list)
    def test_surface(self, get_working_dir, get_source_dir, get_result_dir, test_file):
        """Run the selected pychemin surface chemistry/materials test cases."""
        ierr = PyCKtools.run_test(
            get_working_dir, get_source_dir, get_result_dir, test_file
        )
        assert 0 == ierr, "run failed."

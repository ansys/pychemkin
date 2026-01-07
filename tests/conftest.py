"""PyTest settings for PyChemkin test runs."""

import pytest


def pytest_addoption(parser):
    """Register custom command line options."""
    # select test groups to run
    parser.addoption(
        "--group",
        action="store",
        choices=(
            "all",
            "utilities",
            "equilibrium",
            "batch",
            "engine",
            "PFR",
            "PSR",
            "ERN",
            "",
        ),
        default="all",
        help="run test groups named in GROUPNAME.",
    )
    parser.addoption(
        "--source",
        action="store",
        default="integration_tests",
        help="specify the folder containing the test source files.",
    )
    parser.addoption(
        "--baseline",
        action="store",
        default="baseline",
        help="specify the folder containing the test baselines.",
    )
    parser.addoption(
        "--results",
        action="store",
        default="new_results",
        help="specify the folder containing the test results.",
    )
    parser.addoption(
        "--compare",
        action="store_true",
        help="compare results against the baselines.",
    )


def pytest_configure(config):
    """Register custom markers."""
    # configure marker: @pytest.mark.group(GROUPNAME) for each group class
    config.addinivalue_line("markers", "group(groupname): mark test groups to run.")
    config.option.python_files = ["test_pychemkin_*.py"]


def pytest_runtest_setup(item):
    """Skip PyChemkin test groups that is not specified by the '--group' option."""
    """
    By default, pytest will run all test groups (i.e., --group="all").
    """
    # run group if option "--group" is given and the group name appears
    # in the option argument
    groupnames = []
    for mark in item.iter_markers(name="group"):
        groupnames = mark.args
        if groupnames:
            if item.config.getoption("--group") not in groupnames:
                pytest.skip(f"test group {groupnames!r} skipped.")


def pytest_collection_modifyitems(items):
    """Set the order of test groups and ensure the result comparison runs last."""
    class_order = [
        "TestClassBasic",
        "TestClassUtilities",
        "TestClassEquilibrium",
        "TestClassBatch",
        "TestClassEngine",
        "TestClassPFR",
        "TestClassPSR",
        "TestClassERN",
        "TestCompareResults",
    ]
    class_mapping = {item: item.cls.__name__ for item in items}
    sorted_items = items.copy()
    # Iteratively move tests of each class to the end of the test queue
    for class_ in class_order:
        sorted_items = [it for it in sorted_items if class_mapping[it] != class_] + [
            it for it in sorted_items if class_mapping[it] == class_
        ]
    items[:] = sorted_items


@pytest.fixture
def get_compare(request):
    """Get comparison options."""
    return request.config.getoption("--compare")


@pytest.fixture
def get_source_dir(request):
    """Get alternative source code folder."""
    return request.config.getoption("--source")


@pytest.fixture
def get_baseline_dir(request):
    """Get alternative baseline folder."""
    return request.config.getoption("--baseline")


@pytest.fixture
def get_result_dir(request):
    """Get alternative results folder."""
    return request.config.getoption("--results")


@pytest.fixture
def get_working_dir(request):
    """Get alternative test working folder."""
    return request.config.rootdir

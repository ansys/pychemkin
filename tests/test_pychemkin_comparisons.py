"""Validate PyChemkin test results."""

from pathlib import Path

import pytest

from .tools import PyCKtools


class TestCompareResults:
    """Perform test validation."""

    line_length = 42

    def test_compare_baseline(
        self,
        get_working_dir,
        get_baseline_dir,
        get_result_dir,
        get_compare,
    ):
        """Compare current test results against the baselines."""
        if not get_compare:
            pytest.skip("no --compare option set.")
        # check the folders
        if PyCKtools.TARGET_FOLDER is None:
            new_working = get_working_dir
        else:
            new_working = PyCKtools.TARGET_FOLDER
        new_result_dir = str(Path(new_working) / get_result_dir)
        ierr = PyCKtools.check_folder(new_result_dir)
        assert ierr == 0, f"result folder {new_result_dir} not found."
        baseline_dir = str(Path(get_working_dir) / get_baseline_dir)
        ierr = PyCKtools.check_folder(baseline_dir)
        assert ierr == 0, f"baseline folder {baseline_dir} not found."
        # load the result file
        result_tag = ".result"
        baseline_tag = ".baseline"
        r_file_names = PyCKtools.get_file_names(new_result_dir)
        b_file_names = PyCKtools.get_file_names(baseline_dir)
        # create the comparison log file
        log = Path(new_working) / "compareresults.log"
        logf = log.open(mode="w+")
        count_all_files = 0
        count_missing = 0
        count_skipped = 0
        count_bad = 0
        for rf in r_file_names:
            rf_obj = Path(rf)
            test_name = rf_obj.stem
            extension = rf_obj.suffix
            logf.write(f"\nchecking file {rf}...\n")
            count_all_files += 1
            if extension != result_tag:
                logf.write("-" * self.line_length + "\n")
                logf.write(f"skip file {rf}.\n")
                logf.write("-" * self.line_length + "\n\n")
                count_skipped += 1
            else:
                # prepare for result comparisons
                this_result_file = str(Path(new_result_dir) / rf)
                this_result = PyCKtools.load_results(this_result_file)
                msg = str(test_name) + ": bad result file format."
                assert isinstance(this_result, dict), msg
                # load the baseline file
                base_name = test_name + baseline_tag
                if base_name in b_file_names:
                    this_baseline_file = str(Path(baseline_dir) / base_name)
                    # find corresponding baseline
                    this_baseline = PyCKtools.load_results(this_baseline_file)
                    msg = str(test_name) + ": trouble reading the baseline of test."
                    assert isinstance(this_baseline, dict), msg
                    # get tolerances from the baseline file
                    state_tol = this_baseline.get("tolerance-var", [1.0e-6, 1.0e-2])
                    species_tol = this_baseline.get("tolerance-frac", [1.0e-6, 1.0e-2])
                    rate_tol = this_baseline.get("tolerance-ROP", [1.0e-6, 1.0e-2])
                    # compare stored variables in the results
                    var_list = list(this_result.keys())
                    base_list = list(this_baseline.keys())
                    # baseline contains three more keys for the tolerances
                    msg = str(test_name) + ": result file mismatch."
                    assert PyCKtools.check_list_size(var_list, base_list, 3), msg
                    status = 0
                    # go through all stored variables
                    for var in var_list:
                        # get variable array
                        r_list = this_result.get(var, [0.0])
                        b_list = this_baseline.get(var, [0.0])
                        # check array sizes
                        msg = (
                            str(test_name)
                            + ": results have different number of entries."
                        )
                        assert PyCKtools.check_list_size(r_list, b_list, 0), msg
                        # set tolerances
                        atol, rtol = PyCKtools.get_tolerances(
                            var, state_tol, species_tol, rate_tol
                        )
                        # perform the comparison
                        ierr, bad, diff = PyCKtools.compare_list(
                            r_list, b_list, atol, rtol
                        )
                        #
                        status += ierr
                        if ierr > 0:
                            # there are differences out of tolerances
                            logf.write("-" * self.line_length + "\n")
                            logf.write(f"{test_name}::{var}\n")
                            for i in range(ierr):
                                id = bad[i]
                                msg = [
                                    "index = ",
                                    str(id),
                                    ",  ",
                                    "result value = ",
                                    str(r_list[id]),
                                    ",  ",
                                    "baseline value = ",
                                    str(b_list[id]),
                                    ",  ",
                                    "difference = ",
                                    str(diff[i]),
                                    "\n",
                                ]
                                logf.write(" ".join(msg))
                            logf.write("-" * self.line_length + "\n")
                    if ierr == 0:
                        # no significant difference found
                        logf.write("-" * self.line_length + "\n")
                        logf.write("OK.\n")
                        logf.write("-" * self.line_length + "\n\n")
                    else:
                        count_bad += 1
                    #
                else:
                    # no matching test baseline found
                    # assert 0, f"cannot find baseline for {test_name}."
                    logf.write("-" * self.line_length + "\n")
                    logf.write(f"{test_name}::baseline missing.\n")
                    logf.write("-" * self.line_length + "\n")
                    count_missing += 1
        logf.write("=" * self.line_length + "\n")
        msg = [
            "summary:\n",
            "        total files processed =",
            str(count_all_files),
            "\n",
            "                skipped files =",
            str(count_skipped),
            "\n",
            "        test missing baseline =",
            str(count_missing),
            "\n",
            "     test showing differences =",
            str(count_bad),
            "\n",
        ]
        logf.write(" ".join(msg))
        logf.write("=" * self.line_length + "\n")
        logf.close()

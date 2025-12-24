"""Additional tools for running PyChemkin tests."""

import ast
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np


class PyCKtools:
    """Utilities for testing pychemkin features."""

    FIRST_PASS = 0
    TARGET_FOLDER = ""

    @staticmethod
    def check_folder(thisfolder: str) -> int:
        """Verify the given folder exists."""
        """
        Parameters
        ----------
            thisfolder: string
                current folder name

        Returns
        -------
        error: integer
            error code

        """
        if not Path(thisfolder).exists():
            print(f"Error: folder {thisfolder} does not exist.")
            return 1
        return 0

    @staticmethod
    def create_folder(newfolder: str) -> int:
        """Create or clean up a folder."""
        """
        Parameters
        ----------
            newfolder: string
                new folder name

        Returns
        -------
            error: integer
                error code

        """
        new = Path(newfolder)
        if new.exists():
            # delete any existing files in this target folder
            print(f"Warning: All files in folder {newfolder:s} will be deleted")
            # first delete any existing files
            for oldfile in new.iterdir():
                try:
                    if oldfile.is_symlink():
                        oldfile.unlink()
                    elif oldfile.is_file():
                        oldfile.unlink()
                    elif oldfile.is_dir():
                        shutil.rmtree(str(oldfile))
                except Exception as e:
                    print(f"error: Fail to delete {oldfile.name:s}. Reason: {e:s}")
                    return 1
        else:
            # folder not exist, create new folder
            print(f"... Creating a new sub-folder {newfolder:s} to store run results")
            try:
                new.mkdir()
            except OSError:
                print(f"Error: Fail to create new folder {newfolder:s}")
                return 1
        return 0

    @staticmethod
    def run_test(
        root_dir: str, source_dir: str, result_dir: str, test_file: str
    ) -> int:
        """Run the given PyChemkin test."""
        """
        Parameters
        ----------
            root_dir: string
                rootdir of pytest
            source_dir: string
                the directory containing the PyChemkin source files
            result_dir: string
                the directory containing the test result files
            test_file: string
                name of the PyChemkin test

        Returns
        -------
            ReturnCode: integer
                return code from the subprocess run

        """
        #
        # find the working directory
        current_dir = str(Path.cwd())
        # set sources folder
        source_folder = str(Path(root_dir) / source_dir)
        # check if the source folder exists
        status = PyCKtools.check_folder(source_folder)
        assert 0 == status, "bad source folder name."
        # verify the output folder
        output_folder = str(Path(current_dir) / "outputs")
        # check if the output folder exists
        if PyCKtools.FIRST_PASS == 0:
            status = PyCKtools.create_folder(output_folder)
            assert 0 == status, "fail to get fresh output folder"
        # verify the temporary working folder
        new_working = str(Path(current_dir) / result_dir)
        if PyCKtools.TARGET_FOLDER == "":
            PyCKtools.TARGET_FOLDER = current_dir
        # check if the temporary working folder exists
        if PyCKtools.FIRST_PASS == 0:
            status = PyCKtools.create_folder(new_working)
            assert 0 == status, "fail to get fresh working folder"
            #
            PyCKtools.FIRST_PASS += 1
        # set the test source file
        runfile = test_file + ".py"
        frun = Path(source_folder) / runfile
        # set the test output file
        outfile = test_file + ".out"
        output_folder = Path(current_dir) / "outputs"
        f = output_folder / outfile
        fout = f.open(mode="w")
        # run test
        # change working directory
        os.chdir(new_working)
        try:
            results = subprocess.run(["python", str(frun)], stdout=fout, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command returned non-zero exit status: {e.returncode}")
            print(f"test = {test_file}")
            fout.close()
            # change working directory back
            os.chdir(current_dir)
            assert 0
        except subprocess.TimeoutExpired as e:
            print(f"Command timed out: {e.timeout}")
            print(f"test = {test_file}")
            fout.close()
            # change working directory back
            os.chdir(current_dir)
            assert 0
        except FileNotFoundError as e:
            print(f"Command not found: {e}")
            print(f"test = {test_file}")
            fout.close()
            # change working directory back
            os.chdir(current_dir)
            assert 0
        except OSError as e:
            print(f"OS error occurred: {e}")
            print(f"test = {test_file}")
            fout.close()
            # change working directory back
            os.chdir(current_dir)
            assert 0
        # close the solution output file
        fout.close()
        # clean up unimportant output files
        new = Path(new_working)
        for out_files in new.glob("*.out"):
            Path.unlink(out_files)
        for asc_files in new.glob("*.asc"):
            Path.unlink(asc_files)
        for inp_files in new.glob("*.inp"):
            Path.unlink(inp_files)
        # return the return code from the subprocess run
        # change working directory
        os.chdir(current_dir)
        #
        return results.returncode

    @staticmethod
    def load_results(pyck_result_file: str) -> dict:
        """Read PyChemkin result data from the text file."""
        """
        Parameters
        ----------
            pyck_result_file: string
                name of the result file from simulation

        Returns
        -------
            pyck_result: dict
                PyChemkin test result stored in a dictionary

        """
        with Path.open(pyck_result_file) as f:
            data = f.read()
        pyck_result = ast.literal_eval(data)
        return pyck_result

    @staticmethod
    def check_list_size(list1: list, list2: list, expected_diff: int = 0) -> bool:
        """Verify that the two lists are of the same size."""
        """
        Parameters
        ----------
            list1: List
                the first list of which the size to be compared
            list2: List
                the second list of which the size to be compared
            expected_diff: integer
                the expected length difference of the two lists

        Returns
        -------
            different: bool
                if the two lists have expected length difference

        """
        return abs(len(list1) - len(list2)) == expected_diff

    @staticmethod
    def get_file_names(folder_path: str) -> list[str]:
        """Get a list of all file names in the specified folder."""
        """
        Parameters
        ----------
            folder_path: string
                absolute path to the folder

        Returns
        -------
            list_files: list of strings
                names of actual files in the folder

        """
        file_names = []
        for entry in Path(folder_path).iterdir():
            if entry.is_file():
                file_names.append(entry.name)
        return file_names

    @staticmethod
    def get_tolerances(
        tolerance_name: str,
        state_tol: list[float],
        species_tol: list[float],
        rate_tol: list[float],
    ) -> tuple[float, float]:
        """Find the tolerance values according to the variable type."""
        """
        Parameters
        ----------
            tolerance_name: string
                type of variables to be compared
            state_tol: list of double
                tolerance value for state variables
            species_tol: list of double
                tolerance for species fractions
            rate_tol: list of double
                tolerance for reaction rates

        Returns
        -------
            atol: double
                absolute tolerance to be used
            rtol: double
                relative tolerance to be used

        """
        if "species" in tolerance_name:
            atol = species_tol[0]
            rtol = species_tol[1]
        elif "rate" in tolerance_name:
            atol = rate_tol[0]
            rtol = rate_tol[1]
        else:
            atol = state_tol[0]
            rtol = state_tol[1]
        return atol, rtol

    @staticmethod
    def compare_list(
        r_list: list[int | float],
        b_list: list[int | float],
        atol: float,
        rtol: float,
    ) -> tuple[int, list, list]:
        """Compare the values in the two int or float lists of the same size."""
        """
        Parameters
        ----------
            r_list: list of integers or doubles
                list of variable values from the new solution
            b_list: list of integers or doubles
                list of variable values from the baseline
            atol: double
                absolute tolerance
            rtol: double
                relative tolerance

        Returns
        -------
            status: int
                comparison outcome
            bad_id: list[int]
                list of variable index that the difference > the tolerances
            diff: list[int or float]
                the difference between the values from the two lists

        """
        l_size = len(r_list)
        bad_id = []
        diff = []
        # go down the lists
        for i in range(l_size):
            r_value = r_list[i]
            b_value = b_list[i]
            a_diff = b_value - r_value
            if abs(a_diff) > atol:
                # check relative tolerance
                r_diff = a_diff - b_value * rtol
                if r_diff > 0.0:
                    bad_id.append(i)
                    if np.isclose(b_value, 0.0, atol=1.0e-9):
                        diff.append(abs(a_diff))
                    else:
                        diff.append(abs(a_diff / b_value))
        # check differences
        ierr = len(bad_id)
        return ierr, bad_id, diff

    def init_test_status():
        """Initialize the test status values before running the tests."""
        PyCKtools.TARGET_FOLDER = ""
        PyCKtools.FIRST_PASS = 0

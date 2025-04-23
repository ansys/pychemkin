"""
Additional tools for running PyChemkin tests.
"""
import ast
import glob
import os
import shutil
import subprocess

import numpy as np


class PyCKtools:
    FIRST_PASS = 0
    TARGET_FOLDER = None

    def check_folder(thisfolder):
        """
        Verify the given folder exists.
        """
        if not os.path.exists(thisfolder):
            print(f"Error: folder {thisfolder} does not exist.")
            return 1
        return 0

    def create_folder(newfolder):
        """
        Create or clean up a folder.
        """
        if os.path.exists(newfolder):
            # delete any existing files in this target folder
            print(f"Warning: All files in folder {newfolder:s} will be deleted")
            # first delete any existing files
            for oldfile in os.listdir(newfolder):
                f = os.path.join(newfolder, oldfile)
                try:
                    if os.path.islink(f):
                        os.unlink(f)
                    elif os.path.isfile(f):
                        os.remove(f)
                    elif os.path.isdir(f):
                        shutil.rmtree(f)
                except Exception as e:
                    print(f"error: Fail to delete {oldfile:s}. Reason: {e:s}")
                    return 1
        else:
            # folder not exist, create new folder
            print(f"... Creating a new sub-folder {newfolder:s} to store run results")
            try:
                os.mkdir(newfolder)
            except OSError:
                print(f"Error: Fail to create new folder {newfolder:s}")
                return 1
        return 0

    def run_test(root_dir, source_dir, result_dir, test_file: str) -> int:
        """
        Run the given PyChemkin test.

        Parameters
        ----------
            root_dir: str
                rootdir of pytest
            source_dir: str
                the directory containing the PyChemkin source files
            result_dir: str
                the directory containing the test result files
            test_file: str
                name of the PyChemkin test

        Returns
        -------
            ReturnCode: int
                return code from the subprocess run

        """
        #
        # find the working directory
        current_dir = os.getcwd()
        # set sources folder
        source_folder = os.path.join(root_dir, source_dir)
        # check if the source folder exists
        status = PyCKtools.check_folder(source_folder)
        assert 0 == status, "bad source folder name."
        # verify the output folder
        output_folder = os.path.join(current_dir, "outputs")
        # check if the output folder exists
        if PyCKtools.FIRST_PASS == 0:
            status = PyCKtools.create_folder(output_folder)
            assert 0 == status, "fail to get fresh output folder"
        # verify the temporary working folder
        new_working = os.path.join(current_dir, result_dir)
        if PyCKtools.TARGET_FOLDER is None:
            PyCKtools.TARGET_FOLDER = current_dir
        # check if the temporary working folder exists
        if PyCKtools.FIRST_PASS == 0:
            status = PyCKtools.create_folder(new_working)
            assert 0 == status, "fail to get fresh working folder"
            #
            PyCKtools.FIRST_PASS += 1
        # set the test source file
        frun = test_file + ".py"
        frun = os.path.join(source_folder, frun)
        # set the test output file
        f = test_file + ".out"
        output_folder = os.path.join(current_dir, "outputs")
        f = os.path.join(output_folder, f)
        fout = open(f, "w")
        # run test
        # change working directory
        os.chdir(new_working)
        try:
            results = subprocess.run(["python", frun], stdout=fout, check=True)
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
        for out_files in glob.glob(os.path.join(new_working, "*.out")):
            os.remove(out_files)
        for asc_files in glob.glob(os.path.join(new_working, "*.asc")):
            os.remove(asc_files)
        for inp_files in glob.glob(os.path.join(new_working, "*.inp")):
            os.remove(inp_files)
        # return the return code from the subprocess run
        # change working directory
        os.chdir(current_dir)
        #
        return results.returncode

    def load_results(PyCK_result_file) -> dict:
        """
        Read PyChemkin result data from the text file.

        Returns
        -------
            PyCK_result: dict
                PyChemkin test result stored in a dictionary
        """
        with open(PyCK_result_file) as f:
            data = f.read()
        PyCK_result = ast.literal_eval(data)
        return PyCK_result

    def check_list_size(list1: list, list2: list, expected_diff: int = 0) -> bool:
        """
        Verify that the two lists are of the same size.
        """
        return abs(len(list1) - len(list2)) == expected_diff

    def get_file_names(folder_path):
        """
        Get a list of all file names in the specified folder.
        """
        file_names = []
        for entry in os.scandir(folder_path):
            if entry.is_file():
                file_names.append(entry.name)
        return file_names

    def get_tolerances(
        tolerance_name: str,
        state_tol: list[float],
        species_tol: list[float],
        rate_tol: list[float],
    ):
        """
        Find the tolerance values according to the variable type.
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

    def compare_list(
        r_list, b_list, atol: float, rtol: float
    ) -> tuple[int, list, list]:
        """
        Compare the values in the two int or float lists of the same size.

        Returns
        -------
            status: int
                comparison outcome
            bad_ID: list[int]
                list of variable index that the difference fails to satisfy the tolerances
            diff: list[int or float]
                the difference between the values from the two lists
        """
        l_size = len(r_list)
        bad_ID = []
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
                    bad_ID.append(i)
                    if np.isclose(b_value, 0.0, atol=1.0e-9):
                        diff.append(abs(a_diff))
                    else:
                        diff.append(abs(a_diff / b_value))
        # check differences
        iErr = len(bad_ID)
        return iErr, bad_ID, diff

    def init_test_status():
        """
        Initialize the test status values before running the tests.
        """
        PyCKtools.TARGET_FOLDER = None
        PyCKtools.FIRST_PASS = 0

"""
Test script to run all PyChemkin examples under the examples folder
"""
import os
import shutil
import subprocess
import sys

#


class RunPyChemkinTests:
    # defaults
    comparebenchmark = False

    #
    def __init__(self, thisroot):
        # default directories
        self.rootdir = thisroot
        self.sourcedir = os.path.join(thisroot, "examples")
        self.targetdir = os.path.join(thisroot, "results")
        self.benchmarkdir = None
        self.testfolder = os.path.join(thisroot, "results")

    def instructions(self, infile):
        # check test instruction file
        if os.path.isfile(infile):
            # open the instructions from file
            try:
                instruct = open(infile, "r")
            except FileNotFoundError:
                print("File not found!")
            except PermissionError:
                print("You don't have permission to access this file.")
            except OSError as e:
                print("An OS error occurred:", e)
            # read instructions from file
            lines = instruct.readlines()
            n = len(lines)
            count = 0
            for line in lines:
                count += 1
                if count == 1:
                    # first line: full path to the source folder where the examples/tests *.py files are located
                    self.sourcedir = line
                    self.sourcedir = self.sourcedir.rstrip("\n")
                    print(f"... The source files are in folder: {self.sourcedir}")
                    if not os.path.exists(self.sourcedir):
                        print(
                            f"Error: The source folder {self.sourcedir:s} does not exist"
                        )
                        exit()
                elif count == 2:
                    # second line: subdirectory name where all the new results will be stored
                    self.targetdir = line
                    self.targetdir = self.targetdir.rstrip("\n")
                    print(
                        f"... The results will be stored in sub-folder: {self.targetdir}"
                    )
                    # full path of the target folder
                    f = os.path.join(self.rootdir, self.targetdir)
                    error = createfolder(f)
                    if error == 0:
                        self.targetdir = f
                        self.testfolder = f
                    else:
                        exit()
                elif count == 3:
                    # baseline location
                    f = line
                    f = f.rstrip("\n")
                    print(f"... The benchmarks are in folder: {f}")
                    if not os.path.exists(f):
                        self.comparebenchmark = False
                        print(f"Warning: baseline folder {f:s} does not exist")
                        print("         baseline comparisons turned off")
                    else:
                        self.benchmarkdir = f
                        self.comparebenchmark = True
                else:
                    break

            if count > 3:
                print(f"... additional {n-3} line(s) ignored")
        else:
            # cannot find the instruction file
            print(f"Error: Test instruction file {infile} does not exist")
            exit()

    def runtest(self, pyfile):
        # navigate from the working directory to the test directory
        os.chdir(self.testfolder)
        print(f"... current working directory: {os.getcwd()}")
        # set the test source file
        frun = pyfile + ".py"
        frun = os.path.join(self.sourcedir, frun)
        # set the test output file
        f = pyfile + ".out"
        f = os.path.join(self.testfolder, f)
        fout = open(f, "w+")
        # run test
        status = 0
        try:
            subprocess.run(["python", frun], stdout=fout, check=True)
        except subprocess.CalledProcessError as e:
            print("Command returned non-zero exit status:", e.returncode)
        except subprocess.TimeoutExpired as e:
            print("Command timed out:", e.timeout)
        except FileNotFoundError as e:
            print("Command not found:", e)
        except OSError as e:
            print("OS error occurred:", e)
        # close the solution output file
        fout.close()
        return status

    #


def createfolder(newfolder):
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


# main
if __name__ == "__main__":
    # get test instruction file name from the command line argument
    # check number of arguments
    n = len(sys.argv)
    if n != 2:
        print("Error: Requires 1 command line argument: the test instruction file name")
        print("       example: >>> runtests test.txt")
        print("       the instruction file contains the following information:")
        print(
            "       line #1:<full path to the source directory> where the test py files are located."
        )
        print(
            "       line #2:<the target directory> name of sub-directory where all results will be stored."
        )
        print(
            "       line #3:<full path to the baseline directory> where baselines are stored. (optional)"
        )
        exit()
    # get the test instruction file name
    thisroot = os.getcwd()
    print(f"... Working directory: {thisroot}")
    # new test object
    MyTest = RunPyChemkinTests(thisroot)
    infile = os.path.join(thisroot, sys.argv[1])
    MyTest.instructions(infile)
    # print(f"source folder {MyTest.sourcedir}")
    # print(f"target folder {MyTest.targetdir}")
    # go through all *.py files in the source folder
    fcount = 0
    failed = 0
    failedtests = []
    # go through the examples one by one
    for pyfile in os.listdir(MyTest.sourcedir):
        fullpyfile = os.path.join(MyTest.sourcedir, pyfile)
        if os.path.isfile(fullpyfile):
            fname = os.path.basename(pyfile).split(".")[0]
            fext = os.path.basename(pyfile).split(".")[1]
            if fext.lower() == "py":
                # run py files only
                fcount += 1
                print(f"test {fcount:d}: {pyfile}")
                # create a test folder under the target folder
                MyTest.testfolder = os.path.join(MyTest.targetdir, fname)
                error = createfolder(MyTest.testfolder)
                if error != 0:
                    print(f"Error: creating a test folder for test {pyfile}")
                    exit()
                #
                error = 0
                # run the test
                error = MyTest.runtest(fname)
                if error != 0:
                    print(f"error code = {error}")
                if error != 0:
                    failed += 1
                    failedtests.append(pyfile)
            else:
                # skip non python files
                pass

    print(f"... finish running tests in {MyTest.sourcedir}")
    print(f"... completed tests: {fcount}")
    if failed > 0:
        print(f"... number of failed test(s): {failed}")
        for f in failedtests:
            print(f)

    # get back to the original working directory
    os.chdir(thisroot)
    # compare results (in the output files) with the corresponding baselines
    if MyTest.comparebenchmark:
        pass

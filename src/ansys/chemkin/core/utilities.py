# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""pychemkin utilities."""

import os
from pathlib import Path
import re
import secrets
import shutil
from typing import NoReturn, Union

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.logger import logger

ck_rng = None  # random number generator object


def log_critical_error(msg_parts: list[str]) -> None:
    """Log a critical error message."""
    logger.critical(Color.SPACE.join(msg_parts))


def log_critical_message(msg_parts: list[str]) -> None:
    """Log a critical message without additional formatting."""
    log_critical_error(msg_parts)


def log_error_message(msg_parts: list[str]) -> None:
    """Log an error message without terminating execution."""
    logger.error(Color.SPACE.join(msg_parts))


def log_info_message(msg_parts: list[str]) -> None:
    """Log an informational message without additional formatting."""
    logger.info(Color.SPACE.join(msg_parts))


def error_and_exit(msg_parts: list[str]) -> NoReturn:
    """Log an error message and terminate execution."""
    log_error_message(msg_parts)
    raise SystemExit()


def critical_and_exit(msg_parts: list[str]):
    """Log a critical message and terminate execution."""
    log_critical_message(msg_parts)
    raise SystemExit()


def _log_warning_message(msg_parts: list[str]) -> None:
    """Log a warning message with consistent message joining."""
    logger.warning(Color.SPACE.join(msg_parts))


def where_element_in_array_1d(
    arr: Union[npt.NDArray[np.double], npt.NDArray[np.int32]], target: float
) -> tuple[int, npt.NDArray[np.int32]]:
    """Find number of occurrence and indices of the target value in an array."""
    """Find the number of occurrence and the element index in
    the 1D arr array that matches the target value. Using numpy.argwhere
    might be more efficient. However, the numpy method returns a list of
    lists of occurrence indices while this might be necessary
    for general applications, it is an overkill for simple 1D array cases.

    Parameters
    ----------
        arr: integer or double array
            the reference 1D integer or double array
        target: integer or double scalar
            the target value to be matched

    Returns
    -------
        number_of_occurrences: integer
        occurrence_index: integer array

    """
    count = 0
    # check arr array size
    arr_size = len(arr)
    if arr_size == 0:
        # nothing in arr
        return count, np.zeros(0, dtype=np.int32)
    temp_index = np.zeros(arr_size, dtype=np.int32)
    value = type(arr[0])(target)
    # find all the matching occurrences
    for m in range(arr_size):
        if arr[m] == value:
            temp_index[count] = m
            count += 1
    if count == 0:
        # target is not in arr
        where_index = np.zeros(0, dtype=np.int32)
    else:
        where_index = temp_index[:count]
    return count, where_index


def bisect(ileft: int, iright: int, x: float, xarray) -> int:
    """Use bisectional method to find the largest index in the array."""
    """Use bisectional method to find the largest index in the xarray
    of which its value is small or equal to the target x value.

    Parameters
    ----------
        ileft: integer
            index of xarray that represents the current lower bound
        iright: integer
            index of xarray that represents the current upper bound
        x: double
            target value
        xarray: double array
            a sorted array containing all x values in strictly ascending
            order x[i] < x[i+1]

    Returns
    -------
        itarget: integer
            the largest index in the xarray of which its value is small
            or equal to the target x value

    """
    if (iright - ileft) > 1:
        ihalf = int((ileft + iright) / 2)
        if xarray[ihalf] > x:
            iright = ihalf
        else:
            ileft = ihalf
        itarget = bisect(ileft, iright, x, xarray)
        # print(f"lower bound = {ileft}, upper bound = {iright}, target = {itarget}")
    else:
        itarget = ileft
    return itarget


def find_interpolate_parameters(
    x: float, xarray: npt.NDArray[np.double]
) -> tuple[int, float]:
    """Find the indices branket the given value."""
    """Find indices ileft and iright such that
       xarray[ileft] <= x <= xarray[iright] where iright = ileft + 1

    Parameters
    ----------
        x: double
            target value
        xarray: double array
            a sorted array containing all x values in strictly
            ascending order x[i] < x[i+1]

    Returns
    -------
    itarget: integer
        the largest index in the xarray of which its value is small or
        equal to the target x value
    ratio: double
        the distance ratio = (x - xarray[ileft])/(xarray[ileft+1] - xarray[ileft])

    """
    iarraysize = len(xarray)
    if x == xarray[0]:
        # x = xarray[0]
        itarget = 0
        ratio = 0.0e0
        return itarget, ratio
    if x == xarray[iarraysize - 1]:
        # x = xarray[max]
        itarget = iarraysize - 2
        ratio = 1.0e0
        return itarget, ratio
    if (x - xarray[0]) * (x - xarray[iarraysize - 1]) > 0.0e0:
        # x value is out of bound
        msg = [
            Color.PURPLE,
            "the target value x=",
            str(x),
            "does not fall between",
            str(xarray[0]),
            "and",
            str(xarray[iarraysize - 1]),
            Color.END,
        ]
        log_error_message(msg)
    # bisect method
    ileft = 0
    iright = iarraysize - 1
    itarget = bisect(ileft, iright, x, xarray)
    ratio = (x - xarray[itarget]) / (xarray[itarget + 1] - xarray[itarget])
    return itarget, ratio


def interpolate_array(
    x: float, x_array: npt.NDArray[np.double], y_array: npt.NDArray[np.double]
) -> float:
    """Find the y-value corresponding to the x-value in data pairs."""
    """Find the value in the y_array from the interpolation parameters ileft and ratio
        y = (1-ratio)* y_array[ileft] + ratio * y_array[ileft+1]
        where ileft and ratio are determined from the target x value and the xarray

    Parameters
    ----------
        x: double
            target value
        x_array: double array
            a sorted array containing all x values in strictly
            ascending order x[i] < x[i+1]
        y_array: double array
            dependent variable array

    Returns
    -------
        y: double
            the interpolated dependent variable value corresponding the given x value

    """
    # find the interpolation parameters
    ileft, ratio = find_interpolate_parameters(x, x_array)
    # perform the interpolation to find the y value from the yarray
    y = (1.0e0 - ratio) * y_array[ileft]
    y += ratio * y_array[ileft + 1]
    return y


# stoichiometric
#
def _nonzero_element_in_array_1d(
    arr: Union[npt.NDArray[np.int32], npt.NDArray[np.double]], threshold: float = 0.0
) -> tuple[int, npt.NDArray[np.int32]]:
    """Find the number of occurrence and the indices of the non-zero members."""
    """Find the number of occurrence and the indices of the non-zero (> 0) element
    in the array arr. Using numpy.nonzero might be more efficient. However,
    the numpy method returns a list of lists of occurrence indices while this might
    be necessary for general applications, it is an overkill for simple 1D array cases.

    Parameters
    ----------
        arr: 1-D integer or double array
            the reference array with non-negative integer or double
        threshold: integer or double, optional, default = 0
            the threshold used as the reference value of zero

    Returns
    -------
        nonzero_count: integer
            number_of_occurrences
        nonzero_index: 1-D integer array
            occurrence_index

    """
    # find the number of non-zero counts
    nonzero_count = np.count_nonzero(arr)
    if nonzero_count == 0:
        return int(nonzero_count), np.zeros(0, dtype=np.int32)
    nonzero_index = np.zeros(nonzero_count, dtype=np.int32)
    thrd = type(arr[0])(threshold)
    j = 0
    # find all non-zero occurrences
    for m in range(len(arr)):
        if arr[m] > thrd:
            nonzero_index[j] = m
            j += 1
    return int(nonzero_count), nonzero_index


def random(range: Union[None, tuple[float, float]] = None) -> float:
    """Generate a (reproducible) random floating number."""
    """Generate a (reproducible) random floating number value >= 0.0 and < 1.0
    by using the Numpy pseudo-random number generator.
    If the range tuple (a, b) is given, the random number will
    have a value >= a and < b.

    Parameters
    ----------
        range: tuple of floats (a, b) and b > a, default = (0.0, 1.0)
            the range of the random number values

    Returns
    -------
        random: float
            random number

    """
    global ck_rng
    if ck_rng is None:
        # need initialization
        # get the seeding value
        seed = secrets.randbits(128)
        seed -= 54231
        # create a random number generator instance
        ck_rng = np.random.default_rng(seed)

    if range is None:
        # return value [0, 1)
        return ck_rng.random()
    else:
        # return value [a, b)
        width = range[1] - range[0]
        return range[0] + ck_rng.random() * width


def random_pick_integers(
    numb_picks: int, source_integers: list[int]
) -> tuple[list[int], list[int]]:
    """Return a list of random integers from a given integer list."""
    """
    Return a list of unique integers randomly selected from the source integer list.
    The total number of the random integers to be picked, numb_picks, must be
    < the size of the source integer list.

    Parameters
    ----------
        numb_picks: integer
            the total number of random integers to be picked
        source_integers: list of integers
            the source list of the integers to be picked randomly

    Returns
    -------
        picked_list: list of integer, dimension = numb_picks
            list of unique random integers picked from the source list
        unpicked_list: list of integer, dimension = len(source_integers) - numb_picks
            list of integers of the source list that are not picked
    """
    # check
    if numb_picks < 1:
        msg = [
            Color.PURPLE,
            "the total number of random picks must > 1.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()

    max_picks = len(source_integers)
    if max_picks < numb_picks:
        msg = [
            Color.PURPLE,
            "The number of random picks requested",
            str(numb_picks),
            "is more than the number of integers",
            "in the source list",
            str(max_picks),
            ".",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    #
    picked_list: list[int] = []
    # set source pointer
    arrow = max_picks
    # random pick
    for i in range(numb_picks):
        # randomly pick a slot in source_list
        slot = int(random() * float(arrow)) + 1
        # add the value of the slot to the picked list
        picked_list.append(source_integers[slot - 1])
        if slot != arrow:
            # swap the values of the picked slot and the last available slot
            id1 = slot - 1
            id2 = arrow - 1
            replaced = source_integers[id1]
            source_integers[id1] = source_integers[id2]
            source_integers[id2] = replaced
        # reduce the available slot by 1
        arrow -= 1
    # save the unpicked integers in the source list
    unpicked_list = source_integers[0:arrow]
    return picked_list, unpicked_list


def find_file(filepath: str, partialfilename: str, fileext: str) -> str:
    """Find the correct version of the given partial file name."""
    """This is mostly to handle the different years/versions of the
    MFL mechanisms that come with the Ansys Chemkin installation.

    Parameters
    ----------
        filepath: string
            the directory where the file is located
        partialfilename: string
            the leading portion of the file name
        fileext: string
            file extension

    Returns
    -------
        thefile: string
            full path name of the file, = ""
            if no file matches the 'partialname' in the 'filepath'

    """
    thefile = ""
    for file in Path(filepath).iterdir():
        if ("." + fileext) == file.suffix:
            if re.search(partialfilename, file.name):
                thefile = str(file.resolve())
                break
    return thefile


def copy_file(
    sourcepath: str,
    targetpath: str,
    filename: str,
    new_file_name: Union[str, None] = None,
):
    """Copy a file from the source folder to the target folder."""
    """
    Copy a file from the source folder to the target folder.
    The copied file will have the same name as the source file unless
    the new file name is given.

    Parameters
    ----------
        sourcepath: string
            the directory where the file to be copied is located
        targetpath: string
            the destinated folder where the new file will be sent
        filename: string
            the source file name (with file extension)
        new_file_name: string, optional.
            the new of te newly copied file if it is different
            from the original file name
    """
    # Check the destination directory
    if not Path(sourcepath).is_dir():
        msg = [
            Color.PURPLE,
            "The source folder",
            sourcepath,
            "does not exist.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # Check the target directory
    if not Path(targetpath).is_dir():
        msg = [
            Color.PURPLE,
            "The target folder",
            targetpath,
            "does not exist.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    #
    source_item = Path(sourcepath) / filename
    if new_file_name is None:
        newfile = filename
    else:
        newfile = new_file_name
    target_item = Path(targetpath) / newfile

    try:
        shutil.copy(source_item, target_item)
        print(f"Copied '{filename}' to '{targetpath}'")
    except Exception as e:
        msg = [
            Color.PURPLE,
            "failed to copy",
            filename,
            "to folder",
            targetpath,
            ",\n",
            Color.SPACEx6,
            "error message:",
            str(e),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)


def delete_files_by_extension(targetpath: str, extensions: list[str]):
    """Delete files with specified extensions."""
    """
    Delete files with specified extensions in a given directory if they exist.

    Parameters
    ----------
        targetpath: string
            The path to the directory to scan.
        extensions: list of strings
            A list of file extensions (e.g., ['.txt', '.log']) to delete.
    """
    # get all files with given extensions in the target folder
    files: list[Path] = []
    for e in extensions:
        ext = "*" + e
        for f in Path(targetpath).glob(ext):
            files.append(f)

    for f in files:
        # Check if the file's extension is in the list of extensions to delete
        if isinstance(f, Path):
            if f.is_file():
                f.unlink()
                print(f"Deleted: {str(f)}")
            else:
                print(f"Error deleting {str(f)} : is not a file.")


def check_jupyter_notebook() -> bool:
    """Check if Pychemkin is running in any Jupyter environment."""
    """
    Check if Pychemkin is running in any Jupyter environment.

    Returns
    -------
        status: boolean
            True: running in Jupyter environment
            False: not in Jupyter environment
    """
    try:
        import IPython  # pyright: ignore[reportMissingModuleSource]

        ip = IPython.get_ipython()
        if ip is not None and "IPKernelApp" in ip.config:
            return True
    except (ImportError, AttributeError):
        pass
    return False


def interpolate_point(
    x_value: float,
    x_array: Union[list[float], npt.NDArray[np.double]],
    y_array: Union[list[float], npt.NDArray[np.double]],
) -> tuple[int, float]:
    """Interpolate value from a binary data array."""
    """
    Get the y value from the given x value based on the given
    correlation (y = f(x)) defined by the x_array and the y_array.

    Parameters
    ----------
        x_value: double
            the x value to be used to obtain the y value by interpolation
        x_array: 1-D double array
            the data for all x values
        y_array: 1-D double array
            the y values for each x value data point in the x_array

    Returns
    -------
        index_jp: integer
            the x_array index where x(index_jp - 1) < x_value <= x(index_jp)
        y_value: double
            the interpolated y value: y_value = f(x_value)
    """
    # check array lengths
    if len(x_array) != len(y_array):
        msg = [
            Color.PURPLE,
            "The given arrays must have the same length.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check x range:
    x0 = x_array[0]
    xlast = x_array[len(x_array) - 1]
    if (x_value - x0) * (x_value - xlast) > 0.0:
        # x value is not covered by the x_array data points
        msg = [
            Color.PURPLE,
            "The given x value",
            str(x_value),
            "is out of the range defined in the x array [",
            str(x0),
            ",",
            str(xlast),
            "]",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # find the x range and verify x_array is monotonic
    # ascending order: up = True
    change = xlast - x0
    not_good = bool(np.isclose(np.abs(change), 0.0, atol=1.0e-8))
    if not not_good:
        up = change > 0.0
        x_old = xlast
        #
        for j, x in enumerate(x_array):
            if j > 0:
                change = x - x_old
                if change < 0.0:
                    # x value decreases for ascending array
                    not_good = bool(up)
                elif change > 0.0:
                    # x value increases for descending array
                    not_good = bool(not up)
                else:
                    # flat
                    not_good = True
                    break
            x_old = x
    #
    if not_good:
        msg = [
            Color.PURPLE,
            "The given x array must be monotonic.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # interpolation
    for j, x in enumerate(x_array):
        if j == 0:
            change = x_array[j] - x_value
            if np.isclose(np.abs(change), 0.0, atol=1.0e-8):
                return j, y_array[j]
        else:
            dx = x_array[j] - x_array[j - 1]
            dxjp = x_array[j] - x_value
            frac = dxjp / dx
            if frac >= 0.0:
                y_value = (1.0e0 - frac) * y_array[j] + frac * y_array[j - 1]
                return j, y_value
    return -1, 0.0


class WorkingFolders:
    """folder/file utilities for multiple runs."""

    def __init__(self, dir_name: str, root_dir: str):
        """Create and change to the designated working directory."""
        """
        Create and change to the designated working directory to
        run the parameter study cases or to run the optimization cases.

        Parameters
        ----------
            dir_name: string
                name of the working folder
            root_dir: string
                name of the top/root folder of the working folders
        """
        # create or clean up the working directory for current simulation
        self.root_dir = Path(root_dir)
        # check the root folder
        if not self.root_dir.is_dir():
            msg = [
                Color.PURPLE,
                "The 'root' folder",
                root_dir,
                "does not exist.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self.work_dir = self.root_dir / dir_name
        # check the working folders under the root folder
        # for multiple simulation runs
        if self.work_dir.is_dir():
            # directory exists
            for f in self.work_dir.iterdir():
                if f.is_file():
                    f.unlink()
        else:
            # create a new directory
            self.work_dir.mkdir(parents=True, exist_ok=True)
        # change to the working directory
        os.chdir(self.work_dir)

    @property
    def root(self) -> str:
        """Get the "root" directory of the current structure."""
        """
        Get the "root" directory of the current structure.

        Returns
        -------
            root: string
                the upper level folder
        """
        return str(self.root_dir)

    @property
    def work(self) -> str:
        """Get the "work" directory of the current structure."""
        """
        Get the "work" directory of the current structure.

        Returns
        -------
            work: string
                the working level folder
        """
        return str(self.work_dir)

    def done(self):
        """Change back to the root directory."""
        """
        Change back to the root directory.
        """
        # change to the working directory for this run
        os.chdir(self.root_dir)

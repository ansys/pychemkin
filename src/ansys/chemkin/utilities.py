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

"""
pychemkin utilities
"""
from typing import Union

from ansys.chemkin.chemistry import Chemistry
from ansys.chemkin.color import Color
from ansys.chemkin.logger import logger
import numpy as np
import numpy.typing as npt


def where_element_in_array_1D(
    arr: Union[npt.NDArray[np.double], npt.NDArray[np.int32]], target: Union[int, float]
) -> tuple[int, npt.NDArray[np.int32]]:
    """
    Find the number of occurrence and the element index in the 1D arr array that matches the target value.
    Using numpy.argwhere might be more efficient. However, the numpy method returns a list of lists of occurrence indices
    while this might be necessary for general applications, it is an overkill for simple 1D array cases.

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
        return count, []
    temp_index = np.zeros(arr_size, dtype=np.int32)
    value = type(arr[0])(target)
    # find all the matching occurrences
    for m in range(arr_size):
        if arr[m] == value:
            temp_index[count] = m
            count += 1
    if count == 0:
        # target is not in arr
        where_index = []
    else:
        where_index = temp_index[:count]
    return count, where_index


def bisect(ileft: int, iright: int, x: float, xarray) -> int:
    """
    Use bisectional method to find the largest index in the xarray of which its value is small or equal to the target x value

    Parameters
    ----------
        ileft: integer
            index of xarray that represents the current lower bound
        iright: integer
            index of xarray that represents the current upper bound
        x: double
            target value
        xarray: double array
            a sorted array containing all x values in strictly ascending order x[i] < x[i+1]

    Returns
    -------
        itarget: integer
            the largest index in the xarray of which its value is small or equal to the target x value
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
    """
    Find the index ileft that
       xarray[ileft] <= x <= xarray[iright] where iright = ileft + 1

    Parameters
    ----------
        x: double
            target value
        xarray: double array
            a sorted array containing all x values in strictly ascending order x[i] < x[i+1]

    Returns
    -------
    itarget: integer
        the largest index in the xarray of which its value is small or equal to the target x value
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
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # bisect method
    ileft = 0
    iright = iarraysize - 1
    itarget = bisect(ileft, iright, x, xarray)
    ratio = (x - xarray[itarget]) / (xarray[itarget + 1] - xarray[itarget])
    return itarget, ratio


def interpolate_array(
    x: float, x_array: npt.NDArray[np.double], y_array: npt.NDArray[np.double]
) -> float:
    """
    Find the value in the y_array from the interpolation parameters ileft and ratio
        y = (1-ratio)* y_array[ileft] + ratio * y_array[ileft+1]
        where ileft and ratio are determined from the target x value and the xarray

    Parameters
    ----------
        x: double
            target value
        x_array: double array
            a sorted array containing all x values in strictly ascending order x[i] < x[i+1]
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


def create_mixture_recipe_from_fractions(
    chemistry_set: Chemistry, frac: npt.NDArray[np.double]
) -> tuple[int, list[tuple[str, float]]]:
    """
    Build a PyChemkin mixture recipe/formula from a species fraction array (i.e., mixture mole/mass composition).
    This mixture recipe can then be used to create the corresponding Mixture object.

    Parameters
    ----------
    chemistry_set: Chemistry object
        the Chemistry object will be used to create the mixture
    frac: double array
        mole or mass fractions of the mixture

    Returns
    -------
        count: integer
            the size of the recipe list containing [gas species, mole/mass fraction] tuples
        recipe: list of tuples, [(species_symbol, fraction), ... ]
            non-zero mixture composition corresponding to the given mole/mass fraction array
    """
    # initialization
    count = 0
    recipe = []
    # check Chemistry object
    if not isinstance(chemistry_set, Chemistry):
        msg = [
            Color.PURPLE,
            "the first argument must be a Chemistry object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check array size
    numb_species = chemistry_set.KK
    if len(frac) != numb_species:
        msg = [
            Color.PURPLE,
            "the size of the fractions array does not match",
            "the number of species in the chemistry set.\n",
            "the fraction array size should be",
            str(numb_species),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # build the recipe from frac array
    for k in range(numb_species):
        if frac[k] > 0.0e0:
            species_symbol = chemistry_set.KSymbol[k]
            recipe.append((species_symbol, frac[k]))
            count += 1
    return count, recipe


# stoichiometric
#
def _nonzero_element_in_array_1D(
    arr: Union[npt.NDArray[np.int32], npt.NDArray[np.double]], threshold: float = 0.0
) -> tuple[int, npt.NDArray[np.int32]]:
    """
    Find the number of occurrence and the indices of the non-zero (> 0) element in the array arr.
    Using numpy.nonzero might be more efficient. However, the numpy method returns a list of lists of occurrence indices
    while this might be necessary for general applications, it is an overkill for simple 1D array cases.

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
        return nonzero_count, []
    nonzero_index = np.zeros(nonzero_count, dtype=np.int32)
    thrd = type(arr[0])(threshold)
    j = 0
    # find all non-zero occurrences
    for m in range(len(arr)):
        if arr[m] > thrd:
            nonzero_index[j] = m
            j += 1
    return nonzero_count, nonzero_index


def calculate_stoichiometrics(
    chemistryset: Chemistry,
    fuel_molefrac: npt.NDArray[np.double],
    oxid_molefrac: npt.NDArray[np.double],
    prod_index: npt.NDArray[np.int32],
) -> tuple[float, npt.NDArray[np.double]]:
    """
    calculate the stoichiometric coefficients of the complete combustion reaction of the given fuel and oxidizer mixtures.
    Consider the complete combustion of the fuel + oxidizer mixture
    ::
        (fuel species) + alpha*(oxidizer species) <=> nu(1)*prod(1) + ... + nu(numb_prod)*prod(numb_prod)

    The number of unknowns is equal to the number of elements that make of all the fuel and oxidizer species. And the
    number of product species must be one less than the number of unknowns.
    The unknowns
    ::
        alpha is the stoichiometric coefficient multiplier of the oxidizer species
        nu(1), ... nu(numb_prod) are the stoichiometric coefficients of the product species

    The conservation of elements yields a set of linear algebraic equations
    ::
        A x = b
    in which x = [ -alpha | nu(1), ...., nu(numb_prod) ]  (a vector of size numb_elem ) can be obtained.

    Parameters
    ----------
        chemistryset: Chemistry object
            the Chemistry object used to create the fuel and the oxidizer mixtures
        fuel_molefrac: 1-D double array
            mole fractions of the fuel mixture
        oxid_molefrac: 1-D double array
            mole fractions of the oxidizer mixture
        prod_index: 1-D integer array
            the species indices of the complete combustion products

    Returns
    -------
        alpha: double
            oxidizer_coefficient_multiplier
        nu: 1-D double array
            stoichiometric_coefficients_of_products
    """
    # check the Chemistry object
    if not isinstance(chemistryset, Chemistry):
        msg = [
            Color.PURPLE,
            "the first argument must be a Chemistry object.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # get the number of elements and the number of gas species from the chemistry set
    numb_elem = chemistryset.MM
    numb_species = chemistryset.KK
    # find fuel species array size
    kfuel = len(fuel_molefrac)
    # find oxidizer array size
    koxid = len(oxid_molefrac)
    # check fuel composition array
    if numb_species != kfuel:
        msg = [
            Color.PURPLE,
            "the fuel species array size must be",
            str(numb_species),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # check oxidizer composition array
    if numb_species != koxid:
        msg = [
            Color.PURPLE,
            "the oxidizer species array size must be",
            str(numb_species),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    # find number of product species
    numb_prod = len(prod_index)
    # find fuel species index and count
    numb_fuel, fuel_index = _nonzero_element_in_array_1D(fuel_molefrac)
    # find oxidizer species index and count
    numb_oxid, oxid_index = _nonzero_element_in_array_1D(oxid_molefrac)
    # the same species cannot be fuel and oxidizer at the same time
    for i in oxid_index:
        j, j_index = where_element_in_array_1D(fuel_index, i)
        if j != 0:
            msg = [
                Color.YELLOW,
                "species",
                chemistryset.KSymbol[i],
                "is in both the fuel and the oxidizer mixtures.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

    # find the actual number of elements in fuel and oxidizer
    elem_tally = np.zeros(numb_elem, dtype=np.int32)
    # elements in the fuel species
    for k in fuel_index:
        for m in range(numb_elem):
            elem_count = chemistryset.SpeciesComposition(m, k)
            if elem_count > 0:
                elem_tally[m] += elem_count
    # elements in the oxidizer species
    for k in oxid_index:
        for m in range(numb_elem):
            elem_count = chemistryset.SpeciesComposition(m, k)
            if elem_count > 0:
                elem_tally[m] += elem_count
    numb_coreelem, coreelem_index = _nonzero_element_in_array_1D(elem_tally)
    # check the number of product species
    if numb_prod != (numb_coreelem - 1):
        msg = [
            Color.PURPLE,
            "the number of product species must be",
            str(numb_coreelem - 1),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.error(this_msg)
        exit()
    else:
        # check product elements
        # find elements in product species
        elem_prod = np.zeros(numb_elem, dtype=np.int32)
        for k in prod_index:
            for m in range(numb_elem):
                elem_count = chemistryset.SpeciesComposition(m, k)
                if elem_count > 0:
                    elem_prod[m] += elem_count
        numb_prodelem, prodelem_index = _nonzero_element_in_array_1D(elem_prod)
        # check elements in the products and in the fuel and oxidzer mixtures
        elname = ""
        if numb_prodelem == numb_coreelem:
            for m in prodelem_index:
                if m in coreelem_index:
                    pass
                else:
                    elname = chemistryset.element_symbols[m]
                    msg = [
                        Color.PURPLE,
                        "element",
                        elname,
                        "in products is not in fuel or oxidizer mixtures.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
        else:
            msg = [
                Color.PURPLE,
                "the number of product elements must be the same",
                "as the number of elements in fuel and oxidizer\n",
                Color.SPACEx6,
                "the number of elements in products:",
                str(numb_prodelem),
                "\n",
                Color.SPACEx6,
                "the number of elements in the fuel and the oxidizer:",
                str(numb_coreelem),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
    # create arrays of the linear algebraic system
    A = np.zeros((numb_coreelem, numb_coreelem), dtype=np.double)
    b = np.zeros(numb_coreelem, dtype=np.double)
    # construct the (numb_coreelem x 1) b array on the right-hand side
    #   b = [SUM_k(NCF(1,k)*fuel_molefrac(k)), ... SUM_k(NCF(numb_elem,k)*fuel_molefrac(k))]
    for m in range(numb_coreelem):
        b[m] = 0.0e0
        this_elem = coreelem_index[m]
        for k in range(numb_species):
            elem_count = chemistryset.SpeciesComposition(this_elem, k)
            b[m] += elem_count.astype(np.double) * fuel_molefrac[k]
            # first column of A[1:numb_coreelem, 1]
            A[m][0] += elem_count.astype(np.double) * oxid_molefrac[k]
    # construct the sub-matrix on the right of A[1:numb_coreelem, 2:numb_prod]
    for m in range(numb_coreelem):
        this_elem = coreelem_index[m]
        for k in range(numb_prod):
            k_prod = prod_index[k]
            A[m][k + 1] = chemistryset.SpeciesComposition(this_elem, k_prod)
    # solve the linear system: A x = b
    x = np.linalg.solve(A, b)
    alpha = -x[0]
    nu = x[1:numb_coreelem]
    return alpha, nu

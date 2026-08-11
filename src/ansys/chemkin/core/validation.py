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

"""Reusable validation helpers for core modules."""

from pathlib import Path

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.utilities import log_error_message


def validate_chemid(chemid: int) -> None:
    """Validate Chemistry set index."""
    if chemid < 0:
        log_error_message([Color.PURPLE, "invalid chemistry.", Color.END])
        exit()


def validate_temperature(temp: float) -> None:
    """Validate gas temperature."""
    if temp <= 10.0:
        log_error_message([Color.PURPLE, "invalid temperature value.", Color.END])
        exit()


def validate_pressure_temperature(pressure: float, temp: float) -> None:
    """Validate pressure and temperature pair."""
    if pressure <= 0.0 or (pressure * temp) <= 0.0:
        log_error_message(
            [
                Color.PURPLE,
                "invalid pressure and/or temperature value(s).",
                Color.END,
            ]
        )
        exit()


def validate_density(density: float) -> None:
    """Validate density value."""
    if density <= 0.0:
        log_error_message([Color.PURPLE, "invalid density value.", Color.END])
        exit()


def validate_fraction_arrays(frac, wt, mode: str) -> int:
    """Validate fraction and molecular-weight arrays are same size."""
    kgas = len(frac)
    if kgas != len(wt):
        log_error_message(
            [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
        )
        exit()
    return kgas


def validate_species_array_size(expected: int, actual: int, context: str) -> None:
    """Validate a species-array size against expected species count."""
    if expected != actual:
        log_error_message(
            [
                Color.PURPLE,
                "the",
                context,
                "species array size must be",
                str(expected),
                Color.END,
            ]
        )
        exit()


def validate_minimum_value(value: float, minimum: float, message: str) -> None:
    """Validate that a scalar value is greater than or equal to a minimum."""
    if value < minimum:
        log_error_message([Color.PURPLE, message, Color.END])
        exit()


def validate_equal_list_lengths(
    first_label: str,
    first_len: int,
    second_label: str,
    second_len: int,
) -> None:
    """Validate that two list lengths are equal."""
    if first_len != second_len:
        log_error_message(
            [
                Color.PURPLE,
                "the",
                first_label,
                "list must have the same,",
                "length as the",
                second_label,
                "list.\n",
                Color.SPACEx6,
                first_label,
                "list length =",
                str(first_len),
                "\n",
                Color.SPACEx6,
                second_label,
                "list length =",
                str(second_len),
                Color.END,
            ]
        )
        exit()


def validate_fraction_array_size(actual: int, expected: int) -> None:
    """Validate fractions-array size against chemistry species count."""
    if actual != expected:
        log_error_message(
            [
                Color.PURPLE,
                "the size of the fractions array does not match",
                "the number of species in the chemistry set.\n",
                "the fraction array size should be",
                str(expected),
                Color.END,
            ]
        )
        exit()


def validate_file_exists(file_path: str, context_label: str = "file") -> None:
    """Validate that a filesystem path points to an existing file."""
    if not Path(file_path).is_file():
        log_error_message(
            [
                Color.PURPLE,
                "the specified",
                context_label,
                "does not exist:",
                file_path,
                Color.END,
            ]
        )
        exit()

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
Constants used by Chemkin utilities and models.
"""

# == Chemkin module global parameters -- DO NOT MODIFY without asking Chemkin development team members
boltzmann = 1.3806504e-16  # Boltzmann constant [ergs/K] (double scalar)
avogadro = 6.02214179e23  # Avogadro number [1/mole] (double scalar)
Patm = 1.01325e06  # atmospheric pressure [dynes/cm2] (double scalar)
ergs_per_joule = 1.0e7  # ergs per joule [ergs/J] (double scalar)
joules_per_calorie = 4.184e0  # joules per calorie [J/cal] (double scalar)
ergs_per_calorie = (
    joules_per_calorie * ergs_per_joule
)  # ergs per calorie [erg/cal] (double scalar)
ergs_per_eV = 1.602176487e-12  # ergs per eV [erg/volt] (double scalar)
eV_per_K = ergs_per_eV / boltzmann  # eV per K [volt/K] (double scalar)
RGas = boltzmann * avogadro  # universal gas constant R [ergs/mol-K] (double scalar)
RGas_Cal = (
    RGas * 1.0e-7 / joules_per_calorie
)  # universal gas constant R [cal/mol-K] (double scalar)
# == end of global constants

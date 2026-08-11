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

"""Surface chemistry controller for reactor models.

This module isolates surface-chemistry related behavior from ``ReactorModel``.
The controller is composition-based and operates on a host object that exposes
reactor state and helper methods.
"""

from __future__ import annotations

from ctypes import c_int
from typing import Any, Protocol

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core.chemistry import Chemistry, verify_chemset_sizes
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.surface import Surface
from ansys.chemkin.core.utilities import (
    error_and_exit,
    log_error_message,
    log_info_message,
)
from ansys.chemkin.core.validation import validate_minimum_value


class SurfaceChemistryHost(Protocol):
    """Protocol for reactor objects that can host surface chemistry behavior."""

    has_surface_chemistry: bool
    numbspecies: int
    numbmaterials: int
    surface_chemistry: Any
    _chemset_index: c_int
    _surfaceratemultiplier: float
    _numb_surf_area_set: int
    _surf_area_set: list[int]
    _material_area_fraction: list[float]

    @property
    def temperature(self) -> float:
        """Return the current reactor temperature."""

    def setkeyword(self, key: str, value: Any):
        """Set a Chemkin keyword on the host reactor."""


class SurfaceChemistryController:
    """Controller that encapsulates reactor surface-chemistry operations."""

    _FOUR_SPACES = "    "

    def __init__(self, reactor: SurfaceChemistryHost):
        """Initialize the controller with a reactor host."""
        self._reactor = reactor

    def activate_surface_chemistry(self, chem: Chemistry, mode: str = "normal"):
        """Activate surface chemistry on reactor wall materials."""
        chemset_index = c_int(chem.chemid)
        if chemset_index.value != self._reactor._chemset_index.value:
            msg = [
                Color.PURPLE,
                "incompatible chemistry sets:\n",
                Color.SPACEx6,
                "chemistry set index of the reactor mixture =",
                str(self._reactor._chemset_index.value),
                "\n",
                Color.SPACEx6,
                "the given chemistry set index =",
                str(chemset_index.value),
                Color.END,
            ]
            error_and_exit(msg)

        if not verify_chemset_sizes(
            chem_set_index=chem.chemid, num_gas_species=self._reactor.numbspecies
        ):
            msg = [
                Color.PURPLE,
                "incompatible chemistry sets:",
                "number of gas species in the chemistry sets are different.",
                Color.END,
            ]
            error_and_exit(msg)

        if not chem.verify_surface_mechanism():
            chem.no_surface_mechanism_declaration()
            exit()

        if self._reactor.surface_chemistry is not None:
            if isinstance(self._reactor.surface_chemistry, Surface):
                msg = [
                    Color.YELLOW,
                    "reactor already has surface chemistry set up",
                    "and will be overridden.",
                    Color.END,
                ]
                log_info_message(msg)
            del self._reactor.surface_chemistry

        self._reactor.surface_chemistry = Surface(chem)
        self._reactor.numbmaterials = self._reactor.surface_chemistry.number_materials

        mat_list = self._reactor.surface_chemistry.get_material_names()
        if mode == "normal":
            print("List of available surface materials:\n")
            for i, name in enumerate(mat_list):
                material = self._reactor.surface_chemistry.get_surface_material(name)
                print(f"*** Material Index = {i}")
                material.information()

        self._reactor._surfaceratemultiplier = 1.0e0
        self._reactor._numb_surf_area_set = 0
        self._reactor._surf_area_set = []
        self._reactor._material_area_fraction = []

    def verify_surface_chemistry(self) -> bool:
        """Verify that the host reactor supports surface chemistry."""
        if not self._reactor.has_surface_chemistry:
            self.no_surface_mechanism_declaration()
        return self._reactor.has_surface_chemistry

    def no_surface_mechanism_declaration(self):
        """Inform users that surface chemistry is not available."""
        msg = [
            Color.PURPLE,
            "reactor does not have surface chemistry.\n",
            Color.SPACEx6,
            "please instantiate the reactor with Mixture/Inlet",
            "created from a Chemistry Set with surface chemistry.",
            Color.END,
        ]
        log_error_message(msg)

    def require_surface_chemistry(self) -> bool:
        """Return availability and log guidance when surface chemistry is absent."""
        if not self._reactor.has_surface_chemistry:
            self.no_surface_mechanism_declaration()
            return False
        return True

    def no_surface_material(self, mat_name: str):
        """Log an error for an unknown surface material."""
        msg = [
            Color.PURPLE,
            mat_name.rstrip(),
            "is not a valid surface material",
            Color.END,
        ]
        log_error_message(msg)

    @property
    def get_numb_material(self) -> int:
        """Return number of materials in the active surface mechanism."""
        return self._reactor.numbmaterials

    def get_material_names(self) -> list[str]:
        """Return all surface material names."""
        return self._reactor.surface_chemistry.get_material_names()

    def get_site_species_names(self) -> list[list[str]]:
        """Return site-species names for all surface materials."""
        site_names: list[list[str]] = []
        for material_name in self.get_material_names():
            site_names.append(
                self._reactor.surface_chemistry.get_site_symbols(mat_name=material_name)
            )
        return site_names

    def get_bulk_species_names(self) -> list[list[str]]:
        """Return bulk-species names for all surface materials."""
        bulk_names: list[list[str]] = []
        for material_name in self.get_material_names():
            bulk_names.append(
                self._reactor.surface_chemistry.get_bulk_symbols(material_name)
            )
        return bulk_names

    def get_total_site_species(self) -> int:
        """Return total number of site species across all materials."""
        return self._reactor.surface_chemistry.total_site_species

    def get_total_bulk_species(self) -> int:
        """Return total number of bulk species across all materials."""
        return self._reactor.surface_chemistry.total_bulk_species

    @property
    def surface_ratemultiplier(self) -> float:
        """Return the global surface reaction-rate multiplier."""
        return self._reactor._surfaceratemultiplier

    @surface_ratemultiplier.setter
    def surface_ratemultiplier(self, value: float = 1.0e0):
        """Set the global surface reaction-rate multiplier."""
        validate_minimum_value(
            value=value,
            minimum=0.0,
            message="reaction rate multiplier must >= 0.",
        )
        if self.require_surface_chemistry():
            self._reactor._surfaceratemultiplier = value

    def check_material_area_fraction(self, mat_name: str) -> int:
        """Return whether a material-specific area fraction is set."""
        status = -1
        if self.require_surface_chemistry() and self._reactor._numb_surf_area_set > 0:
            index = self._reactor.surface_chemistry.check_surface_material(mat_name)
            if index >= 0:
                try:
                    status = self._reactor._surf_area_set.index(index)
                except ValueError:
                    status = -1
        return status + 1

    def get_material_area_fraction(self, mat_name: str) -> float:
        """Return the area fraction for a given surface material."""
        if not self.require_surface_chemistry():
            return 0.0

        index = self._reactor.surface_chemistry.check_surface_material(mat_name)
        if index < 0:
            return 0.0

        try:
            loc = self._reactor._surf_area_set.index(index)
            return self._reactor._material_area_fraction[loc]
        except ValueError:
            msg = [
                Color.YELLOW,
                "area fraction of material",
                mat_name.rstrip(),
                "has not been assigned",
                Color.END,
            ]
            log_info_message(msg)
            return 1.0e0 / float(self._reactor.numbmaterials)

    def set_material_area_fraction(self, mat_name: str, fraction: float):
        """Set the area fraction for a given surface material."""
        validate_minimum_value(
            value=fraction,
            minimum=0.0,
            message="material area fraction must >= 0.",
        )

        if self.require_surface_chemistry():
            index = self._reactor.surface_chemistry.check_surface_material(mat_name)
            if index >= 0:
                if index not in self._reactor._surf_area_set:
                    self._reactor._numb_surf_area_set += 1
                    self._reactor._surf_area_set.append(index)
                    self._reactor._material_area_fraction.append(fraction)
                else:
                    loc = self._reactor._surf_area_set.index(index)
                    self._reactor._material_area_fraction[loc] = fraction

    def check_material_temperature(self, mat_name: str) -> int:
        """Return whether a material-specific temperature is set."""
        return self._reactor.surface_chemistry.check_material_temperature(
            mat_name=mat_name
        )

    def get_material_temperature(self, mat_name: str) -> float:
        """Return temperature for a given surface material."""
        if not self.require_surface_chemistry():
            return 0.0

        temp = self._reactor.surface_chemistry.get_material_temperature(mat_name)
        if temp < 0.0:
            temp = self._reactor.temperature
        return temp

    def set_material_temperature(self, mat_name: str, temp: float):
        """Set temperature for a given surface material."""
        validate_minimum_value(
            value=temp,
            minimum=100.0,
            message="material temperature must >= 100. [K]",
        )

        if self.require_surface_chemistry():
            self._reactor.surface_chemistry.set_material_temperature(mat_name, temp)

    def get_all_surface_temperature(self) -> npt.NDArray[np.double]:
        """Return temperatures of all materials in a flat array."""
        n_material = self._reactor.surface_chemistry.number_materials
        surf_temp = np.zeros(n_material, dtype=np.double)
        for i, material_name in enumerate(
            self._reactor.surface_chemistry.get_material_names()
        ):
            surf_temp[i] = self.get_material_temperature(material_name)
        return surf_temp

    def get_site_fraction(self, mat_name: str) -> npt.NDArray[np.double]:
        """Return site fractions for a surface material."""
        return self._reactor.surface_chemistry.get_site_frac(mat_name)

    def set_site_fraction(self, mat_name: str, recipe: list[tuple[str, float]]):
        """Set site fractions for a surface material."""
        self._reactor.surface_chemistry.set_site_frac(mat_name, recipe)

    def get_all_site_fractions(self) -> npt.NDArray[np.double]:
        """Return site fractions for all materials as one array."""
        return self._reactor.surface_chemistry.get_all_site_frac()

    def get_bulk_activity(self, mat_name: str) -> npt.NDArray[np.double]:
        """Return bulk activities for a surface material."""
        return self._reactor.surface_chemistry.get_bulk_frac(mat_name)

    def set_bulk_activity(self, mat_name: str, recipe: list[tuple[str, float]]):
        """Set bulk activities for a surface material."""
        self._reactor.surface_chemistry.set_bulk_frac(mat_name, recipe)

    def get_all_bulk_growth_rates(self) -> npt.NDArray[np.double]:
        """Return bulk growth rates for all materials as one array."""
        return self._reactor.surface_chemistry.get_all_bulk_growth_rates()

    def set_init_surface_coverage(
        self,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Build default/initial site-fraction and bulk-activity arrays."""
        s_count = 0
        b_count = 0
        mat_sites: list[int] = []
        mat_bulks: list[int] = []

        for material in self._reactor.surface_chemistry.materials.values():
            isite = material.number_site_species
            ibulk = material.number_bulk_species
            mat_sites.append(isite)
            mat_bulks.append(ibulk)
            s_count += isite
            b_count += ibulk

        s_init = np.zeros(max(1, s_count), dtype=np.double)
        b_init = np.zeros(max(1, b_count), dtype=np.double)

        is_end = 0
        ib_end = 0
        for i, material_name in enumerate(
            self._reactor.surface_chemistry.material_names
        ):
            is_start = is_end
            is_end += mat_sites[i]
            if mat_sites[i] > 0:
                if self._reactor.surface_chemistry.check_site_frac_set(material_name):
                    s_init[is_start:is_end] = (
                        self._reactor.surface_chemistry.get_site_frac(material_name)
                    )
                else:
                    s_init[is_start:is_end] = 1.0 / mat_sites[i]

            ib_start = ib_end
            ib_end += mat_bulks[i]
            if mat_bulks[i] > 0:
                if self._reactor.surface_chemistry.check_bulk_act_set(material_name):
                    b_init[ib_start:ib_end] = (
                        self._reactor.surface_chemistry.get_bulk_frac(material_name)
                    )
                else:
                    b_init[ib_start:ib_end] = 1.0

        return s_init, b_init

    def set_surface_chemistry_keywords(self):
        """Emit surface-chemistry related keyword lines via host.setkeyword."""
        self._reactor.setkeyword(key="SFAC", value=self.surface_ratemultiplier)

        for material_name in self._reactor.surface_chemistry.get_material_names():
            mat_tag = self._FOUR_SPACES + material_name

            if self.check_material_temperature(mat_name=material_name) > 0:
                temp = self.get_material_temperature(mat_name=material_name)
                self._reactor.setkeyword(key="TSRF" + mat_tag, value=temp)

            if self.check_material_area_fraction(mat_name=material_name) > 0:
                area_frac = self.get_material_area_fraction(mat_name=material_name)
                self._reactor.setkeyword(key="AFRA" + mat_tag, value=area_frac)

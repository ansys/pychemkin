# Copyright (c) 2026 Synopsys, Inc. and ANSYS, Inc. and/or its affiliates.
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

"""Surface chemistry behavior mixed into the Chemistry class."""

import copy
import ctypes
from ctypes import POINTER, c_char_p
from typing import Any, Protocol

from ansys.chemkin.core import chemkin_wrapper as ck_wrapper
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.surface_components import Material
from ansys.chemkin.core.utilities import (
    log_error_message,
    log_info_message,
)

_symbol_length = 16  # Chemkin element/species symbol length
MAX_SPECIES_LENGTH = _symbol_length + 1  # Chemkin element/species symbol length + 1
LP_c_char = ctypes.POINTER(ctypes.c_char)  # pointer to C type character array


class _ValueHolder(Protocol):
    """Protocol for ctypes scalar wrappers exposing a ``value`` attribute."""

    value: int


class _SurfaceChemistryHost(Protocol):
    """Structural contract implemented by ``Chemistry`` for this mixin."""

    label: str
    matsymbol: list[str]
    material_map: dict[str, int]
    materials: dict[str, Any]
    numb_mat: int
    _materialnamedone: int
    _chemset_index: Any
    _num_materials: _ValueHolder
    _num_gas_species: _ValueHolder
    _num_max_site_species: _ValueHolder
    _num_max_bulk_species: _ValueHolder
    _num_gas_reactions: _ValueHolder
    _num_max_surf_reactions: _ValueHolder
    _num_max_phases: _ValueHolder
    max_total_species: int
    max_total_reactions: int

    def verify_surface_mechanism(self) -> bool: ...
    def no_surface_mechanism_declaration(self) -> None: ...

    @property
    def material_names(self) -> list[str]: ...


class SurfaceChemistryMixin:
    """Surface chemistry related methods for the Chemistry class."""

    def set_surface_chemistry(self: _SurfaceChemistryHost):
        """Set up the surface chemistry data."""
        """
        Set up the surface chemistry data of the Chemistry Set.
        Get the dimensions of the surface phases, species,
        and reactions. Instantiate the surface Material objects.
        """
        # number of surface materials of the Chemistry Set
        self.numb_mat = self._num_materials.value
        # dimension for variable arrays related to all species (of any material)
        # in the Chemistry Set
        self.max_total_species = (
            self._num_gas_species.value
            + self._num_max_site_species.value
            + self._num_max_bulk_species.value
        )
        # dimension for variable arrays related to all reactions (of any material)
        # in the Chemistry Set
        self.max_total_reactions = (
            self._num_gas_reactions.value + self._num_max_surf_reactions.value
        )
        # set up surface materials
        # initialization
        self.material_map.clear()  # surface material map {str, int}
        self.materials.clear()  # surface material objects {str, Material object}
        self.matsymbol.clear()  # surface material names [str]
        self._materialnamedone = 0
        # get material names
        _ = self.material_names
        # set up surface materials
        for i, key in enumerate(self.matsymbol):
            # material index is 0-based
            self.material_map[key] = i
            m = Material(mat_index=i, chem=self, label=key)
            self.materials[key] = copy.deepcopy(m)

    def no_surface_mechanism_declaration(self: _SurfaceChemistryHost):
        """Inform users the Chemistry Set does not have surface chemistry."""
        msg = [
            Color.YELLOW,
            "Chemsitry Set",
            self.label,
            "does not contain surface chemistry.",
            Color.END,
        ]
        log_info_message(msg)

    @property
    def number_materials(self: _SurfaceChemistryHost) -> int:
        """Get the number of surface materials."""
        """
        Get the number of surface materials in the Chemistry Set.

        Return
        ------
            numb_materials: number of surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_materials.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_number_phases(self: _SurfaceChemistryHost) -> int:
        """Get the maximum number of phases."""
        """
        Get the maximum number of phases among the surface materials
        of the chemistry set including the gas phase.
        A surface mechanism may contain one or more surface materials, and the
        maximum number of the surface phases (site phases + bulk phases)
        among the materials plus the gas phase will be returned.

        Return
        ------
            max_numb_phases: integer
                the maximum number of phases (including the gas phase)
                of the Chemistry Set
        """
        if self.verify_surface_mechanism():
            return self._num_max_phases.value
        else:
            self.no_surface_mechanism_declaration()
            return 1

    @property
    def max_number_sites(self: _SurfaceChemistryHost) -> int:
        """Get the maximum number of surface site species."""
        """
        Get the maximum number of surface site species among the
        surface materials of the surface chemistry.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned.

        Return
        ------
            max_numb_site_species: integer
                the maximum number of site species of the surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_max_site_species.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_number_bulks(self: _SurfaceChemistryHost) -> int:
        """Get the maximum number of surface bulk species."""
        """
        Get the maximum number of surface bulk species among the
        surface materials of the surface chemistry.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned.

        Return
        ------
            max_numb_bulk_species: integer
                the maximum number of bulk species of the surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_max_bulk_species.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_total_number_species(self: _SurfaceChemistryHost) -> int:
        """Get the maximum number of all species."""
        """
        Get the maximum number of all species (gas, site and bulk) among the surface
        materials of the Chemistry Set. This number
        includes the gas phase species, surface site species, and the bulk species.
        A surface mechanism may contain one or more surface materials, and the
        maximum total species number among the materials will be returned.
        This value can be used to size the arrays holding species
        mole fractions/activities, thermodynamic properties on any material.

        Return
        ------
            max_total_species: integer
                the maximum number of all species among the surface materials
                in the Chemistry Set.
        """
        if self.verify_surface_mechanism():
            return self.max_total_species
        else:
            self.no_surface_mechanism_declaration()
            return self._num_gas_species.value

    @property
    def max_number_surface_reactions(self: _SurfaceChemistryHost) -> int:
        """Get the maximum number of surface reactions."""
        """
        Get the maximum number of surface reactions among the
        surface materials of the surface chemistry.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned. This value can be used
        to size the arrays holding surface reaction rates on any material.

        Return
        ------
            max_numb_surface_rxn: integer
                the maximum number of surface reactions of the surface materials
        """
        if self.verify_surface_mechanism():
            return self._num_max_surf_reactions.value
        else:
            self.no_surface_mechanism_declaration()
            return 0

    @property
    def max_total_number_reactions(self: _SurfaceChemistryHost) -> int:
        """Get the maximum number of gas and surface reactions."""
        """
        Get the maximum number of all (gas and surface) reactions among the
        surface materials in the Chemistry Set.
        A surface mechanism may contain one or more surface materials, and the
        maximum number among the materials will be returned. This value can be used
        to size the arrays holding all reaction rates on any material.

        Return
        ------
            max_total_rxn: integer
                the maximum number of surface reactions of the surface materials
        """
        if self.verify_surface_mechanism():
            return self.max_total_reactions
        else:
            self.no_surface_mechanism_declaration()
            return self._num_gas_reactions.value

    @property
    def material_names(self: _SurfaceChemistryHost) -> list[str]:
        """Get list of surface material names."""
        """
        Get list of surface material names of the Chemistry Set.

        Returns
        -------
            matsymbol: list of strings
                list of surface material names in the surface mechanism
        """
        if not self.verify_surface_mechanism():
            # no surface chemistry data
            self.no_surface_mechanism_declaration()
            return []
        # process surface chemistry data
        if self._materialnamedone == 0:
            # recycle existing data
            buff_m = (LP_c_char * self._num_materials.value)()
            char_buffers = []
            for i in range(0, self._num_materials.value):
                char_buf = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
                char_buffers.append(char_buf)
                buff_m[i] = ctypes.cast(char_buf, LP_c_char)
            pp_m = ctypes.cast(buff_m, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetMaterialNames(self._chemset_index, pp_m)
            if ierr == 0:
                self.material_map.clear()
                for index in range(0, len(buff_m)):
                    mat_bytes = ctypes.cast(buff_m[index], c_char_p).value
                    mat_val = "" if mat_bytes is None else mat_bytes.decode()
                    self.material_map[mat_val.rstrip()] = index
                self._materialnamedone = 1
            else:
                # failed to get species symbols
                msg = [Color.PURPLE, "failed to get surface material names.", Color.END]
                log_error_message(msg)
                exit()
            del buff_m
            del char_buffers
        self.matsymbol[:] = self.material_map
        return self.matsymbol

    def get_material_index(self: _SurfaceChemistryHost, mat_name: str) -> int:
        """Get the surface material index."""
        """
        Get the surface material index by its name.

        Parameters
        ----------
            mat_name: str
                surface material name

        Returns
        -------
            index: integer
                index of the named surface material
        """
        index = -1
        if not self.verify_surface_mechanism():
            # no surface chemistry data
            self.no_surface_mechanism_declaration()
        #
        if mat_name in self.matsymbol:
            index = self.material_map.get(mat_name, -1)
        return index

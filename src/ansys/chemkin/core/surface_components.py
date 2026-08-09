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

"""Surface material and phase data models."""

import copy
import ctypes
from ctypes import POINTER, c_char_p, c_double, c_int
from typing import TYPE_CHECKING, Union

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core import chemkin_wrapper as ck_wrapper
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.constants import P_ATM
from ansys.chemkin.core.logger import logger

if TYPE_CHECKING:
    from ansys.chemkin.core.chemistry import Chemistry

_symbol_length = 16  # Chemkin element/species symbol length
MAX_SPECIES_LENGTH = _symbol_length + 1  # Chemkin element/species symbol length + 1
LP_c_char = ctypes.POINTER(ctypes.c_char)  # pointer to C type character array


def set_temp_array(
    temp: float,
    temp_ele: Union[float, None] = None,
    temp_ion: Union[float, None] = None,
) -> npt.NDArray[np.double]:
    """Set up the temperature array for Chemkin-CFD-API calls."""
    # set electron temperature [K]
    if temp_ele is None:
        temp_e = temp
    else:
        temp_e = temp_ele
    # set ion temperature [K]
    if temp_ion is None:
        temp_i = temp
    else:
        temp_i = temp_ion
    # set the temperature used in Chemkin-CFD-API calls
    t_array = np.ascontiguousarray([temp, temp_e, temp_i], dtype=np.double)
    return t_array


class Material:
    """Surface Chemkin material object."""

    def __init__(self, mat_index: int, chem: "Chemistry", label: str = ""):
        """Create a surface material object."""
        """
        Create a Surface Chemkin material object which may contain multiple
        surface site and bulk phases and reactions that take place on
        these surface phases.
        """
        # check
        if not chem.verify_surface_mechanism():
            msg = [
                Color.RED,
                "The Chemistry set does not contain surface chemistry.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            exit()
        # chemistry set index
        self._chemset_index = chem.chemid
        self._material_index = mat_index
        self._dispersed = False
        self._coal = False
        # material name
        if len(label) == 0:
            self.label = "Material_" + str(mat_index)
        else:
            self.label = label
        # dimensions
        chemset_index = c_int(self._chemset_index)
        material_id = c_int(mat_index)
        self.num_element = chem.mm
        self.num_gas_species = chem.kk
        self.num_gas_reactions = chem.ii_gas
        num_sites = c_int(0)
        num_bulks = c_int(0)
        num_phases = c_int(0)
        num_surf_rxns = c_int(0)
        ierr = ck_wrapper.chemkin.KINGetMaterialSizes(
            chemset_index,
            material_id,
            num_sites,
            num_bulks,
            num_phases,
            num_surf_rxns,
        )
        if ierr == 0:
            # number of phases on this material (including the gas phase)
            self.num_phases = num_phases.value
            # number of site phases on this material
            self.num_site_phase = 0
            # number of bulk phases on this material
            self.num_bulk_phase = 0
            # number of site species on this material
            self.num_site_species = num_sites.value
            # number of bulk species on this material
            self.num_bulk_species = num_bulks.value
            # number of surface reactions on this material
            self.num_surf_reactions = num_surf_rxns.value
            # including the gas phase
            self.total_species = chem.kk + self.num_site_species + self.num_bulk_species
            self.total_reactions = self.num_gas_reactions + self.num_surf_reactions
            #
            self.surf_species_map: dict[str, int] = {}
        else:
            # failed to obtain material information from chemkin
            self.num_phases = 0
            self.num_site_phase = 0
            self.num_bulk_phase = 0
            self.num_site_species = 0
            self.num_bulk_species = 0
            self.num_surf_reactions = 0
            self.total_species = 0
            self.total_reactions = 0
            #
            msg = [Color.PURPLE, "failed to process surface material data.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # process additional material data
        ibase0 = c_int(0)  # array index is 1-based
        # mapping of the surface species symbols (for all site and bulk species)
        # and the global index (starting from the gas species) on this material
        # {species_symbol, global_index}
        self.site_species_map: dict[str, int] = {}
        self.bulk_species_map: dict[str, int] = {}
        # get element symbols from the chemistry set
        self.ESymbols = chem.element_symbols
        # get number of phases on this material
        num_site_phase = c_int(0)
        num_bulk_phase = c_int(0)
        nphase = c_int(0)
        ierr = ck_wrapper.chemkin.KINGetPhaseSizes(
            chemset_index,
            material_id,
            nphase,
            num_site_phase,
            num_bulk_phase,
        )
        self.num_site_phase = num_site_phase.value
        self.num_bulk_phase = num_bulk_phase.value
        self.first_site_phase = 0
        self.last_site_phase = 0
        self.first_bulk_phase = 0
        self.last_bulk_phase = 0
        if self.num_site_phase > 0:
            self.first_site_phase = 1
            self.last_site_phase = self.first_site_phase + self.num_site_phase - 1
        if self.num_bulk_phase > 0:
            self.first_bulk_phase = self.last_site_phase + 1
            self.last_bulk_phase = self.first_bulk_phase + self.num_bulk_phase - 1
        # get number of surface/bulk species and their indices on this material
        self.nspecies_of_phase = np.zeros(self.num_phases, dtype=np.int32)
        self.first_species_id_of_phase = np.zeros(self.num_phases, dtype=np.int32)
        self.last_species_id_of_phase = np.zeros(self.num_phases, dtype=np.int32)
        ierr = ck_wrapper.chemkin.KINGetPhaseSpeciesIndices(
            chemset_index,
            material_id,
            self.nspecies_of_phase,
            self.first_species_id_of_phase,
            self.last_species_id_of_phase,
            ibase0,
        )
        # get number of surface reactions on this material
        # get phase names on this material
        self._phasenamesdone = 0
        self.PhaseSymbol: list[str] = []
        self.phase_map: dict[str, int] = {}
        self.phases: dict[str, SurfacePhase] = {}
        phase_type = ""
        # instantiate surface phases
        self.phase_names
        # global index of the site species on this material
        site_global_index = self.num_gas_species
        # global index of the bulk species on this material
        bulk_global_index = self.num_gas_species + self.num_site_species
        # surface species symbols
        # flag indicating the status of species setup of the phase
        self._sitespeciesnamesdone = 0
        self._bulkspeciesnamesdone = 0
        # site species symbols of this material
        self.site_symbols: list[str] = []
        # bulk species symbols of this material
        self.bulk_symbols: list[str] = []
        # species composition flag
        self.ncf_done = 0
        # surface species molecular weight arrays
        self.sitewt_done = 0
        self.bulkwt_done = 0
        self.site_wt = np.zeros(self.num_site_species, dtype=np.double)
        self.bulk_wt = np.zeros(self.num_bulk_species, dtype=np.double)
        #
        if self.num_site_phase > 0:
            # get site species symbols of this material
            self.site_symbols = self.site_species_names
            if len(self.site_symbols) > 0:
                self._sitespeciesnamesdone = 1
                for i, name in enumerate(self.site_symbols):
                    self.site_species_map[name] = site_global_index + i
            # compute site (phase) density [mole/cm2]
            self.site_density = self.get_site_density()
            # get site species molecular weights [g/mole]
            self.site_wt = self.get_site_molar_weights()
        else:
            self.site_density = np.zeros(0, dtype=np.double)
        #
        if self.num_bulk_phase > 0:
            # bulk species symbols of this material
            self.bulk_symbols = self.bulk_species_names
            if len(self.bulk_symbols) > 0:
                self._bulkspeciesnamesdone = 1
                for i, name in enumerate(self.bulk_symbols):
                    self.bulk_species_map[name] = bulk_global_index + i
            # compute bulk (species) density [g/cm3]
            self.bulk_density = self.get_bulk_density()
            # get site species molecular weights [g/mole]
            self.bulk_wt = self.get_bulk_molar_weights()
        else:
            self.bulk_density = np.zeros(0, dtype=np.double)
        # set up surface phases on this materials
        sp_count = 0
        bp_count = 0
        # loop over all phases
        for i, key in enumerate(self.PhaseSymbol):
            # number of species on this phase
            nspec = self.nspecies_of_phase[i]
            # phase index in Chemkin is 1-based
            # phase #0 is the gas phase
            if i == 0:
                # skip the gas phase
                continue
            elif i <= self.last_site_phase:
                # site phase
                phase_type = "site"
                site_global_index += nspec
                sp_count += 1
            else:
                # bulk phase
                phase_type = "bulk"
                bulk_global_index += nspec
                bp_count += 1
            # instantiate the surface phase
            s = SurfacePhase(phase_index=i, phase_type=phase_type, label=key)
            # set the number of species on this surface phase
            s.number_species = self.nspecies_of_phase[i]
            # set (1-based) global index of the first species
            s.first_species_index = self.first_species_id_of_phase[i]
            # set (1-based) global index of the last species
            s.last_species_index = self.last_species_id_of_phase[i]
            if phase_type == "site":
                # set site phase density [mole/cm2]
                s.site_density = self.site_density[i]
            # add to the phases dictionary
            self.phases[key] = copy.deepcopy(s)
            #
            del s

        # additional surface material properties
        self.surftemp = 300.0
        self.activearea = 0.0

    @property
    def material_index(self) -> int:
        """Get the surface material index."""
        return self._material_index

    @property
    def number_phases(self) -> int:
        """Get the number of phases."""
        """
        Get the number of phases (including the gas phase) of the material.

        Returns
        -------
            numb_phases: integer
                number of phases of the material

        """
        return self.num_phases

    @property
    def number_site_phases(self) -> int:
        """Get the number of site phases."""
        """
        Get the number of site phases of the material.

        Returns
        -------
            numb_site_phase: integer
                number of site phases of the material

        """
        return self.num_site_phase

    @property
    def number_bulk_phases(self) -> int:
        """Get the number of bulk phases."""
        """
        Get the number of bulk phases of the material.

        Returns
        -------
            numb_bulk_phase: integer
                number of bulk phases of the material

        """
        return self.num_bulk_phase

    @property
    def number_total_species(self) -> int:
        """Get the number of species in all phases."""
        """
        Get the number of species in all phases (gas, site, and bulk).

        Returns
        -------
            total_species: integer
                total number of species associated with the material
        """
        return self.total_species

    @property
    def number_site_species(self) -> int:
        """Get the total number of site species."""
        """
        Get the total number of site species (of all site phases) of the material.

        Returns
        -------
            numb_site_species: integer
                total number of site species of the material

        """
        return self.num_site_species

    @property
    def number_bulk_species(self) -> int:
        """Get the total number of bulk species."""
        """
        Get the total number of bulk species (of all bulk phases) of the material.

        Returns
        -------
            numb_bulk_species: integer
                total number of bulk species of the material

        """
        return self.num_bulk_species

    @property
    def first_site_phase_index(self) -> int:
        """Get phase index of the first site phase."""
        """
        Get phase index of the first site phase of the material.

        Returns
        -------
            first_site_phase: integer
                phase index (0-based) of the first site phase of the material

        """
        if self.number_site_phases > 0:
            return self.first_site_phase
        else:
            return -1

    @property
    def last_site_phase_index(self) -> int:
        """Get phase index of the last site phase."""
        """
        Get phase index of the last site phase of the material.

        Returns
        -------
            last_site_phase: integer
                phase index (0-based) of the last site phase of the material

        """
        if self.number_site_phases > 0:
            return self.last_site_phase
        else:
            return -1

    @property
    def first_bulk_phase_index(self) -> int:
        """Get phase index of the first bulk phase."""
        """
        Get phase index of the first bulk phase of the material.

        Returns
        -------
            first_bulk_phase: integer
                phase index (0-based) of the first bulk phase of the material

        """
        if self.number_bulk_phases > 0:
            return self.first_bulk_phase
        else:
            return -1

    @property
    def last_bulk_phase_index(self) -> int:
        """Get phase index of the last bulk phase."""
        """
        Get phase index of the last bulk phase of the material.

        Returns
        -------
            last_bulk_phase: integer
                phase index (0-based) of the last bulk phase of the material

        """
        if self.number_bulk_phases > 0:
            return self.last_bulk_phase
        else:
            return -1

    @property
    def number_surface_reactions(self) -> int:
        """Get the number of surface reactions."""
        """
        Get the number of surface reactions of the material.

        Returns
        -------
            ii_surf: integer
                number of surface reactions of the material

        """
        return self.num_surf_reactions

    @property
    def number_total_reactions(self) -> int:
        """Get the number of all reactions."""
        """
        Get the number of all reactions (gas and surface) of the material.

        Returns
        -------
            ii_total: integer
                total number of reactions associated with the material

        """
        return self.total_reactions

    def information(self):
        """List the size information of the material."""
        # surface material information
        print(f"Material Name: {self.label}")
        # gas-phase information
        print("----")
        print(f"  Number of gas species = {self.num_gas_species}")
        print(f"  Number of gas-phase reactions = {self.num_gas_reactions}")
        # surface phases
        print("----")
        print(f"  Number of surface phases = {self.num_phases - 1}")
        # site phase information
        print("----")
        print(f"  Number of Site phases = {self.num_site_phase}")
        for n, p in enumerate(self.PhaseSymbol):
            if n == 0:
                # skip the gas phase
                continue
            if self.phases[p].phase_type == "site":
                print(f"    Site phase: {n} {p} global ID: {self.phase_map[p]}")
        # site species
        print(f"  Number of Site species on phase = {self.num_site_species}")
        if self.num_site_species > 0:
            print("----")
            for k, map in enumerate(self.site_species_map.items()):
                print(
                    f"    Site species index = {k} {map[0]} "
                    f"global ID: {map[1]} WT = {self.site_wt[k]}"
                )
        # bulk phase information
        print("----")
        print(f"  Number of Bulk phases = {self.num_bulk_phase}")
        for n, p in enumerate(self.PhaseSymbol):
            if n == 0:
                # skip the gas phase
                continue
            if self.phases[p].phase_type == "bulk":
                print(f"    Bulk phase: {n} {p} global ID: {self.phase_map[p]}")
        # bulk species
        print(f"  Number of Bulk species = {self.num_bulk_species}")
        if self.num_bulk_species > 0:
            print("----")
            for k, map in enumerate(self.bulk_species_map.items()):
                print(
                    f"    Bulk species index = {k} {map[0]} "
                    f"global ID: {map[1]} WT = {self.bulk_wt[k]}"
                )
        # surface reaction information
        print("----")
        print(f"  Number of surface reactions = {self.num_surf_reactions}")

    def get_site_specindex(self, symbol: str = "") -> tuple[int, int]:
        """Get the local and the global indices of a site species."""
        """
        Get the phase-wise and the global (material-wise) indices of a site species.

        Parameters
        ----------
            symbol: string
                species symbol

        Returns
        -------
            local_index: integer
                (0-based) species index of the phase
            global_index: integer
                (0-based) global index of the matieral (including the gas phase)
        """
        global_index = -1
        local_index = -1
        if self.num_site_species > 0:
            if len(symbol) > 0:
                global_index = self.site_species_map.get(symbol, -1)
                if global_index > 0:
                    first_id = self.first_site_species_index
                    local_index = global_index - first_id
        return global_index, local_index

    def get_bulk_specindex(self, symbol: str = "") -> tuple[int, int]:
        """Get the local and the global indices of a bulk species."""
        """
        Get the phase-wise and the global (material-wise) indices of a bulk species.

        Parameters
        ----------
            symbol: string
                species symbol

        Returns
        -------
            local_index: integer
                (0-based) species index of the phase
            global_index: integer
                (0-based) global index of the matieral (including the gas phase)
        """
        global_index = -1
        local_index = -1
        if self.num_bulk_species > 0:
            if len(symbol) > 0:
                global_index = self.bulk_species_map.get(symbol, -1)
                if global_index > 0:
                    first_id = self.first_bulk_species_index
                    local_index = global_index - first_id
        return global_index, local_index

    def get_surf_specindex(
        self, symbol: str = "", mode: str = "normal"
    ) -> tuple[int, int, str]:
        """Get the local and the global indices and type of a surface species."""
        """
        Get the phase-wise and the global (material-wise) species indices and
        the phase type of the surface species given by the symbol.

        Parameters
        ----------
            symbol: string
                species symbol
            mode: str, {'normal', 'silent'}
                "silent": no print out

        Returns
        -------
            global_index: integer
                (0-based) global index of the matieral (including the gas phase)
            local_index: integer
                (0-based) species index of the phase
            species_type: string, ['site', 'bulk']
                type of the species
        """
        species_type = ""
        # try site phase species
        global_index, local_index = self.get_site_specindex(symbol)
        if global_index >= 0:
            species_type = "site"
        else:
            # try bulk phase species
            global_index, local_index = self.get_bulk_specindex(symbol)
            if global_index >= 0:
                species_type = "bulk"
            else:
                if mode.lower() == "normal":
                    # failed to find the species
                    msg = [
                        Color.PURPLE,
                        "failed to find species",
                        symbol,
                        "on surface material",
                        self.label,
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                else:
                    return -1, -1, ""
        return global_index, local_index, species_type

    @property
    def phase_names(self):
        """Get list of phase names."""
        """
        Get list of phase names on the surface material.

        Returns
        -------
            Phasesymbol: list of strings
                list of phase names on the surface material
        """
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        if self._phasenamesdone == 0:
            # recycle existing data
            buff_p = (LP_c_char * self.num_phases)()
            for i in range(0, self.num_phases):
                buff_p[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_p = ctypes.cast(buff_p, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetPhaseNames(chem_id, mat_id, pp_p)
            if ierr == 0:
                self.phase_map.clear()
                for index in range(0, len(buff_p)):
                    str_val = ctypes.cast(buff_p[index], c_char_p).value.decode()
                    self.phase_map[str_val] = index
                self._phasenamesdone == 1
            else:
                # failed to get species symbols
                msg = [Color.PURPLE, "failed to get phase names.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_p
        # convert string type
        mylist = list(self.phase_map.keys())
        self.PhaseSymbol.clear()
        for s in mylist:
            self.PhaseSymbol.append(s)
        del mylist
        return self.PhaseSymbol

    @property
    def site_species_names(self):
        """Set list of site species names."""
        """
        Set list of site species names on the surface phase.

        Returns
        -------
            Phasesymbol: list of strings
                list of site species names on the surface site phase
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        #
        if self._sitespeciesnamesdone == 0:
            buff_s = (LP_c_char * self.num_site_species)()
            for i in range(0, self.num_site_species):
                buff_s[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_s = ctypes.cast(buff_s, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetSiteNames(chem_id, mat_id, pp_s)
            if ierr == 0:
                site_symbols = []
                for index in range(0, len(buff_s)):
                    s_val = ctypes.cast(buff_s[index], c_char_p).value.decode()
                    site_symbols.append(s_val)
            else:
                # failed to get surface site species symbols
                msg = [
                    Color.PURPLE,
                    "failed to get surface site species names.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_s
        else:
            site_symbols = self.site_symbols
        return site_symbols

    @property
    def bulk_species_names(self):
        """Set list of bulk species names."""
        """
        Set list of bulk species names on the surface phase.

        Returns
        -------
            Phasesymbol: list of strings
                list of bulk species names on the surface bulk phase
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        #
        if self._bulkspeciesnamesdone == 0:
            buff_b = (LP_c_char * self.num_bulk_species)()
            for i in range(0, self.num_bulk_species):
                buff_b[i] = ctypes.create_string_buffer(MAX_SPECIES_LENGTH)
            pp_b = ctypes.cast(buff_b, POINTER(LP_c_char))
            ierr = ck_wrapper.chemkin.KINGetBulkNames(chem_id, mat_id, pp_b)
            if ierr == 0:
                bulk_symbols = []
                for index in range(0, len(buff_b)):
                    b_val = ctypes.cast(buff_b[index], c_char_p).value.decode()
                    bulk_symbols.append(b_val)
            else:
                # failed to get surface site species symbols
                msg = [
                    Color.PURPLE,
                    "failed to get surface bulk species names.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            del buff_b, pp_b
        else:
            bulk_symbols = self.bulk_symbols
        return bulk_symbols

    def get_site_density(self) -> npt.NDArray[np.double]:
        """Get the site phase densities."""
        """
        Get the densities of all site phases of the material.

        Returns
        -------
            density: 1-D double array, dimension = number of phases of the material
                site phase density [mole/cm2]
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        site_density = np.zeros(self.num_phases, dtype=np.double)
        ierr = ck_wrapper.chemkin.KINGetSiteDensity(chem_id, mat_id, site_density)
        if ierr != 0:
            # failed to get surface site species density
            msg = [Color.PURPLE, "failed to get surface site density.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return site_density

    def get_bulk_density(self) -> npt.NDArray[np.double]:
        """Get bulk species densities."""
        """
        Get the densities of all bulk species of the material.

        Returns
        -------
            density: 1-D double array,
                dimension = number of all bulk species of the material
                bulk species density [g/cm3]
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        density = np.zeros(self.num_bulk_species, dtype=np.double)
        # "fake" pressure and temperature since density of
        # solid bulk species is ~ constant
        temp = c_double(300.0)
        pres = c_double(P_ATM)
        ierr = ck_wrapper.chemkin.KINGetBulkDensity(
            chem_id, mat_id, pres, temp, density
        )
        if ierr != 0:
            # failed to get surface bulk species density
            msg = [Color.PURPLE, "failed to get surface bulk density.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        return density

    def get_site_molar_weights(self) -> npt.NDArray[np.double]:
        """Get the site species molecular weights."""
        """
        Get the molecular weights of all site species of the material.

        Returns
        -------
            wt: 1-D double array,
                dimension = number of all site species of the material
                molecular weights [g/mole]
        """
        if self.num_site_species == 0:
            wt = np.zeros(1, dtype=np.double)
            msg = [
                Color.MAGENTA,
                "no site species on material",
                self.label,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            return wt
        if self.sitewt_done == 0:
            chem_id = c_int(self._chemset_index)
            mat_id = c_int(self._material_index)
            wt = np.zeros(self.num_site_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetSiteMolecularWeights(chem_id, mat_id, wt)
            if ierr != 0:
                # failed to get surface site species molecular weights
                msg = [
                    Color.PURPLE,
                    "failed to get site species molecular weights.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            # check for non-positive molar weights (possible for surface chemistry)
            for k, w in enumerate(wt):
                if w > 0.0:
                    continue
                msg = [
                    Color.YELLOW,
                    "site species",
                    self.site_symbols[k],
                    "has non-positive molecular weight.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            self.sitewt_done = 1
        else:
            wt = copy.deepcopy(self.site_wt)
        return wt

    def get_bulk_molar_weights(self) -> npt.NDArray[np.double]:
        """Get the bulk species molecular weights."""
        """
        Get the molecular weights of all bulk species of the material.

        Returns
        -------
            wt: 1-D double array,
                dimension = number of all site species of the material
                molecular weights [g/mole]
        """
        if self.num_bulk_species == 0:
            wt = np.zeros(1, dtype=np.double)
            msg = [
                Color.MAGENTA,
                "no bulk species on material",
                self.label,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            return wt
        if self.bulkwt_done == 0:
            chem_id = c_int(self._chemset_index)
            mat_id = c_int(self._material_index)
            wt = np.zeros(self.num_bulk_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetBulkMolecularWeights(chem_id, mat_id, wt)
            if ierr != 0:
                # failed to get bulk species molecular weights
                msg = [
                    Color.PURPLE,
                    "failed to get bulk species molecular weights.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            # check for non-positive molar weights (possible for surface chemistry)
            for k, w in enumerate(wt):
                if w > 0.0:
                    continue
                msg = [
                    Color.YELLOW,
                    "bulk species",
                    self.bulk_symbols[k],
                    "has non-positive molecular weight.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
            self.bulkwt_done = 1
        else:
            wt = copy.deepcopy(self.bulk_wt)
        return wt

    @property
    def first_site_species_index(self) -> int:
        """Get the global species index of the first site species."""
        """
        Get the global species index of the first site species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_site_species <= 0:
            return -1
        else:
            return self.site_species_map.get(self.site_symbols[0], -1)

    @property
    def first_bulk_species_index(self) -> int:
        """Get the global species index of the first bulk species."""
        """
        Get the global species index of the first bulk species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_bulk_species <= 0:
            return -1
        else:
            return self.bulk_species_map.get(self.bulk_symbols[0], -1)

    @property
    def last_site_species_index(self) -> int:
        """Get the global species index of the last site species."""
        """
        Get the global species index of the last site species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_site_species <= 0:
            return -1
        else:
            id = self.num_site_species - 1
            return self.site_species_map.get(self.site_symbols[id], -1)

    @property
    def last_bulk_species_index(self) -> int:
        """Get the global species index of the last bulk species."""
        """
        Get the global species index of the last bulk species of the material.

        Returns
        -------
            spec_index: integer
                (0-based) global species index
        """
        if self.num_bulk_species <= 0:
            return -1
        else:
            id = self.num_bulk_species - 1
            return self.bulk_species_map.get(self.bulk_symbols[id], -1)

    def get_site_occupancy(self) -> npt.NDArray[np.int32]:
        """Get site species occupancy."""
        """
        Get site species occupancy.

        Returns
        -------
            occupancy: 1-D integer array,
                dimension = total number of site species of the material
                site occupancy of the site species
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        occ = np.zeros(self.total_species, dtype=np.int32)
        ierr = ck_wrapper.chemkin.KINGetSpeciesOccupancy(chem_id, mat_id, occ)
        # get the global index of the first site species
        kstart = self.first_site_species_index
        # get the global index of the first bulk species (= last site species index + 1)
        kstop = kstart + self.num_site_species
        if kstart < 0:
            ierr = 10
        if ierr == 0:
            site_occupancy = np.zeros(self.num_site_species, dtype=np.int32)
            i = 0
            for k in range(kstart, kstop):
                site_occupancy[i] = occ[k]
                i += 1
            del occ
            return site_occupancy
        else:
            # failed to get bulk species molecular weights
            msg = [
                Color.PURPLE,
                "failed to get site species occupancy.\n",
                Color.SPACEx6,
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_site_species_h(
        self, temp: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get site species enthalpy."""
        """
        Get site species enthalpy.

        Parameters
        ----------
            temp: double, optional
                surface material temperature [K]

        Returns
        -------
            h: 1-D double array,
                dimension = total number of site species of the material
                enthalpy of the site species [ergs/mol]
        """
        if self.num_site_species <= 0:
            msg = [Color.PURPLE, "no site species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        h_all = np.zeros(self.total_species, dtype=np.double)
        #
        if temp is None:
            t = max(self.surftemp, 3.0e2)
            tt = set_temp_array(temp=t)
        else:
            tt = set_temp_array(temp=temp)
        ierr = ck_wrapper.chemkin.KINGetAllSpeciesEnthalpy(chem_id, mat_id, tt, h_all)
        if ierr == 0:
            # extract enthalpy values for the site species only [ergs/mol]
            kstart = self.first_site_species_index
            kstop = kstart + self.num_site_species
            h = h_all[kstart:kstop].copy()
            #
            return h
        else:
            # failed to compute enthalpies
            msg = [
                Color.PURPLE,
                "failed to compute site species enthalpies.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_bulk_species_h(
        self, temp: Union[float, None] = None, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get bulk species enthalpy."""
        """
        Get bulk species enthalpy.

        Parameters
        ----------
            temp: double, optional
                surface material temperature [K]
            pres: double, optional
                bulk pressure [dynes/cm2], it might be different from the gas pressure

        Returns
        -------
            h: 1-D double array,
                dimension = total number of bulk species of the material
                enthalpy of the bulk species [ergs/mol]
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        #
        if temp is None:
            t = max(self.surftemp, 3.0e2)
            tt = set_temp_array(temp=t)
        else:
            tt = set_temp_array(temp=temp)
        #
        if pres is None:
            # bulk pressure is not given
            h_all = np.zeros(self.total_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetAllSpeciesEnthalpy(
                chem_id, mat_id, tt, h_all
            )
            if ierr == 0:
                # extract enthalpy values for the bulk species only [ergs/mol]
                kstart = self.first_bulk_species_index
                kstop = kstart + self.num_bulk_species
                h = h_all[kstart:kstop].copy()
                #
                return h
            else:
                # failed to compute enthalpies
                msg = [
                    Color.PURPLE,
                    "failed to compute bulk species enthalpies.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        else:
            # bulk pressure is given
            pp = c_double(pres)
            h = np.zeros(self.num_bulk_species, dtype=np.double)
            ierr = ck_wrapper.chemkin.KINGetBulkSpeciesEnthalpy(
                chem_id, mat_id, pp, tt, h
            )
            if ierr == 0:
                # convert [ergs/gm] to [ergs/mol]
                for k in range(len(h)):
                    if self.bulk_wt[k] <= 0.0e0:
                        h[k] = 0.0e0
                    else:
                        h[k] = h[k] * self.bulk_wt[k]
                return h
            else:
                # failed to compute enthalpies
                msg = [
                    Color.PURPLE,
                    "failed to compute bulk species enthalpies.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()

    def get_bulk_species_cp(
        self, temp: Union[float, None] = None, pres: Union[float, None] = None
    ) -> npt.NDArray[np.double]:
        """Get bulk species specific heat capacity."""
        """
        Get bulk species specific heat capacity.

        Parameters
        ----------
            temp: double, optional
                surface material temperature [K]
            pres: double, optional
                bulk pressure [dynes/cm2], it might be different from the gas pressure

        Returns
        -------
            cp: 1-D double array,
                dimension = total number of bulk species of the material
                specific heat capacity of the bulk species [ergs/mol-K]
        """
        if self.num_bulk_species <= 0:
            msg = [Color.PURPLE, "no bulk species found.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        cp = np.zeros(self.num_bulk_species, dtype=np.double)
        #
        if pres is None:
            pp = c_double(P_ATM)
        else:
            pp = c_double(pres)
        #
        if temp is None:
            t = max(self.surftemp, 3.0e2)
            tt = set_temp_array(temp=t)
        else:
            tt = set_temp_array(temp=temp)
        #
        ierr = ck_wrapper.chemkin.KINGetBulkSpeciesSpecificHeat(
            chem_id, mat_id, pp, tt, cp
        )
        if ierr == 0:
            # convert [ergs/gm-K] to [ergs/mol-K]
            for k in range(len(cp)):
                if self.bulk_wt[k] <= 0.0e0:
                    cp[k] = 0.0e0
                else:
                    cp[k] = cp[k] * self.bulk_wt[k]
        else:
            # failed to compute enthalpies
            msg = [
                Color.PURPLE,
                "failed to compute bulk species specific heat capacity.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return cp

    def material_species_composition(
        self, elemindex: int = -1, specindex: int = -1
    ) -> int:
        """Get elemental composition of a species."""
        """
        Get elemental composition of a species.

        Parameters
        ----------
            elemindex: integer
                (0-based) index of the element
            specindex: integer
                (0-based) index of the species
                (all species of the material including the gas phase species)

        Returns
        -------
            count: integer
                number of the element in the given species
        """
        if self.ncf_done == 0:
            chem_id = c_int(self._chemset_index)
            mat_id = c_int(self._material_index)
            # initialize the NCF matrix
            dim = (self.num_element, self.total_species)
            self.elementalcomp = np.zeros(dim, dtype=np.int32, order="F")
            # load the NCF matrix
            ierr = ck_wrapper.chemkin.KINGetSurfaceSpeciesComposition(
                chem_id, mat_id, self.elementalcomp
            )
            if ierr != 0:
                msg = [
                    Color.PURPLE,
                    "failed to compute elemental compositions.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                self.ncf_done = 1

        # check element index
        if elemindex < 0 or elemindex >= self.num_element:
            msg = [Color.PURPLE, "element index is out of bound.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # check species index
        if specindex < 0 or specindex >= self.total_species:
            msg = [Color.PURPLE, "species index is out of bound.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        return self.elementalcomp[elemindex][specindex]

    @property
    def surface_temperature(self) -> float:
        """Get surface material temperature."""
        """
        Get the temperature of the material.

        Returns
        -------
            temp: double
                material surface temperature [K]
        """
        return self.surftemp

    @surface_temperature.setter
    def surface_temperature(self, temp: float):
        """Set surface material temperature."""
        """
        Set the surface temperature of the material.

        Parameters
        ----------
            temp: double
                material surface temperature [K]
        """
        if temp > 1.0e2:
            self.surftemp = temp
        else:
            msg = [Color.PURPLE, "surface temperature must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def surface_area(self) -> float:
        """Get surface material area."""
        """
        Get the active surface area of the material.

        Returns
        -------
            area: double
                material surface area [cm2]
        """
        return self.activearea

    @surface_area.setter
    def surface_area(self, area: float):
        """Set surface material area."""
        """
        Set the active surface area of the material.

        Parameters
        ----------
            area: double
                material surface area [cm2]
        """
        if area > 0.0e0:
            self.activearea = area
        else:
            msg = [Color.PURPLE, "surface area must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_surface_reaction_parameters(
        self,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get the Arrhenius reaction rate parameters."""
        """
        Get the Arrhenius reaction rate parameters of all surface reactions
        of the material.

        Returns
        -------
            a_factor: 1-D double array
                Arrhenius pre-exponent A-Factor of reaction in [mole-cm3-sec-K]
            beta: 1-D double array
                Arrhenius temperature exponent [-]
            act_energy: 1-D double array
                activation temperature [K]
        """
        reactionsize = self.num_surf_reactions
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # pre-exponent A factor of all surface reactions in the mechanism
        # in cgs units [mole-cm-sec-K]
        a_factor = np.zeros(shape=reactionsize, dtype=np.double)
        # temperature exponent of all reactions [-]
        beta = np.zeros_like(a_factor, dtype=np.double)
        # activation energy/temperature of all reactions [K]
        act_energy = np.zeros_like(a_factor, dtype=np.double)
        # flag for using A factor as the sticking coefficient
        is_stick = np.zeros_like(a_factor, dtype=np.int32)
        # get the reaction parameters
        ierr = ck_wrapper.chemkin.KINGetSurfaceReactionRateParameters(
            chem_id, mat_id, a_factor, beta, act_energy, is_stick
        )
        if ierr != 0:
            a_factor[:] = 0.0e0
            beta[:] = 0.0e0
            act_energy[:] = 0.0e0
        return a_factor, beta, act_energy

    def set_surface_reaction_afactor(self, reaction_index: int, a_factor: float):
        """Set the Arrhenius A-Factor of the given surface reaction."""
        """
        (Re)set the Arrhenius A-Factor of the given surface reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the surface reaction of which
                the A-Factor to be reset
            a_factor: double
                new A-Factor value in [mole-cm-sec-K]
        """
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        if a_factor < 0.0e0:
            msg = [Color.PURPLE, "A-Factor must >= 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(-reaction_index)  # negative index to "put" A-factor value
        ierr = ck_wrapper.chemkin.KINSetAFactorForASurfaceReaction(
            chem_id, mat_id, ireac, c_double(a_factor)
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to set Arrhenius A-Factor,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_surface_reaction_afactor(self, reaction_index: int) -> float:
        """Get the Arrhenius A-Factor of the given surface reaction."""
        """
        Get the Arrhenius A-Factor of the given surface reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the reaction

        Returns
        -------
            AFactor: double
                Arrhenius A-Factor of the given reaction in [mole-cm3-sec-K]
        """
        # initialization
        a_factor = c_double(0.0e0)
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        # get the A-factor value
        ierr = ck_wrapper.chemkin.KINSetAFactorForASurfaceReaction(
            chem_id, mat_id, ireac, a_factor
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find Arrhenius A-Factor,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return a_factor.value

    def set_surface_reaction_act_energy(self, reaction_index: int, act_energy: float):
        """Set the Arrhenius activation energy of the given surface reaction."""
        """
        (Re)set the Arrhenius activation energy of
        the given surface reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the surface reaction of which
                the activation energy to be reset
            act_energy: double
                new Aactivation eneergy value in [K]
        """
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # if act_energy < 0.0e0:
        #     msg = [Color.PURPLE, "activation energy must >= 0.", Color.END]
        #     this_msg = Color.SPACE.join(msg)
        #     logger.error(this_msg)
        #     exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        # use negative index to "put" activation energy value
        ireac = c_int(-reaction_index)
        ierr = ck_wrapper.chemkin.KINSetActEnergyForASurfaceReaction(
            chem_id, mat_id, ireac, c_double(act_energy)
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to set Arrhenius activation energy,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def get_surface_reaction_act_energy(self, reaction_index: int) -> float:
        """Get the Arrhenius activation energy."""
        """
        Get the Arrhenius activation energy of the given surface
        reaction of the material.

        Parameters
        ----------
            reaction_index: integer
                (1-based) index of the reaction

        Returns
        -------
            act_energy: double
                Arrhenius activation energy of the given reaction in [K]
        """
        # initialization
        act_energy = c_double(0.0e0)
        # check inputs
        if reaction_index > self.num_surf_reactions or reaction_index < 1:
            msg = [
                Color.PURPLE,
                "reaction index is out of bound,",
                "range = [1 ~ " + str(self.num_surf_reactions) + "].",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        # get the activation energy value
        ierr = ck_wrapper.chemkin.KINSetActEnergyForASurfaceReaction(
            chem_id, mat_id, ireac, act_energy
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find Arrhenius activation energy,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        return act_energy.value

    def get_surface_reaction_string(self, reaction_index: int) -> str:
        """Get the reaction string of the surface reaction."""
        """
        Get the reaction string of the surface reaction
        of the material specified by the reaction index.
        The surface reaction index starts from 1.

        Parameters
        ----------
            reaction_index: integer
                (1-based) surface reaction index

        Returns
        -------
            reactionstring: string
                reaction string of the given surface reaction
        """
        # initialization
        reactionstring = ""
        if reaction_index > self.num_surf_reactions:
            msg = [
                Color.PURPLE,
                "reaction index must be <",
                str(self.num_surf_reactions),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        elif reaction_index <= 0:
            msg = [Color.PURPLE, "reaction index must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        chem_id = c_int(self._chemset_index)
        mat_id = c_int(self._material_index)
        # convert the reaction parameters
        ireac = c_int(reaction_index)
        i_string_size = c_int(0)
        # get reaction string (might have to be increased to 2048 for 26R1)
        rstring = bytes(" " * 1024, "utf-8")
        ierr = ck_wrapper.chemkin.KINGetSurfaceReactionString(
            chem_id, mat_id, ireac, i_string_size, rstring
        )
        if ierr != 0:
            msg = [
                Color.PURPLE,
                "failed to find reaction string,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # convert C string back to string
        # print(rstring.decode()[0:iStringSize.value])
        reactionstring = rstring.decode()[0 : i_string_size.value]
        del rstring
        return reactionstring


class SurfacePhase:
    """Surface phase object."""

    def __init__(self, phase_index: int, phase_type: str, label: str = ""):
        """Create a surface phase object."""
        """The surface phase can be either a site phase or a bulk phase."""
        self.phase_index = phase_index
        if len(label) > 0:
            self.label = label
        else:
            self.label = "phase_" + str(phase_index)
        if phase_type in ["site", "bulk"]:
            self.phase_type = phase_type
        else:
            msg = [
                Color.PURPLE,
                "invalid phase type,",
                "phase type must be either 'site' or 'bulk'",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        self._numb_species = 0
        self.first_spec_index = 0
        self.last_spec_index = 0
        self._density = 0.0e0  # site density [mole/cm2] for site phases only

    @property
    def number_species(self) -> int:
        """Get the number of species."""
        """
        Get the number of species on this surface phase.

        Returns
        -------
            numb_species: integer
                number of species
        """
        return self._numb_species

    @number_species.setter
    def number_species(self, count: int):
        """Set the number of species."""
        """
        Set the number of species on this surface phase.

        Parameters
        ----------
            count: integer
                number of species
        """
        if count > 0:
            self._numb_species = count
        else:
            msg = [Color.PURPLE, "number of species must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def site_density(self) -> float:
        """Get surface site density."""
        """
        Get surface site density if this is a surface site phase.

        Returns
        -------
            density: double
                site phase density [mole/cm2]
        """
        return self._density

    @site_density.setter
    def site_density(self, den: float):
        """Set the site density."""
        """
        Set the site density if this is a surface site phase.

        Parameters
        ----------
            den: double
                site phase density [mole/cm2]
        """
        if den > 0.0e0:
            self._density = den
        else:
            msg = [Color.PURPLE, "surface site phase density must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def first_species_index(self) -> int:
        """Get the global index of the first species."""
        """
        Get the global index (1-based) of the first species of the phase.

        Returns
        -------
            global_index: integer
                the (1-based) global index of the first species
        """
        return self.first_spec_index

    @first_species_index.setter
    def first_species_index(self, id: int):
        """Set the global index of the first species."""
        """
        Set the global index (1-based) of the first species of the surface phase.

        Parameters
        ----------
            id: integer
                the (1-based) global index of the first species
        """
        if id > 0:
            self.first_spec_index = id
        else:
            msg = [Color.PURPLE, "global index must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    @property
    def last_species_index(self) -> int:
        """Get the global index of the last species."""
        """
        Get the global index (1-based) of the last species of the surface phase.

        Returns
        -------
            global_index: integer
                the (1-based) global index of the last species
        """
        return self.last_spec_index

    @last_species_index.setter
    def last_species_index(self, id: int):
        """Set the global index of the last species."""
        """
        Set the global index (1-based) of the last species of the surface phase.

        Parameters
        ----------
            id: integer
                the (1-based) global index of the last species
        """
        if id > 0:
            self.last_spec_index = id
        else:
            msg = [Color.PURPLE, "global index must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

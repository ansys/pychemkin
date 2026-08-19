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

"""Surface Chemkin Chemistry utilities."""

import copy
import ctypes
from typing import TYPE_CHECKING, Dict, List, Union

import numpy as np
import numpy.typing as npt

from ansys.chemkin.core import chemkin_wrapper as ck_wrapper
from ansys.chemkin.core.color import Color
from ansys.chemkin.core.surface_components import Material
from ansys.chemkin.core.utilities import (
    error_and_exit,
    log_error_message,
    log_info_message,
    log_warning_message,
    warning_and_exit,
)
from ansys.chemkin.core.validation import (
    validate_minimum_value,
    validate_species_array_size,
)

if TYPE_CHECKING:
    from ansys.chemkin.core.chemistry import Chemistry


class Surface:
    """Chemkin surface chemistry module."""

    """
    An extension to the Mixture/Strean and the ReactorModel classes to include
    the surface Material objects associated with the Chemistry Set.
    """

    def __init__(self, chem: "Chemistry"):
        """Create a Surface object."""
        """
        Create a Surface object to store and to provide the surface chemistry
        related information and utilities. The surface properties included in
        this module are the surface temperature, the surface site fractions,
        the bulk species activities, and the bulk species linear growth rates
        (from reactor model solution) of all surface material defined in the
        surface mechanism of the Chemistry Set.
        """
        # chemistry set index
        self._chemset_index = chem.chemid
        # instantiate the gas phase mixture object
        # self.gas = Mixture(chem)
        # number of gas phase species
        self.num_gas_species = chem.kk
        # number of gas phase reactions
        self.num_gas_reactions = chem.ii_gas
        # gas species molar weight [g/mole]
        self.gas_wt = chem.wt
        # verify suerface chemistry
        if not chem.verify_surface_mechanism():
            chem.no_surface_mechanism_declaration()
            exit()
        # get surface chemistry sizes
        self.num_material = chem.number_materials
        self._num_max_phases = chem.max_number_phases
        self._num_max_site_species = chem.max_number_sites
        self._num_max_bulk_species = chem.max_number_bulks
        self.max_total_species = chem.max_total_number_species
        self._num_max_surf_reactions = chem.max_number_surface_reactions
        self.max_total_reactions = chem.max_total_number_reactions
        # get surface material names List[str]
        self.material_names = list(chem.matsymbol)
        # construct material map from information in the chemistry set
        # material name mapping: dict{material name, surface material index}
        # surface material index starts from 1
        self.material_map = dict(chem.material_map)
        # copy all Material objects
        # material objects mapping: dict{material name, Material object}
        self.materials = dict(chem.materials)
        # molecular masses: dict{material name, molar weights}
        self._site_spec_wt: Dict[str, npt.NDArray[np.double]] = {}
        self._bulk_spec_wt: Dict[str, npt.NDArray[np.double]] = {}
        self._site_wt_set: Dict[str, bool] = {}
        self._bulk_wt_set: Dict[str, bool] = {}
        # site fraction of the materials: dict{material_name, site fractions}
        self._site_fractions: Dict[str, npt.NDArray[np.double]] = {}
        self.sitefrac_set: Dict[str, bool] = {}
        # bulk activity of the materials: dict{material_name, bulk activities}
        self._bulk_activities: Dict[str, npt.NDArray[np.double]] = {}
        self.bulkact_set: Dict[str, bool] = {}
        # bulk species growth rates [cm/sec]
        self.bulk_growth_rates: Dict[str, npt.NDArray[np.double]] = {}
        self.bulk_growth_rates_set: Dict[str, bool] = {}
        # site phase density [mole/cm2] of the materials:
        # dict{material_name, site density}
        self.site_density: Dict[str, npt.NDArray[np.double]] = {}
        self.siteden_set: Dict[str, bool] = {}
        # loop over all surface materials
        self._total_site_species = 0
        self._total_bulk_species = 0
        for mname in self.material_names:
            # get the corresponding Material object
            m = self.materials[mname]
            # get sizes
            n_sites = m.num_site_species
            n_bulks = m.num_bulk_species
            self._total_site_species += n_sites
            self._total_bulk_species += n_bulks
            #
            if n_sites > 0:
                # get molecular weights
                self._site_spec_wt[mname] = m.site_wt
                self._site_wt_set[mname] = True
                # surface site coverage
                self._site_fractions[mname] = np.zeros(n_sites, dtype=np.double)
                # surface site density
                site_den = m.get_site_density()
                self.site_density[mname] = copy.deepcopy(site_den)
                self.siteden_set[mname] = True
            if n_bulks > 0:
                # get molecular weights
                self._bulk_spec_wt[mname] = m.bulk_wt
                self._bulk_wt_set[mname] = True
                # bulk species activity
                self._bulk_activities[mname] = np.zeros(n_bulks, dtype=np.double)
                # bulk species linear growth rate
                self.bulk_growth_rates[mname] = np.zeros(n_bulks, dtype=np.double)
                self.bulk_growth_rates_set[mname] = False
        # material surface temperature [K]
        # by default, the surface temperature is the same as
        # the reactor gas-phase temperature
        self._numb_surf_temp_set = 0
        self._surf_temp_set: list[int] = []
        self.material_temperature: list[float] = []

    def check_surface_material(self, matname: str) -> int:
        """Verify the the given material name."""
        """
        Verify the the given material name is assigned to a Material object.

        Parameters
        ----------
            matname: string
                material name/symbol

        Returns
        -------
            mat_id: inetger
                material index, = -1: not found
        """
        return self.require_surface_material(matname)

    def require_surface_material(self, matname: str) -> int:
        """Return the material index or terminate if not found."""
        mat_id = self.material_map.get(matname, -1)
        if mat_id < 0:
            self.material_not_found(matname)
        return mat_id

    def material_not_found(self, matname: str):
        """Return the the material name not found error message."""
        msg = [
            Color.PURPLE,
            "material",
            matname,
            "is not a surface material\n",
            Color.SPACEx6,
            "Available surface materials are",
            str(self.material_names),
            Color.END,
        ]
        error_and_exit(msg)

    @property
    def number_materials(self) -> int:
        """Get the number of surface materials."""
        """
        Get the number of surface materials in the Chemistry Set.

        Return
        ------
            numb_materials: number of surface materials
        """
        return self.num_material

    def get_material_names(self) -> List[str]:
        """Get all surface material symbols."""
        """
        Get all surface material names of the mechnaism.

        Returns
        -------
            material_names: list of strings
                list of surface material names of the mechanism
        """
        if self.number_materials > 0:
            return copy.deepcopy(self.material_names)
        else:
            msg = [
                Color.PURPLE,
                "mechanism contains no surface material",
                Color.END,
            ]
            log_error_message(msg)
            exit()

    @property
    def max_number_sites(self) -> int:
        """Get the max number of site species."""
        """
        Get the max number of site species among all the surface
        materials in the surface mechanism.

        Returns
        -------
            max_site: inetger
                max number of site species dimension
        """
        return self._num_max_site_species

    @property
    def max_number_bulks(self) -> int:
        """Get the max number of bulk species."""
        """
        Get the max number of bulk species among all the surface
        materials in the surface mechanism.

        Returns
        -------
            max_bulk: inetger
                max number of bulk species dimension
        """
        return self._num_max_bulk_species

    @property
    def max_number_total_species(self) -> int:
        """Get the max total number of species."""
        """
        Get the max number of total (gas + site + bulk) species among all the surface
        materials in the surface mechanism.

        Returns
        -------
            max_spec: inetger
                max species dimension
        """
        return self.max_total_species

    @property
    def max_number_surface_reactions(self) -> int:
        """Get the max number of surface reactions."""
        """
        Get the max number of surface reaction among all the surface
        materials in the surface mechanism.

        Returns
        -------
            max_reactions: inetger
                max surface reaction dimension
        """
        return self._num_max_surf_reactions

    @property
    def max_number_total_reactions(self) -> int:
        """Get the max number of total reactions."""
        """
        Get the max number of gas + surface reaction among all the surface
        materials in the surface mechanism.

        Returns
        -------
            max_reactions: inetger
                max total reaction dimension
        """
        return self.max_total_reactions

    @property
    def total_site_species(self) -> int:
        """Get total number of site species from all materials."""
        """
        Get the total number of surface site species from all
        surface materials of the surface mechanism.

        Returns
        -------
            total_sites: integer
                total number of site species
        """
        return self._total_site_species

    @property
    def total_bulk_species(self) -> int:
        """Get total number of bulk species from all materials."""
        """
        Get the total number of bulk species from all
        surface materials of the surface mechanism.

        Returns
        -------
            total_bulks: integer
                total number of bulk species
        """
        return self._total_bulk_species

    def get_site_frac(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get surface site species fractions."""
        """
        Get surface site species fractions.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol

        Returns
        -------
            frac: 1-D double array, dimension = number of site species
                site species fractions
        """
        # verify the material name
        self.require_surface_material(mat_name)

        n_sites = self.number_site_species(mat_name)
        if n_sites <= 0:
            msg = [
                Color.MAGENTA,
                "there is no site species on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            log_warning_message(msg)
            return np.array([0.0], dtype=np.double)
        #
        if self.sitefrac_set.get(mat_name, False):
            frac = copy.deepcopy(self._site_fractions[mat_name])
        else:
            msg = [
                Color.PURPLE,
                "site coverage is not set on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            log_info_message(msg)
            frac = np.full(n_sites, 1.0 / n_sites)
        return frac

    def set_site_frac(
        self,
        mat_name: str,
        recipe: Union[List[tuple[str, float]], npt.NDArray[np.double]],
    ):
        """Set the surface site fractions."""
        """
        Set the fractions of all surface site species of the material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol
            recipe: list of tuples, [(species_symbol, fraction), ... ]
                non-zero mixture composition corresponding to
                the given site fraction array
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        if self.number_site_species(mat_name) <= 0:
            msg = [
                Color.MAGENTA,
                "there is no site species on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            warning_and_exit(msg)
        # get material object
        m = self.materials[mat_name]
        if len(recipe) == 0:
            msg = [Color.PURPLE, "site fraction recipe cannot be empty.", Color.END]
            error_and_exit(msg)
        # reset site fraction values to 0
        self._site_fractions[mat_name][:] = 0.0e0
        # get site species symbols
        specieslist = m.site_symbols
        #
        if isinstance(recipe[0], tuple):
            for sp, x in recipe:
                if sp in specieslist:
                    index = specieslist.index(sp)
                else:
                    msg = [
                        Color.PURPLE,
                        sp.rstrip(),
                        "is not a valid site species of material",
                        mat_name.rstrip(),
                        Color.END,
                    ]
                    error_and_exit(msg)
                validate_minimum_value(
                    value=x,
                    minimum=0.0,
                    message="negative site fraction.",
                )
                # set site fraction
                self._site_fractions[mat_name][index] = x
            self.sitefrac_set[mat_name] = True
        elif isinstance(recipe[0], (float, np.double)):
            ksites = len(recipe)
            validate_species_array_size(
                expected=m.num_site_species,
                actual=ksites,
                context="site fraction",
            )
            if ksites == m.num_site_species:
                for k in range(ksites):
                    self._site_fractions[mat_name][k] = max(recipe[k], 0.0e0)
                self.sitefrac_set[mat_name] = True
        else:
            msg = [
                Color.PURPLE,
                "the argument must be:\n",
                Color.SPACEx6,
                "(1) a list of tuples: [('O2', 0.21), ('N2', 0.79)]\n",
                "or\n",
                Color.SPACEx6,
                "(2) a site fraction array of size = <number of site species>",
                Color.END,
            ]
            error_and_exit(msg)

    def check_site_frac_set(self, mat_name: str) -> bool:
        """Check if the surface site species fractions are provided."""
        """
        Check if the surface site species fractions of a material are assigned.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol

        Returns
        -------
            set: bool
                the status of the site fraction assignment
        """
        return self.sitefrac_set.get(mat_name, False)

    def get_all_site_frac(self) -> npt.NDArray[np.double]:
        """Get site fractions of all materials in the mechanism."""
        """
        Get all site species fractions of all surface materials
        in the mechanism into a 1-D array.

        Returns
        -------
            site_frac: 1-D double array, dimension = total number of site species
                site species fractions
        """
        # get total number of site species from all materials
        total_sites = self.total_site_species
        # some site species exist
        if total_sites > 0:
            site_frac = np.zeros(total_sites, dtype=np.double)
            m_list = self.get_material_names()
            s_count = 0
            # loop over all surface materials
            for mname in m_list:
                # number of site species of the material
                n_sites = self.number_site_species(mname)
                if self.check_site_frac_set(mname):
                    # site fractions available
                    s_frac = self.get_site_frac(mname)
                    for k in range(n_sites):
                        j = s_count + k
                        site_frac[j] = s_frac[k]
                s_count += n_sites
        else:
            site_frac = np.zeros(1, dtype=np.double)
            msg = [
                Color.YELLOW,
                "no site species declared in the surface mechanism.",
                Color.END,
            ]
            log_warning_message(msg)
        #
        return site_frac

    def get_bulk_frac(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get surface bulk species activities/moles."""
        """
        Get surface bulk species activities/moles of the material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol

        Returns
        -------
            frac: 1-D double array, dimension = number of bulk species
                bulk activities/moles
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        n_bulks = self.number_bulk_species(mat_name)
        if n_bulks <= 0:
            msg = [
                Color.MAGENTA,
                "there is no bulk species on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            log_warning_message(msg)
            return np.array([0.0], dtype=np.double)
        #
        if self.bulkact_set.get(mat_name, False):
            frac = copy.deepcopy(self._bulk_activities[mat_name])
        else:
            msg = [
                Color.PURPLE,
                "bulk activities are not set on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            log_info_message(msg)
            frac = np.full(n_bulks, 1.0)
        return frac

    def set_bulk_frac(
        self,
        mat_name: str,
        recipe: Union[List[tuple[str, float]], npt.NDArray[np.double]],
    ):
        """Assign the bulk species activities."""
        """
        Set the bulk species activities of the given material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol
            recipe: list of tuples, [(species_symbol, fraction), ... ]
                non-zero mixture composition corresponding to
                the given site fraction array

        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        if self.number_bulk_species(mat_name) <= 0:
            msg = [
                Color.MAGENTA,
                "there is no bulk species on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            warning_and_exit(msg)
        # get material object
        m = self.materials[mat_name]
        n_bulks = m.num_bulk_species
        if len(recipe) == 0:
            msg = [Color.PURPLE, "bulk activity recipe cannot be empty.", Color.END]
            error_and_exit(msg)
        # reset bulk activities values to 0
        self._bulk_activities[mat_name][:] = 0.0e0
        # get bulk species symbols
        specieslist = m.bulk_symbols
        #
        if isinstance(recipe[0], tuple):
            for sp, x in recipe:
                if sp in specieslist:
                    index = specieslist.index(sp)
                else:
                    msg = [
                        Color.PURPLE,
                        sp.rstrip(),
                        "is not a valid bulk species of material",
                        mat_name.rstrip(),
                        Color.END,
                    ]
                    error_and_exit(msg)
                validate_minimum_value(
                    value=x,
                    minimum=0.0,
                    message="negative bulk activity.",
                )
                # set bulk activity
                self._bulk_activities[mat_name][index] = x
            self.bulkact_set[mat_name] = True
        elif isinstance(recipe[0], (float, np.double)):
            kbulks = len(recipe)
            validate_species_array_size(
                expected=n_bulks,
                actual=kbulks,
                context="bulk activity",
            )
            if kbulks == n_bulks:
                for k in range(kbulks):
                    self._bulk_activities[mat_name][k] = max(recipe[k], 0.0e0)
                self.bulkact_set[mat_name] = True
        else:
            msg = [
                Color.PURPLE,
                "the argument must be:\n",
                Color.SPACEx6,
                "(1) a list of tuples: [('O2', 0.21), ('N2', 0.79)]\n",
                "or\n",
                Color.SPACEx6,
                "(2) a site fraction array of size = <number of site species>",
                Color.END,
            ]
            error_and_exit(msg)

    def check_bulk_act_set(self, mat_name: str) -> bool:
        """Check if the bulk species activities are provided."""
        """
        Check if the bulk species activities of a material are assigned.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol

        Returns
        -------
            set: bool
                the status of the bulk activity assignment
        """
        return self.bulkact_set.get(mat_name, False)

    def get_bulk_growth_rates(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get bulk species linear growth rates."""
        """
        Get surface bulk species linear growth rates of the material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol

        Returns
        -------
            rates: 1-D double array, dimension = number of bulk species
                bulk species linear growth rates [cm/sec]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        n_bulks = self.number_bulk_species(mat_name)
        if n_bulks <= 0:
            msg = [
                Color.MAGENTA,
                "there is no bulk species on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            log_warning_message(msg)
            return np.array([0.0], dtype=np.double)
        #
        if self.bulk_growth_rates_set.get(mat_name, False):
            rates = copy.deepcopy(self.bulk_growth_rates[mat_name])
        else:
            msg = [
                Color.PURPLE,
                "bulk growth rates are not set on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            log_info_message(msg)
            rates = np.full(n_bulks, 0.0)
        return rates

    def set_bulk_growth_rates(
        self,
        mat_name: str,
        rates: npt.NDArray[np.double],
    ):
        """Assign bulk species growth rate."""
        """
        Assign the values of the bulk species growth rate of a surface material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol
            rates: 1-D double array, dimension = numb of bulk species
                bulk species linear growth rates [cm/sec]

        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        n_bulks = self.number_bulk_species(mat_name)
        if n_bulks <= 0:
            msg = [
                Color.MAGENTA,
                "there is no bulk species on surface material",
                mat_name.rstrip(),
                Color.END,
            ]
            log_warning_message(msg)
            exit()
        #
        validate_species_array_size(
            expected=n_bulks,
            actual=len(rates),
            context="bulk growth rate",
        )
        for k in range(n_bulks):
            self.bulk_growth_rates[mat_name][k] = rates[k]
        self.bulk_growth_rates_set[mat_name] = True

    def check_bulk_growth_rate_set(self, mat_name: str) -> bool:
        """Check if the bulk species growth rates are provided."""
        """
        Check if the bulk species linear growth rates of a material are obtained.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol

        Returns
        -------
            set: bool
                the status of the bulk growth rate data
        """
        return self.bulk_growth_rates_set.get(mat_name, False)

    def get_all_bulk_growth_rates(self) -> npt.NDArray[np.double]:
        """Get bulk growth rates of all materials in the mechanism."""
        """
        Get all bulk species linear growth rates of all surface materials
        in the mechanism into a 1-D array.

        Returns
        -------
            bulk_growth_rates: 1-D double array,
                dimension=total number of bulk species
                bulk species growth rates [cm/sec]
        """
        # get total number of bulk species from all materials
        total_bulks = self.total_bulk_species
        # some bulk species exist
        if total_bulks > 0:
            bulk_growth_rates = np.zeros(total_bulks, dtype=np.double)
            m_list = self.get_material_names()
            b_count = 0
            # loop over all surface materials
            for mname in m_list:
                # number of bulk species of the material
                n_bulks = self.number_bulk_species(mname)
                if self.bulk_growth_rates_set.get(mname, False):
                    # bulk growth rates available
                    b_growth = self.get_bulk_growth_rates(mname)
                    for k in range(n_bulks):
                        j = b_count + k
                        bulk_growth_rates[j] = b_growth[k]
                b_count += n_bulks
        else:
            bulk_growth_rates = np.zeros(1, dtype=np.double)
            msg = [
                Color.YELLOW,
                "no bulk species declared in the surface mechanism.",
                Color.END,
            ]
            log_warning_message(msg)
        #
        return bulk_growth_rates

    def get_activity_array(
        self, mat_name: str, molfrac: npt.NDArray[np.double]
    ) -> npt.NDArray[np.double]:
        """Get the species activity array."""
        """
        Get the species activity array (gas, site, and bulk) of
        the named surface material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol
            molfrac: 1-D double array, dimension = number of gas species
                gas species mole fractions next to the surface material

        Returns
        -------
            act: 1-D double array, dimension = total number of species of the material
                activity array (gas + sites + bulks)
        """
        # verify the material name
        self.require_surface_material(mat_name)
        # get the corresponding Material object
        m = self.materials[mat_name]
        # get dimensions
        numgas = self.num_gas_species
        numsites = m.num_site_species
        numbulks = m.num_bulk_species
        # set up site fraction array
        if numsites > 0:
            if self.sitefrac_set.get(mat_name, False):
                sitefrac = self.get_site_frac(mat_name)
            else:
                sitefrac = np.full(numsites, 0.0e0)
        else:
            sitefrac = np.full(1, 1.0e0)
        # set up bulk activities/moles array
        if numbulks > 0:
            if self.bulkact_set.get(mat_name, False):
                bulkact = self.get_bulk_frac(mat_name)
            else:
                bulkact = np.full(numbulks, 1.0e0)
        else:
            bulkact = np.full(1, 1.0e0)
        #
        act = Surface.set_activity_array(
            ngas=numgas,
            nsites=numsites,
            nbulks=numbulks,
            molfrac=molfrac,
            site_frac=sitefrac,
            bulk_act=bulkact,
        )
        return act

    def list_surface_coverage(
        self, mat_name: str, option: str = " ", bound: float = 0.0e0
    ):
        """List the surface site fractions."""
        """
        List the surface site fractions of the given material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol
            option: string, {'all, ' '}, default = 'all'
                flag indicates to list 'all' site species or
                just the site species with non-zero fraction
            bound: double
                minimum fraction value for the site species to be printed
        """
        #
        mat_index = self.check_surface_material(mat_name)
        if mat_index < 0:
            exit()
        #
        site_names = self.surface_species_symbols(mat_name)
        site_frac = self.get_site_frac(mat_name)
        if len(site_names) == 0:
            msg = [
                Color.MAGENTA,
                "material",
                mat_name.rstrip(),
                "has no site species.",
                Color.END,
            ]
            log_warning_message(msg)
        else:
            print(f'listing surface coverage on material "{mat_name}"\n')
            #
            if option.lower() == "all":
                for k, s in enumerate(site_frac):
                    print(f"{site_names[k]:18} :  {s:e}")
            else:
                # list non-zero components
                for k, s in enumerate(site_frac):
                    if s > np.max([bound, 0.0e0]):
                        print(f"{site_names[k]:18} :  {s:e}")
            #
            del site_names, site_frac

    def list_bulk_activity(
        self, mat_name: str, option: str = " ", bound: float = 0.0e0
    ):
        """List the bulk species activity."""
        """
        List the bulk species activity of the given material.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol
            option: string, {'all, ' '}, default = 'all'
                flag indicates to list 'all' bulk species or
                just the bulk species with non-zero fraction
            bound: double
                minimum fraction value for the bulk species to be printed
        """
        #
        mat_index = self.check_surface_material(mat_name)
        if mat_index < 0:
            exit()
        #
        bulk_names = self.get_bulk_symbols(mat_name)
        bulk_frac = self.get_bulk_frac(mat_name)
        if len(bulk_names) == 0:
            msg = [
                Color.MAGENTA,
                "material",
                mat_name.rstrip(),
                "has no bulk species.",
                Color.END,
            ]
            log_warning_message(msg)
        else:
            print(f'listing surface coverage on material "{mat_name}"\n')
            #
            if option.lower() == "all":
                for k, b in enumerate(bulk_frac):
                    print(f"{bulk_names[k]:18} :  {b:e}")
            else:
                # list non-zero components
                for k, b in enumerate(bulk_frac):
                    if b > np.max([bound, 0.0e0]):
                        print(f"{bulk_names[k]:18} :  {b:e}")
            #
            del bulk_names, bulk_frac

    def get_surface_material(self, mat_name: str) -> Material:
        """Get the Material object."""
        """
        Get the Material object associated to the given material name.

        Returns
        -------
            surface_material: Material object
                Material object with the given material name
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials.get(mat_name, None)
        if m is None:
            self.material_not_found(mat_name)
        else:
            return copy.deepcopy(m)

    def number_surface_reactions(self, mat_name: str) -> int:
        """Get the number of surface reactions."""
        """
        Get the number of surface reactions of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            ii_surf: integer
                number of surface reactions
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return 0
        #
        m = self.materials[mat_name]
        ii_surf = m.number_surface_reactions
        return ii_surf

    def number_surface_species(self, mat_name: str) -> int:
        """Get the total number of surface species."""
        """
        Get the total number of surface species (site + bulk) of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            kk_surf: integer
                number of total surface species
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return 0
        #
        m = self.materials[mat_name]
        kk_surf = m.number_site_species + self.number_bulk_species
        return kk_surf

    def total_number_species(self, mat_name: str) -> int:
        """Get the total number of species."""
        """
        Get the total number of species (gas + site + bulk) of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            total_species: integer
                total number of species
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return self.num_gas_species
        #
        m = self.materials[mat_name]
        total_species = m.number_total_species
        return total_species

    def number_phase(self, mat_name: str) -> int:
        """Get the number of phases."""
        """
        Get the number of phases (gas + surface) of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            n_phases: integer
                number of phases
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return 1
        #
        m = self.materials[mat_name]
        n_phases = m.number_phases
        return n_phases

    def first_site_phase_index(self, mat_name: str) -> int:
        """Get the index of the first site phase."""
        """
        Get the index of the first site phase of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            index: integer
                index of the first site phase
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        index = m.first_site_phase_index
        return index

    def last_site_phase_index(self, mat_name: str) -> int:
        """Get the index of the last site phase."""
        """
        Get the index of the last site phase of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            index: integer
                index of the last site phase
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        index = m.last_site_phase_index
        return index

    def first_bulk_phase_index(self, mat_name: str) -> int:
        """Get the index of the first bulk phase."""
        """
        Get the index of the first bulk phase of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            index: integer
                index of the first bulk phase
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        index = m.first_bulk_phase_index
        return index

    def last_bulk_phase_index(self, mat_name: str) -> int:
        """Get the index of the last bulk phase."""
        """
        Get the index of the last bulk phase of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            index: integer
                index of the last bulk phase
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        index = m.last_bulk_phase_index
        return index

    def get_phase_names(self, mat_name: str) -> List[str]:
        """Get all phase names."""
        """
        Get all phase names of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            psymbols: list of strings
                phase names of the surface material
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        psymbols = m.phase_names
        return psymbols

    def number_site_species(self, mat_name: str) -> int:
        """Get the total number of site species."""
        """
        Get the total number of site species of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            n_sites: integer
                total number of site species of the material
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return 0
        #
        m = self.materials[mat_name]
        n_sites = m.num_site_species
        return n_sites

    def number_bulk_species(self, mat_name: str) -> int:
        """Get the total number of bulk species."""
        """
        Get the total number of bulk species of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            n_bulks: integer
                total number of bulk species of the material
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return 0
        #
        m = self.materials[mat_name]
        n_bulks = m.number_bulk_species
        return n_bulks

    def first_site_spec_index(self, mat_name: str) -> int:
        """Get the global index of the first site species."""
        """
        Get the global index of the first site species of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            global_index: integer
                global index of the first site species of the material
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        global_index = m.first_site_species_index
        return global_index

    def last_site_spec_index(self, mat_name: str) -> int:
        """Get the global index of the last site species."""
        """
        Get the global index of the last site species of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            global_index: integer
                global index of the last site species of the material
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        global_index = m.last_site_species_index
        return global_index

    def first_bulk_spec_index(self, mat_name: str) -> int:
        """Get the global index of the first bulk species."""
        """
        Get the global index of the first bulk species of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            global_index: integer
                global index of the first bulk species of the material
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        global_index = m.first_bulk_species_index
        return global_index

    def last_bulk_spec_index(self, mat_name: str) -> int:
        """Get the global index of the last bulk species."""
        """
        Get the global index of the last bulk species of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            global_index: integer
                global index of the last bulk species of the material
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        if mat_id < 0:
            return -1
        #
        m = self.materials[mat_name]
        global_index = m.last_bulk_species_index
        return global_index

    def get_site_symbols(self, mat_name: str) -> List[str]:
        """Get all site species symbols."""
        """
        Get all site species symbols of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            ssymbols: list of string
                all site species symbols of the material
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        ssymbols = m.site_species_names
        return ssymbols

    def get_bulk_symbols(self, mat_name: str) -> List[str]:
        """Get all bulk species symbols."""
        """
        Get all bulk species symbols of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            bsymbols: list of string
                all bulk species symbols of the material
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        bsymbols = m.bulk_species_names
        return bsymbols

    def surface_species_symbols(self, mat_name: str) -> List[str]:
        """Get all surface species symbols of the material."""
        """
        Get all surface species symbols of the material.

        Parameters
        ----------
            mat_name: string
            material name/symbol

        Returns
        -------
            surf_symbols: list of string
                all surface species symbols of the material

        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        ssymbols = m.site_species_names
        bsymbols = m.bulk_species_names
        surf_symbols = ssymbols + bsymbols
        return surf_symbols

    def get_surf_specindex(self, specname: str) -> tuple[str, int, int]:
        """Get the local and the global species indices."""
        """
        Get the local/surface material and the global species indices
        of the given surface species.

        Parameters
        ----------
            specname: string
                surface species symbol

        Returns
        -------
            name: string
                the material name where the surface species belongs
            global_index: integer
                the global index of the surface species
            local_index: integer
                (0-based) species index of the phase
        """
        global_index = -1
        local_index = -1
        n = 0
        while global_index < 0 and n < self.num_material:
            # material name
            name = self.material_names[n]
            # get the Material object
            m = self.materials[name]
            global_index, local_index, _ = m.get_surf_specindex(specname, mode="silent")
            n += 1
        if local_index >= 0:
            # print(f"species {specname} found on material {name}")
            # print(f"species type = {phase_type}")
            # print(f"global index = {global_index}, {phase_type}"
            #      f"index = {local_index}")
            pass
        else:
            name = "<not found>"
            global_index = -1
            print(f"species {specname} is not a surface species.")
        #
        return name, global_index, local_index

    def site_phase_density(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the site density of the site phases."""
        """
        Get the site density of the site phases of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            site_den: 1-D double array, dimension = number of site phases
                site density [mole/cm2]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        site_den = m.get_site_density()
        # the first site phase
        pstart = m.first_site_phase_index
        pstop = pstart + m.num_site_phase
        return site_den[pstart:pstop]

    def get_site_species_occupany(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the occupancy."""
        """
        Get the occupancy of the surface site species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            occupancy: 1-D double array, dimension = number of total site species
                site occupancy [-]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        occupancy = m.get_site_occupancy()
        return occupancy

    def get_site_wt(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the molar masses of the site species."""
        """
        Get the molar masses of the site species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            site_wt: 1-D double array, dimension = number of total site species
                molecular weight [g/mole]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        if self._site_wt_set.get(mat_name, False):
            site_wt = self._site_spec_wt[mat_name]
        else:
            m = self.materials[mat_name]
            site_wt = m.get_site_molar_weights()
            self._site_spec_wt[mat_name] = site_wt
            self._site_wt_set[mat_name] = True
        return site_wt

    def get_bulk_wt(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the molar masses of the bulk species."""
        """
        Get the molar masses of the bulk species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            bulk_wt: 1-D double array, dimension = number of total bulk species
                molecular weight [g/mole]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        if self._bulk_wt_set.get(mat_name, False):
            bulk_wt = self._bulk_spec_wt[mat_name]
        else:
            m = self.materials[mat_name]
            bulk_wt = m.get_bulk_molar_weights()
            self._bulk_spec_wt[mat_name] = bulk_wt
            self._bulk_wt_set[mat_name] = True
        return bulk_wt

    def get_surf_wt(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the molar masses of all surface species."""
        """
        Get the molar masses of all surface species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            surf_wt: 1-D double array, dimension = number of total surface species
                molecular weight [g/mole]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        site_wt = self.get_site_wt(mat_name)
        bulk_wt = self.get_bulk_wt(mat_name)
        surf_wt = np.hstack((site_wt, bulk_wt))
        return surf_wt

    def get_site_species_h(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the site species enthalpy."""
        """
        Get the enthalpy of the site species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            site_h: 1-D double array, dimension = number of total site species
                enthalpy [ergs/mole]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        site_h = m.get_site_species_h()
        return site_h

    def get_bulk_species_h(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the bulk species enthalpy."""
        """
        Get the enthalpy of the bulk species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            bulk_h: 1-D double array, dimension = number of total bulk species
                enthalpy [ergs/mole]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        bulk_h = m.get_bulk_species_h()
        return bulk_h

    def get_surface_species_h(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the enthalpy of all surface species."""
        """
        Get the enthalpy of all surface species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            surf_h: 1-D double array, dimension = number of total surface species
                enthalpy [ergs/mole]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        bulk_h = m.get_bulk_species_h()
        site_h = m.get_site_species_h()
        surf_h = np.hstack((site_h, bulk_h))
        return surf_h

    def get_bulk_species_cp(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the bulk species specific heat capacity."""
        """
        Get the specific heat capacity of the bulk species of the material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            bulk_cp: 1-D double array, dimension = number of total bulk species
                specific heat capacity [ergs/mole-K]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        bulk_cp = m.get_bulk_species_cp()
        return bulk_cp

    def get_bulk_species_density(self, mat_name: str) -> npt.NDArray[np.double]:
        """Get the bulk species density."""
        """
        Get the density of the bulk species of the material

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            bulk_den: 1-D double array, dimension = number of total bulk species
                mass density [g/cm3]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        bulk_den = m.get_bulk_density()
        return bulk_den

    @staticmethod
    def set_activity_array(
        ngas: int,
        nsites: int,
        nbulks: int,
        molfrac: npt.NDArray[np.double],
        site_frac: npt.NDArray[np.double],
        bulk_act: npt.NDArray[np.double],
        min_value: float = 0.0e0,
    ) -> npt.NDArray[np.double]:
        """Set up the activity array."""
        """
        Set up the activity array from the given gas mole fractions,
        the surface site fractions, the surface bulk activities/moles.
        It is used to call the reaction rate API.

        Parameters
        ----------
            ngas: integer
                number of gas species
            nsites: integer
                number of all site species of the material
            nbulks: integer
                number of all bulk species of the material
            molfrac: 1-D double array, dimension = number_gas_species
                gas species mole fractions
            site_frac: 1-D double array, dimension = number_site_species
                surface site species fraction [-]
            bulk_act: 1-D double array, dimension = number_bulk_species
                bulk species activity/mole fraction [-]
            min_value: double, optional
                species fraction value threshold

        Returns
        -------
            activity: 1-D double array, dimension = total_species
                the full activity arrays used in the chemkin-VGD-API calls
        """
        # check length
        if ngas != len(molfrac):
            msg = [
                Color.PURPLE,
                "gas mole fraction array size does not match",
                "the number of gas species",
                str(ngas),
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # check site species fractions
        if nsites > 0:
            sfrac = copy.deepcopy(site_frac)
        else:
            sfrac = np.full(1, 1.0e0)
        # check bulk species activities
        if nbulks > 0:
            bact = copy.deepcopy(bulk_act)
        else:
            bact = np.full(1, 1.0e0)
        #
        species_sum = ngas + nsites + nbulks
        # initialize the activity array
        activity = np.zeros(species_sum, dtype=np.double)
        count = 0
        x_sum = 0.0e0
        # set up gas species portion of the activity array
        for k, f in enumerate(molfrac):
            if f > min_value:
                activity[k] = f
            else:
                activity[k] = 0.0e0
            x_sum += activity[k]
            count += 1
        # normalization
        if x_sum > 0.0e0:
            for k in range(0, count):
                activity[k] /= x_sum
        else:
            msg = [
                Color.PURPLE,
                "sum of gas species mole fractions <= 0:",
                str(x_sum),
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # set up surface site species portion of the activity array
        x_sum = 0.0e0
        for k, f in enumerate(sfrac):
            if f >= min_value:
                sfrac[k] = f
            else:
                sfrac[k] = 0.0e0
            x_sum += sfrac[k]
        #
        if x_sum > 0.0e0:
            for f in sfrac:
                activity[count] = f / x_sum
                count += 1
        else:
            # site fractions are not given. set them to the same fraction
            f_eq = 1.0e0 / nsites
            for f in sfrac:
                activity[count] = f_eq
                count += 1
        # set up surface bulk species portion of the activity array
        x_sum = 0.0e0
        for k, f in enumerate(bact):
            bact[k] = max(f, 0.0e0)
            x_sum += bact[k]
        if x_sum > 0.0e0:
            # do not normalize bulk activities
            for f in bact:
                activity[count] = f
                count += 1
        else:
            # use default bulk activities = 1.0
            for f in bact:
                activity[count] = 1.0e0
                count += 1
        #
        del sfrac, bact
        #
        return activity

    def check_material_temperature(self, mat_name: str) -> int:
        """Verify the surface temperature is given."""
        """
        Verify the surface temperature of this material is explicitly
        provided.

        Returns
        -------
            status: integer
                indication of the material temperature specification
                = 0: not provided; > 0:  provided
        """
        status = -1
        if self._numb_surf_temp_set > 0:
            index = self.check_surface_material(mat_name)
            if index >= 0:
                try:
                    status = self._surf_temp_set.index(index)
                except ValueError:
                    status = -1
        return status + 1

    def get_material_temperature(self, mat_name: str) -> float:
        """Get the surface temperature."""
        """
        Get the surface temperature of the given material.

        Parameters
        ----------
            mat_name: string
                surface material name

        Returns
        -------
            temp: double
                surface temperature [K] of the given material
        """
        index = self.check_surface_material(mat_name)
        if index >= 0:
            try:
                loc = self._surf_temp_set.index(index)
                return self.material_temperature[loc]
            except ValueError:
                msg = [
                    Color.YELLOW,
                    "temperature of material",
                    mat_name.rstrip(),
                    "has not been assigned",
                    Color.END,
                ]
                log_info_message(msg)
                return -1.0
        else:
            return 0.0

    def set_material_temperature(self, mat_name: str, temp: float):
        """Set the value of the surface temperature."""
        """
        Set the value of the surface temperature of a material.

        Parameters
        ----------
            mat_name: string
                surface material name
            temp: double
                surface temperature [K] of the given material
        """
        if temp < 100.0:
            # invalid temperature value
            msg = [Color.PURPLE, "material temperature must >= 100. [K]", Color.END]
            error_and_exit(msg)
        else:
            index = self.check_surface_material(mat_name)
            if index >= 0:
                if index not in self._surf_temp_set:
                    self._numb_surf_temp_set += 1
                    self._surf_temp_set.append(index)
                    self.material_temperature.append(temp)
                else:
                    loc = self._surf_temp_set.index(index)
                    self.material_temperature[loc] = temp

    @staticmethod
    def surface_rate_of_production(
        chem_id: int,
        mat_id: int,
        numb_gas: int,
        numb_phase: int,
        numb_sites: int,
        numb_bulks: int,
        p: float,
        t: float,
        molfrac: npt.NDArray[np.double],
        site_frac: npt.NDArray[np.double],
        bulk_act: npt.NDArray[np.double],
        site_den: npt.NDArray[np.double],
        mode: str,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get species molar rate of production."""
        """
        Get species molar rate of production from the given gas mixture
        and surface condition: gas pressure, gas composition, surface temperature,
        and surface material coverage.

        Parameters
        ----------
            chem_id: integer
                chemistry set index associated with the mixture
            mat_id: integer
                surface material index
            numb_gas: integer
                number of gas species
            numb_phase: integer
                number of phases (including the gas pahse)
                associated with the chemistry set
            numb_sites: integer
                number of surface site species of the surface material
            numb_bulks: integer
                number of bulk species of the surface material
            p: double
                mixture pressure in [dynes/cm2]
            t: double
                surface temperature in [K]
            molfrac: 1-D double array, dimension = number_gas_species
                mixture composition given as mole fractions
            site_frac: 1-D double array, dimension = number_site_species
                surface site species fraction [-]
            bulk_act: 1-D double array, dimension = number_bulk_species
                bulk species activity/mole fraction [-]
            site_den: 1-D double array, dimension = number_phases
                surface phase density [mole/cm2]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            rop: 1-D double array, dimension = number_species
                species molar rate of production in [mol/cm2-sec]
            phase_prod_rates: 1-D double array, dimension = number of phases
                surface site phase production rates in [mole/cm2-sec]
        """
        # check inputs
        if chem_id < 0:
            msg = [Color.PURPLE, "invalid chemistry.", Color.END]
            log_error_message(msg)
            exit()
        if p <= 0.0 or (p * t) <= 0.0:
            msg = [
                Color.PURPLE,
                "invalid pressure and/or temperature value(s).",
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # number species
        kgas = len(molfrac)
        if kgas != numb_gas:
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # initialization
        # species rate of production by surface reactions [mole/cm2-sec]
        numb_species = numb_gas + numb_sites + numb_bulks
        rop = np.empty(numb_species, dtype=np.double)
        # surface site phase production rate by surface reactions [mole/cm2-sec]
        phase_prod_rates = np.empty(numb_phase, dtype=np.double)
        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chem_id)
        mat_index = ctypes.c_int(mat_id)
        pp = ctypes.c_double(p)  # pressure scalar
        tt = ctypes.c_double(t)  # temperature scalar
        # construct the activity array
        # (gas mole fraction, site fraction, bulk activity)
        act = Surface.set_activity_array(
            numb_gas,
            numb_sites,
            numb_bulks,
            molfrac,
            site_frac,
            bulk_act,
        )
        # compute mass density from mass fraction
        ierr = ck_wrapper.chemkin.KINGetSurfaceProductionRates(
            ctypes.byref(chemset_index),
            ctypes.byref(mat_index),
            ctypes.byref(pp),
            ctypes.byref(tt),
            act,
            site_den,
            rop,
            phase_prod_rates,
        )
        if ierr == 0:
            return rop, phase_prod_rates
        else:
            # failed to compute species molar rates of production
            msg = [
                Color.PURPLE,
                "failed to compute species molar rates of production.",
                Color.END,
            ]
            log_error_message(msg)
            exit()

    @staticmethod
    def surface_reaction_rates(
        chem_id: int,
        mat_id: int,
        numb_gas: int,
        numb_sites: int,
        numb_bulks: int,
        numb_reaction: int,
        p: float,
        t: float,
        molfrac: npt.NDArray[np.double],
        site_frac: npt.NDArray[np.double],
        bulk_act: npt.NDArray[np.double],
        site_den: npt.NDArray[np.double],
        mode: str,
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get molar reaction rates of the reactive surface."""
        """
        Get molar reaction rates of the reactive surface and the gas mixture next to it
        from the given condition: pressure, surface temperature,
        gas species composition, and surface species coverage.

        Parameters
        ----------
            chem_id: integer
                chemistry set index
            mat_id: integer
                surface material index
            numb_gas: integer
                number of gas species
            numb_sites: integer
                number of surface site species of the surface material
            numb_bulks: integer
                number of bulk species of the surface material
            numb_reaction: integer
                number of gas reactions associated with the chemistry set
            p: double
                mixture pressure in [dynes/cm2]
            t: double
                surface temperature in [K]
            molfrac: 1-D double array, dimension = number_gas_species
                mixture composition given as mole fractions
            site_frac: 1-D double array, dimension = number_site_species
                surface site species fraction [-]
            bulk_act: 1-D double array, dimension = number_bulk_species
                bulk species activity/mole fraction [-]
            site_den: 1-D double array, dimension = number_phases
                surface phase density [mole/cm2]
            mode: string, {'mole', 'mass'}
                flag indicates the frac array is 'mass' or 'mole' fractions

        Returns
        -------
            k_forward: 1-D double array, dimension = number of surface reaction
                forward molar rates of the reactions in [mol/cm2-sec]
            k_reverse: 1-D double array, dimension = number of surface reaction
                reverse molar rates of the reactions in [mol/cm2-sec]
        """
        # check inputs
        if chem_id < 0:
            msg = [Color.PURPLE, "invalid chemistry.", Color.END]
            log_error_message(msg)
            exit()
        if p <= 0.0 or (p * t) <= 0.0:
            msg = [
                Color.PURPLE,
                "invalid pressure and/or temperature value(s).",
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # number species
        kgas = len(molfrac)
        if kgas != numb_gas:
            msg = [
                Color.PURPLE,
                mode,
                "fraction and molar mass arrays",
                "must have the same size =",
                str(kgas),
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # initialization
        k_forward = np.empty(numb_reaction, dtype=np.double)
        k_reverse = np.empty_like(k_forward, dtype=np.double)
        # convert parameters to c pointers
        chemset_index = ctypes.c_int(chem_id)
        mat_index = ctypes.c_int(mat_id)
        pp = ctypes.c_double(p)  # pressure scalar
        tt = ctypes.c_double(t)  # temperature scalar
        # construct the activity array
        # (gas mole fraction, site fraction, bulk activity)
        act = Surface.set_activity_array(
            numb_gas,
            numb_sites,
            numb_bulks,
            molfrac,
            site_frac,
            bulk_act,
        )
        # compute mass density from mass fraction
        ierr = ck_wrapper.chemkin.KINGetSurfaceReactionRates(
            chemset_index, mat_index, tt, pp, act, site_den, k_forward, k_reverse
        )
        if ierr == 0:
            return k_forward, k_reverse
        else:
            # failed to compute reaction rates
            msg = [
                Color.PURPLE,
                "failed to compute surface reaction rates,",
                "error code =",
                str(ierr),
                Color.END,
            ]
            log_error_message(msg)
            exit()

    def set_surface_coverage(
        self, mat_name: str
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get the surface coverage."""
        """
        Get the surface coverage (site phase density, site species fractions,
        and bulk activities) of the given material.

        Parameters
        ----------
            mat_name: string
                material name/symbol

        Returns
        -------
            site_den: 1-D double array, dimension = number_phases
                surface phase density [mole/cm2]
            site_frac: 1-D double array, dimension = number_site_species
                surface site species fraction [-]
            bulk_act: 1-D double array, dimension = number_bulk_species
                bulk species activity/mole fraction [-]
        """
        # get the corresponding Material object
        m = self.materials[mat_name]
        # get surface site phase density [mole/cm2]
        if self.siteden_set.get(mat_name, False):
            sden = self.site_density[mat_name]
        else:
            # get the site density from surface data
            sden = m.get_site_density()
            self.site_density[mat_name] = sden
            self.siteden_set[mat_name] = True
        # get surface site fractions
        if self.sitefrac_set.get(mat_name, False):
            sfrac = self.get_site_frac(mat_name)
        else:
            # surface coverage is not given
            sfrac = np.zeros(max(m.number_site_species, 1), dtype=np.double)
        # get bulk activities
        if self.bulkact_set.get(mat_name, False):
            bfrac = self.get_bulk_frac(mat_name)
        else:
            # bulk activity is not given
            bfrac = np.zeros(max(m.number_bulk_species, 1), dtype=np.double)
        return sden, sfrac, bfrac

    def rop_surf(
        self,
        mat_name: str,
        pres: float,
        surf_temp: float,
        molfrac: npt.NDArray[np.double],
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get species molar rate of production."""
        """
        Get species molar rate of production from the given mixture condition:
        pressure, surface temperature, gas mixture composition, and surface coverage.

        Parameters
        ----------
            mat_name: string
                material name/symbol
            pres: double
                gas mixture pressure in [dynes/cm2]
            surf_temp: double
                surface temperature in [K]
            molfrac: 1-D double array, dimension = number_gas_species
                gas mixture composition in mole fractions

        Returns
        -------
            rop_surf: 1-D double array, dimension = number_species
                species molar rate of production in [mol/cm3-sec]
            site_prodrate: 1-D double array, dimension = number of phases
                surface site phase production rates in [mole/cm2-sec]
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        # check temperature
        if surf_temp <= 250.0:
            msg = [Color.PURPLE, "surface temperature must > 250.0 [K]", Color.END]
            log_error_message(msg)
            exit()
        # check pressure
        if pres <= 0.0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] must > 0.0.",
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # get the corresponding Material object
        m = self.materials[mat_name]
        # get sizes
        n_gas = self.num_gas_species
        num_phase = m.number_phases
        n_sites = m.num_site_species
        n_bulks = m.num_bulk_species
        # get surface overages
        sden, sfrac, bfrac = self.set_surface_coverage(mat_name)
        # compute the species production rates by surface reactions of this material
        rop_surf, site_prodrate = Surface.surface_rate_of_production(
            chem_id=self._chemset_index,
            mat_id=mat_id,
            numb_gas=n_gas,
            numb_phase=num_phase,
            numb_sites=n_sites,
            numb_bulks=n_bulks,
            p=pres,
            t=surf_temp,
            molfrac=molfrac,
            site_frac=sfrac,
            bulk_act=bfrac,
            site_den=sden,
            mode="mole",
        )
        #
        return rop_surf, site_prodrate

    def rxnrates_surf(
        self,
        mat_name: str,
        pres: float,
        surf_temp: float,
        molfrac: npt.NDArray[np.double],
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Get molar rates of the surface reactions."""
        """
        Get molar rates of the surface reactions from the given gas and
        surface conditions: pressure, surface temperature, gas mixture composition,
        and surface coverage.

        Parameters
        ----------
            mat_name: string
                material name/symbol
            pres: double
                gas mixture pressure in [dynes/cm2]
            surf_temp: double
                surface temperature in [K]
            molfrac: 1-D double array, dimension = number of gas species
                gas mixture composition in mole fractions

        Returns
        -------
            k_forward: 1-D double array, dimension = number of surface reaction
                forward molar rates of the reactions in [mol/cm2-sec]
            k_reverse: 1-D double array, dimension = number of surface reaction
                reverse molar rates of the reactions in [mol/cm2-sec]
        """
        # verify the material name
        mat_id = self.check_surface_material(mat_name)
        # check temperature
        if surf_temp <= 250.0:
            msg = [Color.PURPLE, "surface temperature must > 250.0 [K]", Color.END]
            log_error_message(msg)
            exit()
        # check pressure
        if pres <= 0.0:
            msg = [
                Color.PURPLE,
                "mixture pressure [dynes/cm2] must > 0.0.",
                Color.END,
            ]
            log_error_message(msg)
            exit()
        # initialization
        # get the material object
        m = self.materials[mat_name]
        # get sizes
        n_gas = self.num_gas_species
        numsurf_reactions = m.number_surface_reactions
        n_sites = m.num_site_species
        n_bulks = m.num_bulk_species
        # get surface overages
        sden, sfrac, bfrac = self.set_surface_coverage(mat_name)
        #
        k_forward = np.empty(numsurf_reactions, dtype=np.double)
        k_reverse = np.empty_like(k_forward, dtype=np.double)
        #
        # mixture mole fraction given
        k_forward, k_reverse = Surface.surface_reaction_rates(
            chem_id=self._chemset_index,
            mat_id=mat_id,
            numb_gas=n_gas,
            numb_sites=n_sites,
            numb_bulks=n_bulks,
            numb_reaction=numsurf_reactions,
            p=pres,
            t=surf_temp,
            molfrac=molfrac,
            site_frac=sfrac,
            bulk_act=bfrac,
            site_den=sden,
            mode="mole",
        )
        return k_forward, k_reverse

    def get_bulk_molar_growth_rates(
        self, mat_name: str, surfrop: npt.NDArray[np.double]
    ) -> npt.NDArray[np.double]:
        """Get the molar growth rates of the bulk species."""
        """
        Get the molar growth rates of the bulk species of the material
        from the surface species production rates. Use 'rop_surf' method to
        obtain the species production rates before calling this method.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol associated with the given
                species production rates
            surfrop: 1-D double array,
                dimension = total species of the material (including gas species)
                species molar production rates by the surface reactions of
                the material [mole/cm2-sec]

        Returns
        -------
            growth_rates: 1-D double array, dimension = number of bulk species
                molar growth rates of the bulk species [mole/cm2-sec]
        """
        # verify the material name
        self.require_surface_material(mat_name)
        #
        m = self.materials[mat_name]
        # check array size
        n_species = self.total_number_species(mat_name)
        if len(surfrop) != n_species:
            msg = [
                Color.PURPLE,
                "the rate of production array 'surf_ROP' has wrong size\n",
                Color.SPACEx6,
                "expected:",
                str(n_species),
                Color.SPACEx6,
                "actual:",
                str(len(surfrop)),
                Color.END,
            ]
            log_error_message(msg)
            exit()
        bulk1 = m.first_bulk_species_index
        nbulks = m.num_bulk_species
        bulkn = bulk1 + nbulks
        # initialization
        growth_rates = np.zeros(max(nbulks, 1), dtype=np.double)
        if nbulks > 0:
            # extract the bulk species molar growth rates from surfrop [mole/cm2-sec]
            i = 0
            for k in range(bulk1, bulkn):
                growth_rates[i] = surfrop[k]
                i += 1
        return growth_rates

    def get_bulk_mass_growth_rates(
        self, mat_name: str, surfrop: npt.NDArray[np.double]
    ) -> npt.NDArray[np.double]:
        """Get the mass growth rates of the bulk species."""
        """
        Get the mass growth rates of the bulk species of the material from
        the surface species production rates. Use 'rop_surf' method to
        obtain the species production rates before calling this method.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol associated with
                the given species production rates
            surf_rop: 1-D double array,
                dimension = total species of the material (including gas species)
                species molar production rates by the surface reactions
                of the material [mole/cm2-sec]

        Returns
        -------
            growth_rates: 1-D double array, dimension = number of bulk species
                mass growth rates of the bulk species [g/cm2-sec]
        """
        # get molar growth rates
        growth_rates = self.get_bulk_molar_growth_rates(mat_name, surfrop)
        # get bulk species molar mass [g/mole]
        m = self.materials[mat_name]
        nbulks = m.num_bulk_species
        bulk_wt = m.get_bulk_molar_weights()
        if nbulks > 0:
            # convert to mass growth rates [g/cm2-sec]
            for k in range(len(growth_rates)):
                growth_rates[k] *= bulk_wt[k]
        return growth_rates

    def get_bulk_linear_growth_rates(
        self, mat_name: str, surfrop: npt.NDArray[np.double]
    ) -> npt.NDArray[np.double]:
        """Get the linear growth rates of the bulk species."""
        """
        Get the linear growth rates of the bulk species of the material
        from the surface species production rates. Use 'rop_surf' method to
        obtain the species production rates before calling this method.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol associated with the
                given species production rates
            surfrop: 1-D double array,
                dimension = total species of the material (including gas species)
                species molar production rates by the surface reactions
                of the material [mole/cm2-sec]

        Returns
        -------
            linear_growth: 1-D double array, dimension = number of bulk species
                linear growth rates of the bulk species [cm/sec]
        """
        # get molar growth rates
        growth_rates = self.get_bulk_molar_growth_rates(mat_name, surfrop)
        # get bulk species molar mass [g/mole]
        m = self.materials[mat_name]
        nbulks = m.num_bulk_species
        bulk_wt = m.get_bulk_molar_weights()
        bulk_den = m.get_bulk_density()
        # initialization
        linear_growth = np.zeros(max(nbulks, 1), dtype=np.double)
        if nbulks > 0:
            # convert to linear growth rates [cm/sec]
            for k, g in enumerate(growth_rates):
                den = bulk_den[k]
                if den > 0.0e0:
                    linear_growth[k] = g * bulk_wt[k] / den
                else:
                    linear_growth[k] = 0.0e0
        #
        return linear_growth

    def get_gas_production_rates(
        self, mat_name: str, surfrop: npt.NDArray[np.double]
    ) -> npt.NDArray[np.double]:
        """Get the molar production rates of the gas species."""
        """
        Get the molar production rates of the gas species from the surface
        species production rates. Use 'rop_surf' method to obtain
        the species production rates before calling this method.

        Parameters
        ----------
            mat_name: string
                surface material name/symbol associated with
                the given species production rates
            surfrop: 1-D double array,
                dimension = total species of the material (including gas species)
                species molar production rates by the surface reactions
                of the material [mole/cm2-sec]

        Returns
        -------
            gas_prod_rates: 1-D double array, dimension = number of gas species
                gas species production rates from the surface reactions
                of the material [mole/cm2-sec]
        """
        # get molar growth rates
        growth_rates = self.get_bulk_molar_growth_rates(mat_name, surfrop)
        #
        n_gas = self.num_gas_species
        # initialization
        gas_prod_rates = np.zeros(n_gas, dtype=np.double)
        # extract gas species production rates from the surface material [mole/cm2-sec]
        for k in range(n_gas):
            gas_prod_rates[k] = growth_rates[k]
        return gas_prod_rates

    def stefan_mass_flux(self, mat_name: str, surfrop: npt.NDArray[np.double]) -> float:
        """Get the Stefan mass flux due to surface reactions."""
        """
        Get the Stefan mass flux due to surface reactions.
        This mass flux represents the net mass flux created by the adsorption
        and the desorption of the gas species due to surface reactions.
        Divide this mass flux by the gas density above the surface to
        obtain the Stefan velocity. Use 'rop_surf' method to obtain
        the species production rates before calling this method.

        Returns
        -------
            mass_flux: double
                the Stefan mass flux due to surface reactions [g/cm2-sec]
        """
        # get the gas species molar production rates
        gas_prod_rates = self.get_gas_production_rates(mat_name, surfrop)
        #
        gas_wt = self.gas_wt
        mass_flux = 0.0e0
        for k, g in enumerate(gas_prod_rates):
            mass_flux += g * gas_wt[k]
        return mass_flux

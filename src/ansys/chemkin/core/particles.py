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

"""Aerosol size distribution models."""

from typing import NoReturn

import numpy as np
import numpy.typing as npt
from scipy.stats import gamma, lognorm

from ansys.chemkin.core.color import Color
from ansys.chemkin.core.constants import AVOGADRO
from ansys.chemkin.core.surface_components import Material
from ansys.chemkin.core.utilities import error_and_exit


class Aerosol(Material):
    """Aerosol material object."""

    def __init__(self, chemset_index: int, material_index: int, label: str = ""):
        """Create an aerosol material object."""
        # The aerosol material is modeled as a dispersed surface material in Chemkin.
        # Nucleation and mass growth of the aerosol are modeled as surface reactions on
        # the dispersed surface material. The physical processes such as coagulation
        # and aggregation are modeled by the usual aerosol models.
        # The size class of the particles is represented by the number of moles of the
        # bulk species that makes up the core of the particles. The core bulk species
        # is defined in the surface mechanism by the auxilary tag "core" to a bulk
        # species. Alternatively, it can also be defined as the only bulk species
        # product of the nucleation reactions. The native open site species on the
        # particle surface is used to balance the nucleation reactions and is defined
        # by the "native" tag in the surface mechanism. By default, the first surface
        # site species in the surface mechanism is used as the native open site species.
        super().__init__(chemset_index, material_index, label)
        self.material_type = "aerosol"
        self._dispersed = True
        # nucleation reaction information
        self.number_nucleation_reactions: int = 0
        self.nucleation_reaction_indices: list[int] = []
        self.inception_size: list[float] = []
        # the bulk species that makes up the core of the particles
        self.core_bulk_species_index: int = 0
        self.core_bulk_species_symbol: str = ""
        # the default open active site species on the particle surface
        self.native_open_site_species_index: int = 0
        self.native_open_site_species_symbol: str = ""
        # particle size distribution tracking method: "moments" or "sectional"
        self.population_distribution_method: str = ""
        # aggregation parameters
        # coagulation parameters


class PDF:
    """Object representing the population density function (PDF) of aerosols."""

    def __init__(self, core_wt: float, bulk_density: float):
        """Initialize the aerosol pupolation density function object."""
        """
        Instantiate a population density function (PDF) object with the phisical
        peoperties of the aerosol particles. The sizes of the particles are represented
        by the number of the core element (an atom or a molecule) consistitutes the
        particles. The core_wt is the atomic/molar mass of the core element in g/cm3.
        The bulk density [g/cm3] indicates how the elements in the particles are packed
        on average. The particle size class is roughly porportional to the volume of
        the aerosol particle (primary partile).
        
        Parameters
        ----------
            core_wt: float, > 0.
                the molecular mass of the atom that constitutes the core of the
                aerosol particles [g/mole].
            bulk_density: float, > 0.
                the bulk density of the primary particles [g/cm3].
        """
        # particle physical properties
        self._avogadro_number = AVOGADRO
        self.core_molar_mass = core_wt  # g
        self.bulk_density = bulk_density  # g/cm3
        self.mass_per_class = self.core_molar_mass / self._avogadro_number  # g
        self.vol_per_class = self.mass_per_class / self.bulk_density * 1.0e12  # micron3
        # population density function (PDF) parameters
        self.m0 = 0.0
        # mean
        self.mean = 0.0
        # variance
        self.variance = 0.0
        # standard deviation
        self.std_dev = 0.0
        # minimum/ lower bond
        self.mean_shifted = 0.0
        self.x_min = 0.0

    @property
    def vol0(self) -> float:
        """Get volume of the core class unit."""
        return self.vol_per_class  # micron3

    @property
    def mean_class(self) -> float:
        """Get the mean size class of the distribution."""
        return self.mean  # class

    @property
    def mean_volume(self) -> float:
        """Get the volume corresponding to the mean size class."""
        return self.mean * self.vol_per_class  # micron3

    def class_to_volume(self, x: float) -> float:
        """Convert size class to particle volume."""
        """
        Convert the given particle size class to volume in micron^3.
        
        Parameters
        ----------
            x: float, > 0.
            particle size class
        
        Returns
        -------
            vol: float, >= 0.
            particle volume [micron3]
        """
        vol = x * self.vol_per_class
        return vol  # micron3

    def volume_to_class(self, v: float) -> float:
        """Convert particle volume to size class."""
        """
        Convert the given particle volume in micron^3 to size class.
        
        Parameters
        ----------
            v: float, > 0.
            particle volume [micron3]

        Returns
        -------
            x: float, >= 0.
            particle size class (number of cores)
        """
        x = v / self.vol_per_class
        return x

    def class_to_diameter(self, x: float) -> float:
        """Convert particle class to particle diameter."""
        """
        Convert the given particle class to particle diameter in micron.
        
        Parameters
        ----------
            x: float, > 0.
            particle size class
        
        Returns
        -------
            diameter: float, >= 0.
            particle diameter [micron]
        """
        class_vol = self.class_to_volume(x)
        diameter = np.cbrt(6.0 * class_vol / np.pi)  # cm
        return diameter  # micron

    def diameter_to_class(self, diameter: float) -> float:
        """Convert particle diameter to particle class."""
        """
        Convert the given particle diameter in micron to size class.
        
        Parameters
        ----------
            diameter: float, > 0.
            particle diameter [micron]
        
        Returns
        -------
            x: float, >= 0.
            particle size class (number of cores)
        """
        d = diameter  # micron
        class_vol = np.pi * d * d * d / 6.0
        x = self.volume_to_class(class_vol)
        return x

    def set_moments(self, m0: float, m1: float, m2: float, x_min: float) -> NoReturn:
        """Calculate the key properties of the PDF from the moments values."""
        """
        Calculate the key properties of the PDF from the size class moments values.
        
        Parameters
        ----------
            m0: float, > 0.
                the zeroth moment of the particle size class distribution representing
                the total number of particles in the particle population.
            m1: float, > m0.
                the first moment of the PDF representing the total number core
                elements in the particles.
            m2: float, > m1.
                the second moment of the PDF.
            x_min: float, >= 0.
                the minimum size class in the distribution. Normally this is the
                smallest particle size that is formed by the nucleation reactions.
                That is, the minumum inception size class.
        """
        # check the moments values
        if m0 > 0.0:
            if m1 / m0 <= 1.0 or m2 / m0 <= 1.0:
                msg = [
                    Color.PURPLE,
                    "The higher size moments must greater than the lower moments.",
                    Color.END,
                ]
                error_and_exit(Color.SPACE.join(msg))
        else:
            msg = [Color.PURPLE, "The zeroth size moment must > 0.", Color.END]
            error_and_exit(Color.SPACE.join(msg))
        if x_min < 0.0:
            msg = [Color.PURPLE, "The minimum size class must >= 0.", Color.END]
            error_and_exit(Color.SPACE.join(msg))
        # Calculate standard statistical properties
        self.m0 = m0
        self.mean = m1 / m0  # mean value
        self.variance = (m2 / m0) - (self.mean * self.mean)  # variance
        self.x_min = x_min  # lower bound

        # Shift the mean by the lower boundary
        # (minimum size class from the nucleation reactions)
        self.mean_shifted = self.mean - x_min
        if self.mean_shifted <= 0:
            msg = [
                Color.PURPLE,
                "The mean must be strictly greater than the minimum value x_min.",
                Color.END,
            ]
            error_and_exit(Color.SPACE.join(msg))
        # standard deviation
        self.std_dev = np.sqrt(self.variance)
        # provide basic information about the PDF
        print(f"Calculated Mean Class: {self.mean:.4f}")
        print(f"Calculated Mean Diameter: {self.class_to_diameter(self.mean)} (micron)")
        print(f"Calculated Variance: {self.variance:.4f}")
        print(f"Shifted Mean (above x_min): {self.mean_shifted:.4f}\n")

    def get_gamma_pdf(
        self, samples: int = 1000, limit: float = 5.0, mode: str = "class"
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Reconstruct the PDF by assuming the gamma function distribution."""
        """
        Applied the given PDF properties of the particle population to
        the gamma function distribution form.

        Parameters
        ----------
            samples: integer.
                number of data points used to construct the PDF.
            limit: float.
                number of standard deviations above the mean to define
                the upper bound of the PDF.
            mode: string, {"class", "volume"}.
                option to provide PDF either in particle size class or particle volume.

        Returns
        -------
            x: float, array of size ``samples``.
                particle size class [-] or volume [micron3].
            pdf_gamma: float, array of size ``samples``.
                number density.
        """
        # --- GAMMA DISTRIBUTION PARAMETERS ---
        alpha = (self.mean_shifted * self.mean_shifted) / self.variance  # Shape
        beta = self.variance / self.mean_shifted  # Scale
        print(
            f"Gamma Parameters -> Shape (alpha): {alpha:.4f}, Scale (beta): {beta:.4f}"
        )
        # Generate x (size class) values for plotting
        # (from x_min to 5 standard deviations above the mean)
        x = np.linspace(self.x_min, self.mean + limit * self.std_dev, samples)
        # Calculate PDFs (Scipy handles shifting via the 'loc' argument)
        # Multiply by m0 to scale the height if the total area m0 is not equal to 1
        pdf_gamma = self.m0 * gamma.pdf(x, a=alpha, loc=self.x_min, scale=beta)
        # convert the size class to volume
        if mode == "volume":
            x *= self.vol_per_class  # micron3
        return x, pdf_gamma

    def get_lognormal_pdf(
        self, samples: int = 1000, limit: float = 5.0, mode: str = "class"
    ) -> tuple[npt.NDArray[np.double], npt.NDArray[np.double]]:
        """Reconstruct the PDF by assuming the lognormal distribution."""
        """
        Applied the given PDF properties of the particle population to
        the lognormal distribution form.
        
        Parameters
        ----------
            samples: integer.
                number of data points used to construct the PDF.
            limit: float.
                number of standard deviations above the mean to define
                the upper bound of the PDF.
            mode: string, {"class", "volume"}.
                option to provide PDF either in particle size class or particle volume.
        
        Returns
        -------
            x: float, array of size ``samples``.
                particle size class [-] or volume [micron3].
            pdf_lognorm: float, array of size ``samples``.
                number density.
        """
        # --- LOGNORMAL DISTRIBUTION PARAMETERS ---
        # Convert shifted mean and variance to underlying normal parameters
        sigma_log_sq = np.log(1 + (self.variance / (self.mean_shifted**2)))
        sigma_log = np.sqrt(sigma_log_sq)  # Shape parameter for scipy
        mu_log = np.log(self.mean_shifted) - (sigma_log_sq / 2)
        scale_log = np.exp(mu_log)  # Scale parameter for scipy
        print(
            f"Lognormal Parameters -> Mu_log: {mu_log:.4f}, "
            f"Sigma_log: {sigma_log:.4f}\n"
        )
        # Generate x (size class) values for plotting
        # (from x_min to 5 standard deviations above the mean)
        x = np.linspace(self.x_min, self.mean + limit * self.std_dev, samples)
        # Calculate PDFs (Scipy handles shifting via the 'loc' argument)
        # Multiply by m0 to scale the height if the total area m0 is not equal to 1
        pdf_lognorm = self.m0 * lognorm.pdf(
            x, s=sigma_log, loc=self.x_min, scale=scale_log
        )
        # convert size class to volume
        if mode == "volume":
            x *= self.vol_per_class  # micron3
        return x, pdf_lognorm

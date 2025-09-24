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
    Micro mixing model for the multi-zone engine models.
"""

import copy

from ansys.chemkin.color import Color as Color
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.mixture import (
    calculate_mixture_temperature_from_enthalpy,
    interpolate_mixtures,
)
from ansys.chemkin.utilities import random, random_pick_integers
import numpy as np


class MicroMixing:
    """
    Micro mixing process
    """
    def __init__(self):
        """
        Stochastic micro mixing process module
        The event particles are assumed to have the same mass.
        """
        #
        # total number of events/particles
        self.numb_particles = 0
        # perform random micro mixing between particles
        self.delta_time = 0.0e0
        # tubulence scalar mixing time scale [sec]
        self.mixing_time_scale = 0.0e0
        # micro mixing mopdel parameter
        self.mixing_model_parameter = 1.0
        # stability limit for the maximum fraction of particles can be mixed per mixing process
        self.mixing_fraction_limit = 0.3
        # number of particle mixtures
        # by default one particle per mixture
        self._numb_mixtures = 0
        # mapping particles to mixture index {particle index: mixture index}
        # by default, {1:1, 2:2, ...}
        self.particles_map: dict[int, int] = {}
        # mapping mixture index to mixture objects for the particles {mixture index: mixture object}
        self.mixture_map: dict[int, Stream] = {}
        # number of particles per mixture [number of particles per mixture]
        # by default, [1, 1, ...]
        self.mixture_particles: list[int] = []
        # mixtures of which some particles' properties are changed {mixture index: list of changed mixture objects} 
        self.changed_mixtures: dict[int, list[Stream]] = {}
        # unit mass of the particle
        self.particle_unit_mass = 0.0e0

    def set_numb_particles(self, nparticles: int):
        """
        Set the total number of events/particles

        Parameters
        ----------
            nparticles: integer
                total number of events/particles
        """
        if nparticles > 1:
            self.numb_particles = nparticles
        else:
            msg = [
                Color.PURPLE,
                "the total number of particles must > 1.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_mixing_time_step(self, dtime: float):
        """
        Set the time step size for micro mixing process

        Parameters
        ----------
            dtime: double
                time step size (duration) of the micro mixing process [sec]
        """
        if dtime > 0.0e0:
            self.delta_time = dtime
        else:
            msg = [
                Color.PURPLE,
                "the time step size must > 0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_mixing_time_scale(self, tau: float):
        """
        Set the characteristic scalar mixing time scale

        Parameters
        ----------
            tau: double
                characteristic scalar mixing time scale [sec]
        """
        if tau > 0.0e0:
            self.mixing_time_scale = tau
        else:
            msg = [
                Color.PURPLE,
                "the scalar mixing time scale must > 0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_mixing_model_parameter(self, cmix: float):
        """
        Set the parameter of the micro mixing model

        Parameters
        ----------
            cmix: double
                micro mixing model parameter
        """
        if cmix > 0.0e0:
            self.mixing_model_parameter = cmix
        else:
            msg = [
                Color.PURPLE,
                "the model parameter must > 0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_particle_mixtures(self, particle_mixtures: list[Stream]):
        """
        Set the mixtures to represent the particle properties

        Parameters
        ----------
            particle_mixtures: list of Mixture objects
                Mixture to represent certain group of particles
        """
        # number of particle mixtures
        self._numb_mixtures = len(particle_mixtures)
        # check
        if self._numb_mixtures <= 0:
            msg = [
                Color.PURPLE,
                "the mixture list provided is empty.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # set up mixtures
        self.mixture_map.clear()
        # 1-based mixture index
        mixture_index = 1
        for m in particle_mixtures:
            # mapping mixture index to mixture objects for the particles {mixture index: mixture object}
            self.mixture_map[mixture_index] = copy.deepcopy(m)
            mixture_index += 1

    def calculate_particles_per_mixture(self):
        """
        Determine the particle count distribution among the mixtures that represent
        the gas properties of those particles
        """
        # check
        if self._numb_mixtures == 0:
            msg = [
                Color.PURPLE,
                "particle properties must be given as mixtures,",
                "use 'set_particle_mixtures()' to provide particle properties.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        if self._numb_mixtures == self.numb_particles:
            # 1 particle per mixture
            self.particles_map.clear()
            for n in range(self.numb_particles):
                # mapping particles to mixture index {particle index: mixture index}
                # 1-based index
                id = n + 1
                self.particles_map[id] = id
                # number of particles per mixture [number of particles per mixture]
                self.mixture_particles.append(1)
        else:
            # distribute the particles according to the mixture mass fraction
            mixture_mass_frac = []
            # total mass of the mixtures
            sum_mass = 0.0
            for m in self.mixture_map.values():
                # compute mixture mass
                den = m.RHO
                vol = m.volume
                mass = den * vol
                sum_mass += mass
                mixture_mass_frac.append(mass)
            # normalize to obtain mixture mass fractions
            mixture_mass_frac[:] /= sum_mass
            # unit mass of the particle
            self.particle_unit_mass = sum_mass / float(self.numb_particles)
            # find the largest mixture
            largest_zone_index = mixture_mass_frac.index(max(mixture_mass_frac))
            # initialization
            self.mixture_particles.clear()
            zone_index = 0
            sum_particles = 0
            for f in mixture_mass_frac:
                if zone_index == largest_zone_index:
                    # temporarily use "-1" as the place holder for the largest mixture
                    self.mixture_particles.append(-1)
                else:
                    # compute the particle count according to the mixture mass fraction
                    p = int(self.numb_particles * f)
                    self.mixture_particles.append(p)
                    sum_particles += p
                zone_index += 1
            # assign any particle count distribution error to the largest mixture 
            self.mixture_particles[largest_zone_index] = self.numb_particles - sum_particles
            # map the particles
            self.map_mixture_particles()

    def map_mixture_particles(self):
        """
        Map particle index to the corresponding mixture index
        """
        # set up the particles_map: dict[int, int] = {}
        # initialization
        self.particles_map.clear()
        # 1-based mixture index
        mixture_index = 1
        start = 1
        for n in self.mixture_particles:
            end = start + n
            for i in range(start, end):
                self.particles_map[i] = mixture_index
            mixture_index += 1
            start = end

    def get_source_integer(self, min_integer: int, max_integer: int) -> list[int]:
        """
        Create a list of integers as the source of the random pick process

        Parameters
        ----------
            min_integer: integer
                lower bound of the integer list
            max_integer: integer
                upper bound of the integer list

        Returns
        -------
            source_list: list of integer, dimension = (max_integer - min_integer) + 1
                list of integers in the closed interval [min_integer, max_integer]
        """
        # check
        if max_integer <= min_integer:
            msg = [
                Color.PURPLE,
                "the integer bounds are out of order,",
                "the lower bound value",
                str(min_integer),
                "must <",
                "the upper bound value",
                str(max_integer),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        source_list = []
        for i in range(min_integer, max_integer + 1):
            source_list.append(i)
        return source_list

    def update_mixtures(self):
        """
        Update the mixtures after each mixing sub process
        """
        # update the zone mixtures after the properties of some of their particles
        # have been changed due to the micro mixing process
        for mixture_index, changed_mixtures in self.changed_mixtures.items():
            # number of particels that has been modified by the micro mixing process
            particle_count = self.mixture_particles[mixture_index - 1]
            # find mass fraction of each partcle (the particles must have the same mass)
            mass_frac = 1.0e0 / float(particle_count)
            # create a working copy of the zone ixture
            zone_mixture = copy.deepcopy(self.mixture_map.get(mixture_index))
            # compute the net changes of total enthalpy and species mass fractions
            h_sum = 0.0e0
            y_sum = np.zeros_like(zone_mixture.Y, dtype=np.double)
            changed_count = 0
            #
            for m in changed_mixtures:
                changed_count += 1
                # mean molecular weight of the mixture
                mean_wt = m.WTM
                # enthalpy [erg] per particle
                # HML [erg/mol]
                h_sum += m.HML() * mass_frac / mean_wt
                # species mass fractions
                this_y = m.Y
                for k in range(len(this_y)):
                    y_sum[k] += this_y[k] * mass_frac
            # compute the portion of the mixture properties that are not affected by micro mixing
            # total mass fraction of all unchanged particles
            unchanged_frac = float(particle_count - changed_count) * mass_frac
            # the amount of total enthalpy of the unchanged particles
            unchanged_h = zone_mixture.HML() * unchanged_frac / zone_mixture.WTM
            # compute the new mixture total enthalpy (changed + unchanged)
            h_sum += unchanged_h
            # convert from [erg/g] to [erg/mol]
            h_sum *= zone_mixture.WTM
            # compute the new mixture species mass fractions (changed + unchanged)
            zone_y = zone_mixture.Y
            for k in range(len(this_y)):
                y_sum[k] += zone_y[k] * unchanged_frac
            #
            zone_mixture.resetcomposition()
            zone_mixture.Y = y_sum
            # compute the new mixture gas temperature
            # set the guessed temperature
            t_guessed = zone_mixture.temperature
            iErr = calculate_mixture_temperature_from_enthalpy(
                mixture=zone_mixture, mixtureH=h_sum, guesstemperature=t_guessed
            )
            if iErr != 0:
                msg = [
                    Color.PURPLE,
                    "failed to determine the mean mixture temperature,",
                    "error code =",
                    str(iErr),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            # update the zone mixture properties
            self.mixture_map[mixture_index] = copy.deepcopy(zone_mixture)
            # clean up
            del y_sum, zone_mixture       

    def modified_curls(self, delta_time: float, tau: float, cmix: float = 1.0) -> list[Stream]:
        """
        Modified Curl's mixing model

        Parameters
        ----------
            delta_time: double
                time step size (duration) of the micro mixing process [sec]
            tau: double
                characteristic scalar mixing time scale [sec]
            cmix: double, default = 1.0
                micro mixing model parameter

        Returns
        -------
            new_mixtures: list of Mixture objects
                zone mixtures after the micro mixing process
        """
        # check
        self.set_mixing_model_parameter(cmix)
        if self.numb_particles <= 1:
            msg = [
                Color.PURPLE,
                "the total number of particles must > 1.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        # model parameters
        self.set_mixing_model_parameter(cmix)
        self.set_mixing_time_scale(tau)
        self.set_mixing_time_step(delta_time)

        # set up particle mixtures
        if self._numb_mixtures == 0:
            # particle mixture is not provided
            msg = [
                Color.PURPLE,
                "missing particle properties,",
                "use 'set_particle_mixtures()'",
                "to provide particle properties as Mixtures",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            # set up mappings
            self.calculate_particles_per_mixture()
        # compute mixing frequency
        mixing_frequency = self.mixing_model_parameter / self.mixing_time_scale
        # compute the fraction of particles to be mixed
        mixing_frac = self.delta_time * mixing_frequency
        # create the source integer list for the random selection process
        # create source integer list
        min_integer = 1
        max_integer = self.numb_particles
        max_picks = (max_integer - min_integer) + 1
        source_list = self.get_source_integer(min_integer, max_integer)
        # determine the number of subprocesses for the mixing process
        # total number of pairs of particles to be mixed
        total_numb_mixing_pairs = int(mixing_frac * self.numb_particles) + 1
        # max number of mixing pairs per mixing process
        mixing_pairs_limit = int(self.mixing_fraction_limit * self.numb_particles)
        # the mixing process must be divided into subprocesses to maintain stability
        total_count_mixing_subprocess = total_numb_mixing_pairs // mixing_pairs_limit + 1
        # leftover pairs for the last mixing subprocess
        if total_count_mixing_subprocess > 1:
            # number of particle pairs can be mixed per mixing subprocess
            numb_mixing_particles = mixing_pairs_limit
            # leftover pairs for the last mixing subprocess
            numb_leftover_pairs = total_numb_mixing_pairs % mixing_pairs_limit
            sum_pairs = 0
            # loops over total_count_mixing_subprocess - 1
            for loop in range(total_count_mixing_subprocess - 1):
                # randomly pick particle pairs
                numb_picks = numb_mixing_particles * 2
                # perform micro mixing of the selected pairs
                self.particle_mixing_curls(numb_picks, source_list)
                # number of particle pairs mixed
                sum_pairs += numb_mixing_particles
        else:
            # perform the last (or the only) micro mixing subprocess
            sum_pairs = 0
            # leftover pairs for the last mixing subprocess
            numb_leftover_pairs = total_numb_mixing_pairs

        # perform the last (or the only) micro mixing subprocess
        if numb_leftover_pairs > 0:
            # last subprocess
            numb_picks = numb_leftover_pairs * 2
            # perform micro mixing of the selected pairs
            self.particle_mixing_curls(numb_picks, source_list)
        # compile the updated mixtures
        new_mixtures: list[Stream] = []
        for m in self.mixture_map.values():
            new_mixtures.append(copy.deepcopy(m))
        return new_mixtures

    def particle_mixing_curls(self, numb_picks: int, source_list: list[int]):
        """
        Select numb_picks random integers from the source_list and mix them two by two
        to create numb_picks new mixtures
        """
        # pick a series of random integers from the source list
        picked, unpicked = random_pick_integers(numb_picks, source_list)
        # initialization
        index = 0
        self.changed_mixtures.clear()
        # loop over all particle pairs in the mixing subprocess
        while index < numb_picks:
            # get index of the first particle A of the pair
            particle_id_A = picked[index]
            # get mixture index of particle A
            mixture_id_A = self.particles_map.get(particle_id_A, 0)
            mixtureA = self.mixture_map.get(mixture_id_A, None)
            # print(f"part A = {index}  {particle_id_A}  {mixture_id_A}")
            index += 1
            # get index of the second particle B of the pair
            particle_id_B = picked[index]
            # get mixture index of particle B
            mixture_id_B = self.particles_map.get(particle_id_B, 0)
            mixtureB = self.mixture_map.get(mixture_id_B, None)
            # print(f"part B = {index}  {particle_id_B}  {mixture_id_B}")
            index += 1
            #
            if mixture_id_A != mixture_id_B:
                # perform micro mixing of the selected pairs if they belong to different mixtures
                # find random change stride towards the mean of the two mixtures
                # each particle moves half of the stride distance
                change = 0.5 * random()
                # find the mean mixtue of the two particles
                mixtureAVE = interpolate_mixtures(mixtureB, mixtureA, ratio=0.5)
                #print(f"inter T: {mixtureA.temperature} {mixtureB.temperature} {mixtureAVE.temperature}")
                #print(f"inter X: {mixtureA.X[13]} {mixtureB.X[13]} {mixtureAVE.X[13]}")
                # modify the properties of the first particle A of the pair
                mixture_new = interpolate_mixtures(mixtureAVE, mixtureA, change)
                #print(f"change = {change}")
                #print(f"interA T: {mixtureA.temperature} {mixtureAVE.temperature} {mixture_new.temperature}")
                #print(f"interA X: {mixtureA.X[13]} {mixtureAVE.X[13]} {mixture_new.X[13]}")
                # store the modified properties in the changed dictionary
                this_list = self.changed_mixtures.get(mixture_id_A, [])
                this_list.append(copy.deepcopy(mixture_new))
                self.changed_mixtures[mixture_id_A] = copy.deepcopy(this_list)
                # modify the properties of the second particle B of the pair
                mixture_new = interpolate_mixtures(mixtureAVE, mixtureB, change)
                # store the modified properties in the changed dictionary
                this_list = self.changed_mixtures.get(mixture_id_B, [])
                this_list.append(copy.deepcopy(mixture_new))
                self.changed_mixtures[mixture_id_B] = copy.deepcopy(this_list)
        # update the mixtures after each micro mixing subprocess
        self.update_mixtures()
        # clean up
        del picked, unpicked, mixtureAVE, mixture_new, this_list

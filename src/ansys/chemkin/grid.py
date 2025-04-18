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
Chemkin grid quality control parameters.
"""

import copy

from ansys.chemkin.color import Color
from ansys.chemkin.logger import logger
import numpy as np
import numpy.typing as npt


class Grid:
    """
    Grid quality control parameters for Chemkin 1-D steady-state reactor models.
    """

    def __init__(self):
        """
        Grid quality control parameters for Chemkin 1-D steady-state reactor models.
        """
        self.max_numb_grid_points = 250
        self.max_numb_adapt_points = 10
        self.gradient = 0.1
        self.curvature = 0.5
        self.numb_grid_points = 6
        self.starting_x = 0.0
        self.ending_x = 0.0
        self.reaction_zone_center_x = 0.0
        self.reaction_zone_width = 0.0
        self.grid_profile = []
        self.numb_grid_profile = 0

    def set_numb_grid_points(self, numb_points: int):
        """
        Set the number of uniform grid points at the start of simulation.

        Parameters
        ----------
            numb_points: integer, default = 6
                number of initial grid points
        """
        if numb_points > 0:
            self.numb_grid_points = numb_points
        else:
            msg = [Color.PURPLE, "number of points must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_max_grid_points(self, numb_points: int):
        """
        Set the maximum number of grid points allowed during the solution refinement.

        Parameters
        ----------
            numb_points: integer, default = 250
                maximum number of grid points
        """
        if numb_points > 0:
            self.max_numb_grid_points = numb_points
        else:
            msg = [Color.PURPLE, "number of points must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    @property
    def start_position(self) -> float:
        """
        Get the coordinate value of the first grid point, that is, the inlet/entrance.

        Returns
        -------
            position: double
                coordinate of the first grid point [cm]
        """
        return self.starting_x

    @start_position.setter
    def start_position(self, position: float):
        """
        Reset the coordinate value of the first grid point, that is, the inlet/entrance.

        Parameters
        ----------
            position: double, default = 0.0
                coordinate of the first grid point [cm]
        """
        self.starting_x = position

    @property
    def end_position(self) -> float:
        """
        Get the coordinate value of the last grid point, that is, the outlet/exit/gap.

        Returns
        -------
            position: double
                coordinate of the last grid point [cm]
        """
        return self.ending_x

    @end_position.setter
    def end_position(self, position: float):
        """
        Set the coordinate value of the last grid point, that is, the outlet/exit/gap.

        Parameters
        ----------
            position: double
                coordinate of the last grid point [cm]
        """
        if position > self.starting_x:
            self.ending_x = position
        else:
            msg = [Color.PURPLE, "ending position must > starting position.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_reaction_zone_center(self, position: float):
        """
        Set the coordinate value of the reaction/mixing zone center.

        Parameters
        ----------
            position: double
                coordinate of center of the reaction/mixing zone [cm]
        """
        if position < self.starting_x:
            msg = [Color.PURPLE, "zone center must >= starting position.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        elif position > self.ending_x:
            msg = [Color.PURPLE, "ending position must <= ending position.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        else:
            self.reaction_zone_center_x = position

    def set_reaction_zone_width(self, size: float):
        """
        Set the width of the reaction/mixing.

        Parameters
        ----------
            size: double
                width of the reaction/mixing zone [cm]
        """
        if size > 0.0:
            self.reaction_zone_width = size
        else:
            msg = [Color.PURPLE, "zone width must > 0.0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_max_adaptive_points(self, numb_points: int):
        """
        Set the maximum number of adaptive grid points allowed per solution refinement.

        Parameters
        ----------
            numb_points: integer, default = 10
                maximum number of adapted grid points can be added
        """
        if numb_points > 0:
            self.max_numb_adapt_points = numb_points
        elif numb_points == 0:
            msg = [
                Color.PURPLE,
                "number of points must > 0.\n",
                Color.SPACEx6,
                "set 'GRID' and 'CURV' to 1 to turn OFF grid adaption.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        else:
            msg = [Color.PURPLE, "number of points must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def set_solution_quality(self, gradient: float = 0.1, curvature: float = 0.5):
        """
        Set the maximum gradient and curvature ratios in the final solution profile.
        The solver will attempt to add more grid points to improve the resolution
        of the solution profiles till both gradient and curvature ratios are below
        the specified values or till the number of grid points exceeds the maximum
        quantity allowed.

        Parameters
        ----------
            gradient: double, default = 0.1
                the maximum gradient ratio of in the final solution profiles.
            curvature: double, default = 0.5
                the maximum curvature ratio of in the final solution profiles.
        """
        # check gradient ratio value
        if gradient <= 0.0:
            msg = [Color.PURPLE, "gradient ratio must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        elif gradient > 1.0:
            msg = [Color.PURPLE, "gradient ratio must <= 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        else:
            self.gradient = gradient
        # check curvature ratio value
        if curvature <= 0.0:
            msg = [Color.PURPLE, "curvature ratio must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        elif curvature > 1.0:
            msg = [Color.PURPLE, "curvature ratio must <= 1.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        else:
            self.curvature = curvature

    def set_grid_profile(self, mesh: npt.NDArray[np.double]) -> int:
        """
        Specify the grid point coordinates of the initial grid points.

        Parameters
        ----------
            mesh: 1-D double array
                initial grid point positions [cm]

        Returns
        -------
            error code
        """
        ierror = 0
        ngrids = len(mesh)
        #
        if ngrids == 0:
            msg = [Color.PURPLE, "the 'mesh' parameter is empty.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 1
            exit()
        # check first grid point = starting position
        if mesh[0] != self.starting_x:
            msg = [
                Color.PURPLE,
                "the first grid point must match the starting position.\n",
                Color.SPACEx6,
                "the starting position =",
                str(self.starting_x),
                "\n",
                Color.SPACEx6,
                "the first mesh position =",
                str(mesh[0]),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 2
        # check last grid point = ending position
        if mesh[ngrids - 1] != self.ending_x:
            msg = [
                Color.PURPLE,
                "the last grid point must match the ending position.\n",
                Color.SPACEx6,
                "the ending position =",
                str(self.starting_x),
                "\n",
                Color.SPACEx6,
                "the last mesh position =",
                str(mesh[ngrids - 1]),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            ierror = 3

        if ierror == 0:
            self.numb_grid_profile = ngrids
            self.grid_profile = copy.deepcopy(mesh)

        return ierror

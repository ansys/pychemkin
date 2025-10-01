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
    Perfectly stirred reactor cluster, a PSR *only* network, of which the PSRs are solved
    simultaneously.
"""

import copy
from ctypes import c_double, c_int
from typing import Union

from ansys.chemkin import chemkin_wrapper
from ansys.chemkin.chemistry import (
    check_chemistryset,
    chemistryset_initialized,
    force_activate_chemistryset,
    verbose,
    verify_version,
)
from ansys.chemkin.color import Color as Color
from ansys.chemkin.inlet import Stream
from ansys.chemkin.logger import logger
from ansys.chemkin.reactormodel import Keyword
from ansys.chemkin.stirreactors.PSR import perfectlystirredreactor as PSR
from ansys.chemkin.stirreactors.openreactor import openreactor
import numpy as np


class PSRCluster(openreactor):
    """
    A cluster of perfectly-stirred reactors. The reactor network chain must consist of
    PSRs ONLY. The first PSR must have at least ONE external inlet stream, and the external
    must be attached to the last PSR of the network.
    """

    def __init__(self, config: list[PSR], label: Union[str, None] = None):
        """
        A cluster of connected steady-state constant pressure perfectly-stirred reactor (PSR) objects.
        The leading reactor (PSR #1) must have at least 'ONE' external inlet. The last reactor must
        have the ONLY external outlet of the reactor network. Recycling connection
        (connection from a downstream PSR to an upstream PSR) is allowed in this cluster configuration.

        Parameters
        ----------
            config: list of PSR objects
                a list of PSR objects representing their order and connection in the cluster.
            label: string, optional
                reactor network name
        """
        # check minimum version requirement = 2026 R1
        if not verify_version(261):
            exit()
        #
        if label is None:
            label = "PSR_cluster_" + str(len(config))
        # check chain configuration
        ireac = 0
        ierror = 0
        # total number of external inlets into the PSR chain
        self.total_external_inlets = 0
        # total external inlet mass flow rate into the PSR chain
        self.total_mass_flow_rate = 0.0
        # number of inlets per PSR
        self._ninlets = np.zeros(len(config), dtype=np.int32)
        for psr in config:
            # check valid PSR objects
            if not isinstance(psr, PSR):
                msg = [
                    Color.RED,
                    "Object #",
                    str(ireac + 1),
                    "is NOT a PSR object.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierror += 1
                continue
            else:
                # check external inlet for PSR #1 (the leading PSR)
                n = psr.number_external_inlets
                if ireac == 0 and n == 0:
                    msg = [
                        Color.RED,
                        "Leading PSR must have at least 1 external inlet.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.critical(this_msg)
                    ierror += 1
                    continue
                elif ireac == 0:
                    # initialize the entire PSR network as an open PSR
                    openreactor.__init__(
                        self, guessedmixture=psr.reactormixture, label=label
                    )
                    # use the leading PSR to set the problem type and the energy type of the network
                    self._problemtype = psr._problemtype
                    self._energytype = psr._energytype
                else:
                    # PSR problem type and energy type must be consistent within the network
                    if psr._problemtype.value != self._problemtype.value:
                        msg = [
                            Color.RED,
                            "PSR #",
                            str(ireac + 1),
                            "has a different problem type.",
                            Color.END,
                        ]
                        this_msg = Color.SPACE.join(msg)
                        logger.error(this_msg)
                        ierror += 1
                    if psr._energytype.value != self._energytype.value:
                        msg = [
                            Color.RED,
                            "PSR #",
                            str(ireac + 1),
                            "has a different energy type.",
                            Color.END,
                        ]
                        this_msg = Color.SPACE.join(msg)
                        logger.error(this_msg)
                        ierror += 1
                # total number of external inlets to the cluster
                self.total_external_inlets += n
                # net cluster external inlet mass flow rate [g/sec]
                self.total_mass_flow_rate += psr.net_mass_flowrate
                # assign the PSR index to the PSR
                psr.set_reactor_index(ireac + 1)
                # set the number of external inlets of each PSR
                self._ninlets[ireac] = n
            # move on to the next object in the list
            ireac += 1
        # error in PSR chain configuration
        if ierror != 0 or ireac != len(config):
            exit()
        # total number of PSRs in the chain
        self._npsrs = c_int(ireac)
        self.npsrs = ireac
        # mapping of PSR name and reactor index
        # {reactor name : reactor index}
        self.psr_map: dict[str, int] = {}
        # dictionary of the PSRs in the network
        # {PSR index : Reactor object}
        self.psr_objects: dict[int, PSR] = {}
        # set up internal PSR configuration
        # the PSR index is 1-based
        ipsr = 1
        for psr in config:
            self.psr_map[psr.label] = ipsr
            self.psr_objects[ipsr] = psr
            ipsr += 1
        # solution Stream objects for inter-connecting streams between the PSRs
        # array size = number of PSRs in the network
        # {PSR index : outflow Stream object}
        self.psr_solutions: dict[int, Stream] = {}
        # output unit
        self._myLOUT = c_int(159)
        # no zone
        self._nzones = c_int(0)
        # default reactor type settings
        # Perfectly-Stirred Reactor (PSR) model
        self._reactortype = c_int(2)
        # Steady-State PSR only
        self._solvertype = c_int(self.SolverTypes.get("SteadyState", 2))
        # recycling stream connectivity
        self.numb_recycling_streams = 0
        # PSR recycling connection dictionary: {psr index: [(target psr label, outflow fraction), (...)]}
        self.recycling_connections: dict[int, list[tuple[str, float]]] = {}
        # run status
        self.cluster_run_status = -100
        # PSR solution streams
        self._solution_streamarray: list[Stream] = []
        # net/total outlet mass flow rate from the PSR cluster to the surroundings
        self.outlet_mass_flow_rate = 0.0
        # heat exchange between PSRs in the network
        # number of heat exchange reactor pairs
        self.numb_heat_exchange_pairs = 0
        # heat exchange PSR pairs
        self.heat_exchange_pair: list[tuple[int, int]] = []
        # heat exchange apparent/effective heat transfer coefficients [cal/cm2-K-sec]
        self.heat_exchange_coeffcients: list[float] = []
        # heat exchange apparent/effective heat transfer surface area [cm2]
        self.heat_exchange_areas: list[float] = []
        # print cluster information
        print(f"number reactors in the network = {self.npsrs}")
        print(f"total number of external inlets = {self.total_external_inlets}")

    @property
    def numb_PSRs(self) -> int:
        """
        Get the number of PSRs in the network

        Returns
        -------
            numb_psr: integer
                number of PSRs
        """
        return self.npsrs

    def numb_external_inlet(self, name: str) -> int:
        """
        Get the number of external inlet to given reactor

        Parameters
        ----------
            name: string
                reactor label

        Returns
        -------
            numb_inlets: integer
                number of external inlets
        """
        id = self.psr_map.get(name, 0)
        if id == 0:
            msg = [Color.PURPLE, name, "is NOT in the PSR cluster.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            return 0
        return self._ninlets[id - 1]

    @property
    def total_numb_external_inlets(self) -> int:
        """
        Get the total number of external inlet to the reactor cluster

        Returns
        -------
            numb_inlets: integer
                total number of external inlets
        """
        return self.total_external_inlets

    @property
    def total_inlet_mass_flow_rate(self) -> float:
        """
        Get the total external inlet mass flow rate to the reactor cluster

        Returns
        -------
            massflowrate: double
                total external inlet mass flow rate [g/sec]
        """
        return self.total_mass_flow_rate

    def get_reactor_label(self, reactor_index: int) -> str:
        """
        Get the PSR name/label corresponding to the reactor index in the cluster.

        Parameters
        ----------
            reactor_index: integer
                reactor index

        Returns
        -------
            name: string
                reactor name/label
        """
        if self.npsrs > 0:
            # loop over all reactors
            for name, id in self.psr_map.items():
                if reactor_index == id:
                    # return the corresponding reactor name/label
                    return name
        # cannot find a reactor with the given index
        msg = [
            Color.MAGENTA,
            "reactor #",
            str(reactor_index),
            "is NOT found in the network.",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.warning(this_msg)
        return ""

    def set_recycling_stream(
        self, source_psr: str, connections: list[tuple[str, float]]
    ):
        """
        Set outflowing stream connectivity of the given PSR

        Parameters
        ----------
            source_psr: string
                label of the source PSR in the PSR network
            connections: list of tuples of string and double
                list of tuples consisting of PSR label and mass flow rate fraction
        """
        if source_psr not in self.psr_map:
            msg = [Color.PURPLE, source_psr, "is NOT in the PSR cluster.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

        this_psr = self.psr_map[source_psr]
        if this_psr in self.recycling_connections.keys():
            # nullify the existing connections
            msg = [
                Color.YELLOW,
                "Outlfow configuration of",
                source_psr,
                "is already defined,",
                "this new connection will override the existing one.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            #
            self.recycling_connections[this_psr] = []
            self.numb_recycling_streams -= 1
        # check connections
        sum = 0.0
        throughflow = False
        for outflow in connections:
            psr_name = outflow[0]
            frac = outflow[1]
            id = -1
            if psr_name not in self.psr_map:
                msg = [Color.PURPLE, psr_name, "is NOT in the PSR cluster.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                id = self.psr_map[psr_name]
            #
            if id == this_psr:
                msg = [
                    Color.PURPLE,
                    "(",
                    psr_name,
                    ",",
                    str(frac),
                    ")\n",
                    Color.SPACEx6,
                    "self recycling is NOT allowed.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            elif frac < 0.0:
                msg = [
                    Color.PURPLE,
                    "(",
                    psr_name,
                    ",",
                    str(frac),
                    ")\n",
                    Color.SPACEx6,
                    "recycling outflow fraction must >= 0.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                if id == this_psr + 1:
                    msg = [
                        Color.YELLOW,
                        "PSR",
                        source_psr,
                        "through-flow mass fraction =",
                        str(frac),
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
                    throughflow = True
                sum += frac
        #
        diff = 1.0e0 - sum
        if diff < 0.0:
            msg = [
                Color.PURPLE,
                "total outflow mass fraction",
                str(sum),
                "> 1.0",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            if throughflow:
                # through flow fraction is given
                if diff > 1.0e-6:
                    # the total outflow mass fraction is NOT summed to 1
                    msg = [
                        Color.PURPLE,
                        "total outflow mass fraction",
                        str(sum),
                        "< 1.0 with through flow",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
            else:
                msg = [
                    Color.YELLOW,
                    "PSR",
                    source_psr,
                    "through-flow mass fraction =",
                    str(diff),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)
                if this_psr == self.npsrs:
                    msg = [
                        Color.YELLOW,
                        "PSR",
                        source_psr,
                        "the through-flow is the net network outflow.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
        # set the outflow connection of this PSR
        self.recycling_connections[this_psr] = copy.deepcopy(connections)
        self.numb_recycling_streams += 1

    def set_recycling_keywords(self):
        """
        Add outflow recycling connectivity keywords
        """
        tag = "RECY"
        fourspaces = "    "
        for ipsr, connections in self.recycling_connections.items():
            # construct the recycling keyword lines
            keyline = tag + fourspaces + str(ipsr)
            for name, frac in connections:
                id = self.psr_map[name]
                addtarget = str(id) + fourspaces + str(frac)
                this_line = keyline + fourspaces + addtarget
                self.setkeyword(key=this_line, value=True)

    def add_heat_exchange(
        self,
        reactor1: str,
        reactor2: str,
        heat_transfer_coeff: float,
        heat_transfer_area: float,
    ):
        """
        Add heat exchange connection between two PSRs in the network

        Parameters
        ----------
            reactor1: string
                label of the first PSR of the pairing
            reactor2: string
                label of the second PSR of the pairing
            heat_transfer_coeff: double
                apparent/effective heat transfer coefficient between the two PSRs [cal/cm2-K-sec]
            heat_transfer_aree: double
                effective heat transfer surface area between the two PSRs [cm2]
        """
        # find the reactor index from its label
        idr1 = self.psr_map.get(reactor1, 0)
        idr2 = self.psr_map.get(reactor2, 0)
        # create heat exchange pair
        this_pair = (idr1, idr2)
        # check
        iErr = 0
        for count in range(2):
            id = this_pair[count]
            if id <= 0:
                # reactor named does not exist in the network
                iErr += 1
                msg = [
                    Color.PURPLE,
                    "cannot find reactor",
                    self.get_reactor_label(id),
                    "in the network.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
        # check heat tansfer coefficient
        if heat_transfer_coeff < 0.0:
            iErr += 1
            msg = [
                Color.PURPLE,
                "heat transfer coefficient must >= 0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        # check heat transfer area
        if heat_transfer_area < 0.0:
            iErr += 1
            msg = [
                Color.PURPLE,
                "heat transfer area must >= 0.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        if iErr > 0:
            exit()
        # set up the heat exchange pair
        self.numb_heat_exchange_pairs += 1
        self.heat_exchange_pair.append(this_pair)
        # heat transfer parameter between reactor 1 and reactor 2 [cal/cm2-K-sec]
        self.heat_exchange_coeffcients.append(heat_transfer_coeff)
        # heat transfer surface area between reactor 1 and reactor 2 [cm2]
        self.heat_exchange_areas.append(heat_transfer_area)

    def set_heat_exchange_keywords(self):
        """
        Add heat exchange connectivity keywords
        """
        tag = "QXCO"
        twospaces = "  "
        for count in range(self.numb_heat_exchange_pairs):
            p = self.heat_exchange_pair[count]
            idr1 = p[0]
            idr2 = p[1]
            h = self.heat_exchange_coeffcients[count]
            area = self.heat_exchange_areas[count]
            # construct the recycling keyword lines
            keyline = tag + twospaces + str(idr1) + twospaces + str(idr2)
            addtarget = str(h) + twospaces + str(area)
            this_line = keyline + twospaces + addtarget
            self.setkeyword(key=this_line, value=True)

    def __process_keywords(self) -> int:
        """
        Process input keywords for the reactor model

        Returns
        -------
            Error code: integer
        """
        iErr = 0
        iErrc = 0
        iErrInputs = 0
        # set_verbose(True)
        ipsr = 1
        # loop over all PSRs in the network
        for psr in self.psr_objects.values():
            # process the keywords of the PSR
            iErrc = psr.cluster_process_keywords()
            # check status
            if iErrc != 0:
                msg = [
                    Color.PURPLE,
                    "PSR number",
                    str(ipsr),
                    "\n",
                    Color.SPACEx6,
                    "failed to set up basic reactor keywords,",
                    "error code =",
                    str(iErrc),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                iErr += iErrc
                continue
            ipsr += 1
        ipsr -= 1
        iErrInputs = abs(ipsr - self.npsrs)
        iErr += iErrInputs
        # set maximum number of inlets per PSR
        max_PSR_inlets = int(np.max(self._ninlets))
        max_PSR_inlets += self.npsrs
        iErrInputs = chemkin_wrapper.chemkin.KINAll0D_SetMaxInletSize(
            self._chemset_index, c_int(max_PSR_inlets)
        )
        # flow connectivity other than the through-flow
        if self.numb_recycling_streams > 0:
            self.set_recycling_keywords()
        # heat exchange keywords
        if self.numb_heat_exchange_pairs > 0:
            self.set_heat_exchange_keywords()
        # set additional keywords
        # create input lines from additional user-specified keywords
        iErrInputs, nlines = self.createkeywordinputlines()
        if iErrInputs == 0:
            # process additional keywords in _keyword_index and _keyword_lines
            for s in self._keyword_lines:
                # convert string to byte
                line = bytes(s, "utf-8")
                # set additional keyword one by one
                iErrKey = chemkin_wrapper.chemkin.KINAll0D_SetUserKeyword(line)
                iErrInputs += iErrKey
            if iErrInputs == 0:
                if verbose():
                    msg = [
                        Color.YELLOW,
                        str(nlines),
                        "additional input lines are added.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
            else:
                msg = [
                    Color.PURPLE,
                    "failed to create additional input lines,",
                    "error code =",
                    str(iErrInputs),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
        else:
            msg = [
                Color.PURPLE,
                "failed to process additional keywords, error code =",
                str(iErrInputs),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
        return iErr + iErrInputs

    def __run_model(self) -> int:
        """
        Run the reactor model after the keywords are processed

        Returns
        -------
            Error code: integer
        """
        # run the simulation without keyword inputs
        iErr = chemkin_wrapper.chemkin.KINAll0D_Calculate(self._chemset_index)
        return iErr

    def run(self) -> int:
        """
        Generic Chemkin run reactor model method

        Returns
        -------
            Error code: integer
        """
        # initialize the PSR network model as an equivalent PSR
        # set up basic PSR parameters
        #
        # activate the Chemistry set associated with the Reactor instance
        force_activate_chemistryset(self._chemset_index.value)
        #
        iErr = chemkin_wrapper.chemkin.KINAll0D_Setup(
            self._chemset_index,
            self._reactortype,
            self._problemtype,
            self._energytype,
            self._solvertype,
            self._npsrs,
            self._ninlets,
            self._nzones,
        )
        if iErr == 0:
            # setup reactor model working arrays
            iErr = chemkin_wrapper.chemkin.KINAll0D_SetupWorkArrays(
                self._myLOUT, self._chemset_index
            )
            iErr *= 10
        if iErr != 0:
            msg = [
                Color.PURPLE,
                "failed to initialize the PSR network model",
                self.label,
                "\n",
                Color.SPACEx6,
                "error code =",
                str(iErr),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        #
        # get ready to run the reactor network model
        # initialize Chemkin-CFD-API
        msg = [
            Color.YELLOW,
            "running PSR network model",
            self.__class__.__name__,
            self.label,
            "...\n",
            Color.SPACEx6,
            "initialization =",
            str(check_chemistryset(self._chemset_index.value)),
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        if not check_chemistryset(self._chemset_index.value):
            # Chemkin-CFD-API is not initialized: reinitialize Chemkin-CFD-API
            msg = [Color.YELLOW, "initializing Chemkin ...", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            retVal = chemkin_wrapper.chemkin.KINInitialize(
                self._chemset_index, c_int(0)
            )
            if retVal != 0:
                msg = [
                    Color.RED,
                    "Chemkin-CFD-API initialization failed;",
                    "code =",
                    str(retVal),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                exit()
            else:
                chemistryset_initialized(self._chemset_index.value)

        # keyword processing
        msg = [
            Color.YELLOW,
            "processing and generating keyword inputs ...",
            Color.END,
        ]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        if Keyword.noFullKeyword:
            # use API calls
            retVal = (
                self.__process_keywords()
            )  # each PSR to perform its own keyword processing
        else:
            # use full keywords
            msg = [
                Color.RED,
                "full keyword option not available for PSR network model.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            retVal = 100
        if retVal != 0:
            msg = [
                Color.RED,
                "generating the keyword inputs,",
                "error code =",
                str(retVal),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)
            return retVal
        logger.debug("processing keywords complete")

        # run reactor network model
        msg = [Color.YELLOW, "running reactor network simulation ...", Color.END]
        this_msg = Color.SPACE.join(msg)
        logger.info(this_msg)
        # suppress text output to file
        if self.suppress_output:
            iErr = chemkin_wrapper.chemkin.KINAll0D_SuppressOutput()
        if Keyword.noFullKeyword:
            # use API calls
            retVal = self.__run_model()
        # update run status
        self.setrunstatus(code=retVal)
        msg = ["simulation completed,", "status =", str(retVal), Color.END]
        if retVal == 0:
            msg.insert(0, Color.GREEN)
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            msg.insert(0, Color.RED)
            this_msg = Color.SPACE.join(msg)
            logger.critical(this_msg)

        return retVal

    def get_cluster_run_status(self) -> int:
        """
        Get cluster run status

        Returns
        -------
            status: integer
                run status, 0=all reactor success; -100=not run; other=failed
        """
        self.cluster_run_status = self.getrunstatus(mode="silent")
        return self.cluster_run_status

    def get_reactor_stream(self, reactor_name: str) -> Stream:
        """
        Get the solution Stream object of the given reactor name/label.

        Parameters
        ----------
            reactor_name: string
                reactor name

        Returns
        -------
            solution_stream: Stream object
                solution of the reactor specified
        """
        # validate solution
        if self.get_cluster_run_status() != 0:
            msg = [
                Color.MAGENTA,
                "reactor network has NOT been solved successfully.\n",
                Color.SPACEx6,
                "please adjust reactor parameters and",
                "rerun the reactor network.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        # check reactor
        id = self.psr_map.get(reactor_name, 0)
        if id == 0:
            # cannot find a reactor with the given name
            msg = [
                Color.MAGENTA,
                "reactor",
                reactor_name,
                "is NOT found in the network.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        # prepare reactor solution
        return self._solution_streamarray[id - 1]

    def get_cluster_solutionstatus(self) -> bool:
        """
        Get the status of the post-process

        Returns
        -------
            status: boolean
                True = solution mixtures is ready,
                False = solution mixtures are yet to be processed
        """
        status = False
        if len(self._solution_streamarray) > 0:
            status = True
        return status

    def process_cluster_solution(self) -> int:
        """
        Post-process solution to extract the raw solution variable data from
        PSR cluster simulation results

        Returns
        -------
            run_status: integer
                error code: 0=success; -100=not run; other=failed
        """
        iErr = 0
        # validate solution
        if self.get_cluster_run_status() != 0:
            msg = [
                Color.MAGENTA,
                "reactor network has NOT been solved successfully.\n",
                Color.SPACEx6,
                "please adjust reactor parameters and",
                "rerun the reactor network.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        # check existing raw data
        if self.get_cluster_solutionstatus():
            msg = [
                Color.YELLOW,
                "PSR network solution has been processed before,",
                "any existing solution data will be deleted from the memory.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)

        # reset mixture solution parameters
        self._solution_streamarray.clear()

        # extract the solution mixture from the PSRs
        # create a species mass fraction array to hold the steady-state solution
        frac = np.zeros(self.numbspecies, dtype=np.double)
        #
        for ireac, psr in self.psr_objects.items():
            msg = [
                Color.YELLOW,
                "post-processing raw solution data ...\n",
                Color.SPACEx6,
                "PSR #",
                str(ireac),
                ":",
                psr.label,
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            ipsr = c_int(ireac)
            # create a Stream object to hold the mixture properties of current solution
            smixture = copy.deepcopy(psr.reactormixture)
            # get raw solution data
            temp = c_double(0.0)
            pres = c_double(0.0)
            iErr = chemkin_wrapper.chemkin.KINAll0D_GetSolution_perPSR(
                ipsr, temp, pres, frac
            )
            if iErr != 0:
                msg = [
                    Color.RED,
                    "failed to extract solution from",
                    "PSR #",
                    str(ireac),
                    ":",
                    psr.label,
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                return iErr

            # steady-state presure solution [dynes/cm2]
            smixture.pressure = pres.value
            # steady-state temperature solution [K]
            smixture.temperature = temp.value
            # set mixture composition
            if self._speciesmode == "mass":
                # mass fractions
                smixture.Y = frac
            else:
                # mole fractions
                smixture.X = frac
            # get reactor outlet mass flow rate [g/sec]
            # this is the total mass flow rate of PSR (n) and may not be the same as
            # the mass flow rate of the through flow that goes to PSR (n+1)
            exitmassflowrate = c_double(0.0)
            iErr = chemkin_wrapper.chemkin.KINAll0D_GetExitMassFlowRate_perPSR(
                ipsr, exitmassflowrate
            )
            if iErr == 0:
                smixture.mass_flowrate = max(0.0, exitmassflowrate.value)
            else:
                smixture.mass_flowrate = 0.0
                msg = [
                    Color.RED,
                    "failed to get the total outlet mass flow rate,",
                    "error code =",
                    str(iErr),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.critical(this_msg)
                return iErr
            # update Stream solution
            self._solution_streamarray.insert(ireac - 1, copy.deepcopy(smixture))
            # celan up
            del smixture

        # net outlet mass flow rate from the PSR cluster
        netexitmassflowrate = c_double(0.0)
        iErr = chemkin_wrapper.chemkin.KINAll0D_GetExitMassFlowRate(netexitmassflowrate)
        self.outlet_mass_flow_rate = netexitmassflowrate.value
        # clean up
        del frac

        return iErr

    def get_cluster_outlet_flowrate(self) -> float:
        """
        Return the net/total outlet mass flow rate to the surroundings.
        In the absence of surface chemistry, this value must be the same the
        'total_inlet_mass_flow_rate' of the cluster.

        Returns
        -------
            out_massflowrate: double
                total outlet mass flow rate [g/sec]
        """
        # validate solution
        if self.get_cluster_run_status() != 0:
            msg = [
                Color.MAGENTA,
                "reactor network has NOT been solved successfully.\n",
                Color.SPACEx6,
                "please adjust reactor parameters and",
                "rerun the reactor network.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        return self.outlet_mass_flow_rate

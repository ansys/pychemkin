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
    Hybrid reactor network comprised of a mix of open reactors such as PSR and PFR.
"""

import copy
from typing import Union

from ansys.chemkin.chemistry import Chemistry, verbose
from ansys.chemkin.color import Color as Color
from ansys.chemkin.flowreactors.PFR import PlugFlowReactor as PFR
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.inlet import adiabatic_mixing_streams, clone_stream, compare_streams
from ansys.chemkin.logger import logger
from ansys.chemkin.stirreactors.PSR import perfectlystirredreactor as PSR


class ReactorNetwork:
    """
    The hybrid reactor network allows internal recycling stream and reactor outflow splitting.
    The reactors are solved individually in terms. For network with complex internal
    connections, "Tearing points" can be manually defined and "tear stream" method is
    applied to solve the entire network iteratively.
    """

    def __init__(self, chem: Chemistry):
        """
        Create a hybrid reactor network object in which the reactors are solved individually.

        Parameters
        ----------
            chem: Chemistry set object
                Chemistry set
        """
        # check parameters
        if not isinstance(chem, Chemistry):
            msg = [
                Color.RED,
                'the parameter must be a "Chemistry Set" object',
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()
        else:
            self.network_chem = chem
        # number of reactors in the network
        self.numb_reactors = 0
        # index of the last reactor of the network
        self.last_reactor = 0
        # reactor network outlet index
        self._exit_index = -10
        # reactor network outlet symbol
        self._exit_name = "EXIT>>"
        # number of reactor network outlets
        self.numb_external_outlet = 0
        # dictionary of reactors having external outlet flow
        # {outlet index : reactor index}
        self.external_outlets: dict[int, int] = {}
        # external outlet streams
        # {outlet index : Stream object}
        self.external_outlet_streams: dict[int, Stream] = {}
        # mapping of reactor name and reactor index
        # {reactor name : reactor index}
        self.reactor_map: dict[str, int] = {}
        # dictionary of the reactors joined the network
        # {reactor index : Reactor object}
        self.reactor_objects: dict[int, Union[PSR, PFR]] = {}
        # solution Stream objects for inter-connecting streams between the reactors
        # array size = number of reactors in the network
        # {"source" reactor index : outflow Stream object}
        self.reactor_solutions: dict[int, Stream] = {}
        # outlet flow connectivity to other reactors in the network
        # {"source" reactor index : {[("target reactor index", outflow split fraction), ... ]}
        self.outflow_targets: dict[int, list[tuple[int, float]]] = {}
        # flag to set internal outflow connection
        self.outflow_altered = True
        # incoming flow from outside of the reactor network
        # {reactor index : number of external inlets}
        self.external_connections: dict[int, int] = {}
        # incoming flow collected from other reactors in the network
        # {"target" reactor index : [("source reactor index", outflow split fraction), ... ]}
        self.inflow_sources: dict[int, list[tuple[int, float]]] = {}
        # stream object representing all internal flows to the reactor (inckluding the through flow)
        # {reactor index : Stream object}
        self.internal_inflow: dict[int, Stream] = {}
        # flag indicating the inletnal inlet stream of the reactor is "connected"
        # {reactor index : True/False}
        self.internal_inflow_ready: dict[int, bool] = {}
        # number of "tearing points"
        self.numb_tearpoints = 0
        # index of the "tear stream" source reactor
        # array size = number of "tearing points" in the network
        self.tearpoint: list[int] = []
        # maximum number of tear stream iteration count
        self.max_tearloop_count = 200
        # tearing stream convergence tolerance
        self.tolerance = 1.0e-6
        # relaxation factor for tear loop iteration
        self.relaxation_factor = 1.0
        # convergence status of tear iteration
        self.tear_converged = False
        # run status
        self.network_run_status = -100

    def get_reactor_label(self, reactor_index: int) -> str:
        """
        Get the name/label of the reactor corresponding to the reactor index in the
        reactor network.

        Parameters
        ----------
            reactor_index: integer
                reactor index

        Returns
        -------
            name: string
                reactor name/label
        """
        if self.numb_reactors > 0:
            # loop over all reactors
            for name, id in self.reactor_map.items():
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

    def add_reactor(self, reactor: Union[PSR, PFR]):
        """
        Add a reactor to the network in order.
        Plan the order (sequence) of reactor addition carefully. The order of the reactors is
        somewhat important as it might affect the convergence rate of the network, especially with
        the presence of any "tearing stream".

        Parameters
        ----------
            reactor: open reactor (PSR or PFR) object
                the reactor object to be added to the network
        """
        # get the reactor name/label
        reactor_label = reactor.label
        # check reactor has already joined the network
        if reactor_label in self.reactor_map:
            # reactor already exists in the network
            msg = [Color.MAGENTA, "reactor already in the network.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            return
        else:
            # check: the first reactor must have at least one external inlet
            # by definition PFR has ONE external inlet
            if isinstance(reactor, PSR):
                # check if the reactor is a PSR
                if self.numb_reactors == 0 and reactor.numbexternalinlets <= 0:
                    msg = [
                        Color.PURPLE,
                        "reactor number 1",
                        reactor_label,
                        "has NO external inlet.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
            # add new reactor
            self.numb_reactors += 1
            # set the current reactor as the last reactor in the network
            self.last_reactor = self.numb_reactors
            # add to the map
            self.reactor_map[reactor_label] = self.numb_reactors
            # add to the reactor dictionary
            self.reactor_objects[self.numb_reactors] = reactor
            # any external inlets
            self.external_connections[
                self.numb_reactors
            ] = reactor.number_external_inlets
            #
            if verbose():
                msg = [
                    Color.YELLOW,
                    "new reactor",
                    reactor_label,
                    "is reactor number",
                    str(self.numb_reactors),
                    "of the network.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.info(this_msg)

    def add_reactor_list(self, reactor_list: list[Union[PSR, PFR]]):
        """
        Add a list of reactors to the network in order.
        Plan the order (sequence) of reactor addition carefully. The order of the reactors is
        somewhat important as it might affect the convergence rate of the network, especially with
        the presence of any "tearing stream".

        Parameters
        ----------
            reactor: list of open reactor (PSR or PFR) object
                the reactor objects to be added to the network
        """
        # add the reactors one by one in the order given in the list
        for rxtor in reactor_list:
            self.add_reactor(rxtor)

    def show_reactors(self):
        """
        Show the reactor labels in the network.
        """
        if self.numb_reactors <= 0:
            msg = [
                Color.YELLOW,
                "reactor network contains no reactor.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
        else:
            for name, id in self.reactor_map.items():
                print(f"  Reactor No. {id}: {name}")

    @property
    def number_reactors(self) -> int:
        """
        Get the number of reactors in the network

        Returns
        -------
            numb_reactors: integer
                number of reactors in the network
        """
        return self.numb_reactors

    @property
    def number_external_outlets(self) -> int:
        """
        Get the number of external outlets from the network

        Returns
        -------
            numb_outlets: integer
                number of external outlets from the network
        """
        return self.numb_external_outlet

    def show_internal_outflow_connections(self):
        """
        Show the ouflow connections to other reactors in the network.
        """
        if len(self.outflow_targets) > 0:
            for reactor_id, outflows in self.outflow_targets.items():
                print(
                    f"reactor number {reactor_id}: {self.get_reactor_label(reactor_id)}"
                )
                print("-- outflow connections --")
                for v in outflows:
                    if v[0] == self._exit_index:
                        # this is an external outlet
                        print("  *external outlet flow*")
                    else:
                        if v[0] == reactor_id + 1:
                            # it is a through flow to the downstream reactor
                            print("  *through flow*")
                        print(
                            f"  reactor no. {v[0]} {self.get_reactor_label(v[0])}: {v[1]}"
                        )
                print("-" * 10)
        else:
            # no outlet flow connection defined: chain network
            for name, id in self.reactor_map.items():
                # construct the out flow connection (which contains the through flows only)
                downstream = id + 1
                self.outflow_targets[id] = [(downstream, 1.0)]
                print(f"reactor number {id}: {name}")
                print("-- outflow connections --")
                print("  *through flow*")
                print(
                    f"  reactor no. {downstream} {self.get_reactor_label(downstream)}: 1.0"
                )
                print("-" * 10)

    def show_internal_inflow_connections(self):
        """
        Show the incoming flow connections from other reactors in the network.
        """
        if len(self.inflow_sources) > 0:
            for reactor_id, inflows in self.inflow_sources.items():
                print(
                    f"reactor number {reactor_id}: {self.get_reactor_label(reactor_id)}"
                )
                print("-- internal inlet sources --")
                for v in inflows:
                    if v[0] == reactor_id + 1:
                        print("  *through flow*")
                    print(
                        f"  reactor no. {v[0]} {self.get_reactor_label(v[0])}: {v[1]}"
                    )
                print("-" * 10)
        else:
            # no inlet flow connection defined
            for name, id in self.reactor_map.items():
                # construct the out flow connection (which contains the through flows only)
                downstream = id + 1
                self.outflow_targets[id] = [(downstream, 1.0)]
                print(f"reactor number {id}: {name}")
                print("-- internal inlet sources --")
                print("  *no internal inlet stream from other reactors*")
                print("-" * 10)

    def add_outflow_connections(
        self, source_label: str, outflow_split: list[tuple[str, float]]
    ):
        """
        Add outflow connections to other reactors in the network. The connection is given
        by a tuple consistinig of the target reactor name and the mass flow rate split fraction.
        The connection to the immediate downstream reactor (through flow) is optional.

        Parameters
        ----------
            source_label: string
                name/label of the source reactor
            outflow_split: list of tuples of (target reactor name, split fraction)
                outflow connections from the source reactor.
                target reactor name: string
                split fraction: double, <= 1.0
        """
        # initialization
        ierror = 0
        # sum of the outflow split fractions
        total_frac = 0.0
        # through flow connection is NOT defined
        thruflow = False
        # given connection table
        connect_table = []
        # get the reactor index from its name
        reactor_index = self.reactor_map[source_label]
        # reactor index of the immediate downstream reactor
        downstream = reactor_index + 1
        # check reactor label
        if source_label in self.reactor_map:
            # check outflow connection data
            for v in outflow_split:
                # split is a list of tuples of (target reactor name, split fraction)
                if v[0] == self._exit_name:
                    # this reactor has an external outlet
                    self.set_external_outlet(reactor_index)
                    id = self._exit_index
                else:
                    id = self.reactor_map.get(v[0], 0)
                frac = v[1]
                # convert the split format to (target reactor index, split fraction)
                connect_table.append((id, frac))
                #
                if id == self._exit_index:
                    # external outlet
                    pass
                elif id > self.numb_reactors or id <= 0:
                    # target reactor is not part of the network
                    msg = [
                        Color.PURPLE,
                        "target reactor",
                        v[0],
                        "is NOT in the network.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    ierror += 1
                elif id == reactor_index:
                    # recycle stream back to the source reactor is not allowed
                    msg = [
                        Color.PURPLE,
                        "outflow connection to self",
                        source_label,
                        "is not allowed.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    ierror += 1
                elif id == downstream:
                    # through flow fraction is given in the data
                    thruflow = True

                if frac < 0.0 or frac > 1.0:
                    # split fraction value out of bound
                    msg = [
                        Color.PURPLE,
                        "outflow split fraction to",
                        v[0],
                        "must 0 <= and <= 1.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    ierror += 1
                else:
                    # find the total outflow fraction
                    total_frac += frac

            # check split fraction sum
            if total_frac > 1.0:
                msg = [
                    Color.PURPLE,
                    "sum of outflow fraction =",
                    str(total_frac),
                    "which is > 1.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                ierror += 1

            if ierror > 0:
                # exit on outflow configuration error
                del connect_table
                exit()

            # check existing connection data
            if reactor_index in self.outflow_targets:
                #
                msg = [
                    Color.MAGENTA,
                    "outflow connection exists for",
                    source_label,
                    "\n",
                    Color.SPACEx6,
                    "the existing connection data will be reset.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                # clear the existing connection
                self.outflow_targets[reactor_index].clear()
            # add the split table to the connection table
            if not thruflow and total_frac < 1.0:
                # through flow is not defined in the connection data
                # add through flow connection
                connect_table.append((downstream, 1.0 - total_frac))
                self.outflow_targets[reactor_index] = connect_table
                del connect_table
            else:
                if total_frac <= 0.0:
                    msg = [
                        Color.PURPLE,
                        "sum of outflow fraction =",
                        str(total_frac),
                        "which is ~ 0.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                # add the split table to the connection table
                new_table = []
                for con in connect_table:
                    # split is a list of tuples of (target reactor index, split fraction)
                    # normalize the split fractions
                    new_frac = con[1] / total_frac
                    new_table.append((con[0], new_frac))
                self.outflow_targets[reactor_index] = copy.deepcopy(new_table)
                del connect_table
                del new_table
        else:
            # reactor does not exist in the network
            msg = [
                Color.PURPLE,
                "reactor NOT in the network.\n",
                Color.SPACEx6,
                'please "add" reactor',
                source_label,
                "to the network first",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)

    def clear_connections(self):
        """
        Clear the internal connection configurations.
        """
        # clear internal outlet split
        self.outflow_targets.clear()
        # clear internal inlet source
        self.inflow_sources.clear()
        # reset the connection flag
        self.outflow_altered = True
        # clear the network external outlet
        self.external_outlets.clear()
        self.numb_external_outlet = 0

    def remove_reactor(self, name: str):
        """
        Remove the named reactor from the network.

        Parameters
        ----------
            name: string
                reactor name/label
        """
        print(f"reactor map {self.reactor_map}")
        id = self.reactor_map.get(name, 0)
        # check reactor
        if id <= 0:
            # reactor does not exist in the network
            msg = [
                Color.MAGENTA,
                "reactor",
                name,
                "is NOT in the network.",
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            return
        else:
            # remove flow connections
            # loop over all reactor
            for i, connect_table in self.outflow_targets.items():
                #
                if i == id:
                    # skip the reactor to be removed
                    continue
                else:
                    # skip the reactor to be removed
                    count = 0
                    for v in connect_table:
                        if v[0] == id:
                            # find the record of the named reactor as the target
                            this_index = count
                            break
                        else:
                            count += 1
                    # remove the reference from the target list
                    del connect_table[this_index]

            # remove outflow connection of the named reactor
            if id in self.outflow_targets.keys():
                del self.outflow_targets[id]
            if id in self.inflow_sources.keys():
                del self.inflow_sources[id]
            self.outflow_altered = True
            # remove from the network external outlet list
            if id in self.external_outlets.values():
                n = list(self.external_outlets.values()).index(id)
                nexit = list(self.external_outlets.keys())[n]
                del self.external_outlets[nexit]
                self.numb_external_outlet -= 1
                # in case this is the last reactor of a chain network
                # add an external outlet to the second to last reactor
                if self.numb_external_outlet <= 0:
                    self.set_external_outlet(id - 1)
            # remove external inlet connection
            if id in self.external_connections.keys():
                del self.external_connections[id]
            # remove solution
            if id in self.reactor_solutions.values():
                del self.reactor_solutions[id]
            # update the last reactor index
            if id == self.last_reactor:
                self.last_reactor -= 1
            # remove reactor object
            del self.reactor_objects[id]
            # remove reactor mapping
            del self.reactor_map[name]
            #
            self.numb_reactors -= 1
            print(f"reactor {name} is removed from the network.")
            print(f"updated map = {self.reactor_map}")

    def set_reactor_outflow(self):
        """
        Set up and verify the the outlet flow connections from this reactor
        to the target reactors in the network.
        """
        # configure the outlet flow connections to other reactors in the network
        print("<<<reactor network configuration>>>")
        if len(self.outflow_targets) == 0:
            # chain reactors: R1 -> R2 -> R3 ...
            # there is no outlet flow splitting in the network
            # construct the outflow connections
            for n in self.reactor_map.values():
                # check if it is the last reactor in the network
                # i = list(self.reactor_map.values()).index(n)
                # name = list(self.reactor_map.keys())[i]
                name = self.get_reactor_label(n)
                print(f" - processing reactor # {n} {name}")
                if n == self.last_reactor:
                    # there is no more downstream reactor
                    self.set_external_outlet(n)
                    print(
                        f" - added reactor # {n} {name}"
                        + f" as the external outlets # {self.numb_external_outlet}"
                    )
                    # 100% outlet mass flow goes oout of the network
                    self.outflow_targets[n] = [(self._exit_index, 1.0)]
                else:
                    # index of the immediate downstream reactor
                    downstream = n + 1
                    # 100% outlet mass flow goes to the next reactor (through flow only)
                    self.outflow_targets[n] = [(downstream, 1.0)]
                    print(
                        " - added thru flow to "
                        + f"reactor # {downstream} {self.get_reactor_label(downstream)}"
                    )
                print(f" reactor # {n} {name}")
                print(f" updated outlet target list {self.outflow_targets.get(n, [])}")
        else:
            # network with outlet flow splitting but without recycling stream
            # R1 -> R2 -> R4 -> R5
            #  \--> R3 ---^
            # the connection has been verified when the connections were added
            print(f" *total number of reactors = {self.numb_reactors}")
            print(f" *the last reactor # {self.last_reactor}")
            print(f" *total number of external outlets = {self.numb_external_outlet}")
            print("=" * 20)
            for n in self.reactor_map.values():
                # check if it is the last reactor in the network
                name = self.get_reactor_label(n)
                print(f"   reactor # {n} {name}")
                print(
                    f"   updated outlet target list {self.outflow_targets.get(n, [])}"
                )
        # collect the incoming streams from other reactors in the network
        self.set_inflow_connections()
        #
        if verbose():
            print("<<<network internal connectivity>>>")
            for id in self.reactor_map.values():
                print(f" * reactor index {id}")
                if id in self.inflow_sources.keys():
                    print(f"   < inlet sources: {self.inflow_sources[id]}")
                else:
                    print("   = no internal inflow*")
                print(f"   > outflow split: {self.outflow_targets[id]}")
            print("=" * 20)

    def set_inflow_connections(self):
        """
        Set up the sources of the internal network inlet stream to the reactor.
        """
        # initialize the collection
        self.inflow_sources.clear()
        # loop over all target reactors
        for n in self.reactor_map.values():
            # initialize the internal inlet stream set up flag
            self.internal_inflow_ready[n] = False
            # loop over all connections from all reactors in the network
            for id, targets in self.outflow_targets.items():
                if id != n:
                    # check if this reactor is a target of other reactors
                    for v in targets:
                        if v[0] == n:
                            # part of the outflow from this reactor goes to the target reactor
                            this_list = list(self.inflow_sources.get(n, []))
                            this_list.append((id, v[1]))
                            self.inflow_sources[n] = this_list

    def set_external_outlet(self, reactor_index: int):
        """
        Add a new network external outlet to the reactor.

        Parameters
        ----------
            reactor_index: integer
                reactor index
        """
        # increase the number of netork external outlets
        self.numb_external_outlet += 1
        # add the reactor index to the external outlet dictionary
        self.external_outlets[self.numb_external_outlet] = reactor_index

    def calculate_incoming_streams(self, reactor_index: int) -> Union[Stream, None]:
        """
        Calculate the combined internal incoming streams from other reactors in the network.

        Parameters
        ----------
            reactor_index: integer
                index of the current (target) reactor

        Returns
        -------
            incoming_stream: Stream object
                the total internal stream going into the current PSR
        """
        # initialization
        initialized = False
        incoming_stream = None
        # check incoming stream sources for this reactor
        if reactor_index in self.inflow_sources.keys():
            # this reactor has incoming stream from internal sources
            this_source = self.inflow_sources.get(reactor_index, [])
            for v in this_source:
                # get the source list
                id = v[0]
                frac = v[1]
                #
                if id in self.reactor_solutions.keys():
                    # check the internal inlet stream for this reactor
                    if initialized:
                        # the total internal inlet to the current reactor is already created
                        # merge the internal flow from the source PSR into the total internal inlet stream to the current PSR
                        this_stream = copy.deepcopy(
                            self.reactor_solutions.get(id, None)
                        )
                        if this_stream is None:
                            # failed to find the solution mixture of the source reactor
                            msg = [
                                Color.PURPLE,
                                "solution mixture from reactor",
                                self.get_reactor_label(id),
                                "cannot be found.",
                                Color.END,
                            ]
                            this_msg = Color.SPACE.join(msg)
                            logger.error(this_msg)
                            exit()
                        #
                        merging_flowrate = this_stream.mass_flowrate * frac
                        this_stream.mass_flowrate = merging_flowrate
                        new_stream = adiabatic_mixing_streams(
                            this_stream, incoming_stream
                        )
                        clone_stream(new_stream, incoming_stream)
                        del this_stream
                        del new_stream
                    else:
                        # create the internal inlet stream
                        name = "from_network_internal"
                        incoming_stream = Stream(self.network_chem, label=name)
                        clone_stream(self.reactor_solutions[id], incoming_stream)
                        incoming_stream.mass_flowrate = (
                            self.reactor_solutions[id].mass_flowrate * frac
                        )
                        initialized = True
                else:
                    # skip the reactor does not have solution yet
                    print(
                        f"building internal inlet for reactor {self.get_reactor_label(reactor_index)}"
                    )
                    print(f"reactor {self.get_reactor_label(id)} not having solution")
                    print(f"source: {v}")
        else:
            # this reactor does not have any incoming internal stream
            incoming_stream = None

        return incoming_stream

    def set_internal_inlet(self, reactor_index: int) -> int:
        """
        Create or update the merged inlet stream to the reactor from the rest of
        the reactors in the network.

        Parameters
        ----------
            reactor_index: integer
                reactor index

        Returns
        -------
            status: integer
                error code
        """
        id = reactor_index
        status = 0
        # create/update the internal inlet stream properties
        inlet_stream = self.calculate_incoming_streams(id)
        #
        if inlet_stream is None:
            # try to find incoming stream for an isolated reactor
            status = 1
            #
            if id not in self.external_connections.keys():
                # the reactor does not have external inlet, either
                msg = [
                    Color.PURPLE,
                    "run failure: reactor",
                    self.get_reactor_label(id),
                    "\n",
                    Color.SPACEx6,
                    "is not connected to other reactors",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
        else:
            # update internal flow inlet
            self.internal_inflow[id] = copy.deepcopy(inlet_stream)
        del inlet_stream
        return status

    def create_internal_inlet(self, reactor_index: int):
        """
        Create a new inlet stream that combines all incoming streams from the other
        reactor network to the current reactor.

        Parameters
        ----------
            reactor_index: integer
                index of the current reactor
        """
        # calculate the total incoming streams from other reactors in the network
        status = self.set_internal_inlet(reactor_index)
        # check in the internal inlets stream has already been created
        if status == 0 and not self.internal_inflow_ready[reactor_index]:
            # need to add the internal inlet stream to the reactor
            this_rxtor = self.reactor_objects[reactor_index]
            this_rxtor.set_inlet(self.internal_inflow[reactor_index])
            # update flag
            self.internal_inflow_ready[reactor_index] = True

    def get_network_run_status(self) -> int:
        """
        Get network run status

        Returns
        -------
            status: integer
                run status, 0=all reactor success;-100=not run; other=failed
        """
        sum_status = 0
        for rxtor in self.reactor_objects.values():
            status = rxtor.getrunstatus(mode="silent")
            sum_status += status

        if sum_status == 0:
            # all reactors have been run with return code = 0
            self.network_run_status = 0
        elif sum_status > 0:
            # all reactors have been run but some with return code > 0
            self.network_run_status = sum_status
        else:
            # some reactors have NOT been run
            pass
        return self.network_run_status

    def run(self) -> int:
        """
        Solve the hybrid reactor network by solving the individual reactors in
        the sequence as they are added to the network. If there is any "tear stream"
        in the network, the solution process will be repeated till the properties of the
        "tear stream" are converged.

        Returns
        -------
            run status: integer
        """
        run_status = 0
        # construct the incoming flow connections from the outflow connection table
        if self.outflow_altered:
            self.set_reactor_outflow()
        #
        if self.numb_tearpoints == 0:
            # there is no recycling tear stream in the network
            run_status = self.run_without_tearstream()
        else:
            # there are tear streams in the network
            run_status = self.run_with_tearstream()
        return run_status

    def get_reactor_Stream(self, reactor_name: str) -> Stream:
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
        if self.get_network_run_status() != 0:
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
        id = self.reactor_map.get(reactor_name, 0)
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
        return self.reactor_solutions[id]

    def set_external_streams(self):
        """
        Set up external outlet streams.
        """
        # initialization
        self.external_outlet_streams.clear()
        # validate solution
        if self.get_network_run_status() != 0:
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

        # loop over all external outlets
        for ioutlet, nrxtor in self.external_outlets.items():
            #
            frac = 1.0
            if nrxtor in self.outflow_targets.keys():
                # get the actual fraction of flow going to the external exit
                for v in self.outflow_targets[nrxtor]:
                    if v[0] == self._exit_index:
                        frac = v[1]
            # update the external outlet properties
            this_outlet = copy.deepcopy(self.reactor_solutions[nrxtor])
            self.external_outlet_streams[ioutlet] = this_outlet
            self.external_outlet_streams[ioutlet].mass_flowrate *= frac
            del this_outlet
            # display external outlet stream information
            if verbose():
                print("-" * 10)
                print(f"network external outlet # {ioutlet} from reactor # {nrxtor}")
                print(
                    f"   mass flow rate = {self.reactor_solutions[nrxtor].mass_flowrate * frac} [g/sec]"
                )
                print(
                    f"   temperature = {self.reactor_solutions[nrxtor].temperature} [K]"
                )
                print("-" * 10)

    def get_external_stream(self, stream_index: int) -> list[Stream]:
        """
        Get the list of external outlet Stream objects.

        Parameters
        ----------
            stream_index: integer
                the index of the external outlet

        Returns
        -------
            external_stream: Stream object
                external outlet stream properties
        """
        # check external oulet setup
        if self.numb_external_outlet <= 0:
            msg = [Color.MAGENTA, "no external outlet defined.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()
        # validate solution
        if self.get_network_run_status() != 0:
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
        # prepare the outlet stream
        return self.external_outlet_streams[stream_index]

    def run_without_tearstream(self) -> int:
        """
        Run the individual reactors in the network one by one without using
        tear stream iteration.

        Returns
        -------
            run_status: integer
                error code
        """
        # initialization
        status = 0
        # run the reactors in order
        for id, rxtor in self.reactor_objects.items():
            # construct through flow inlet
            # check if the reactor has any incoming connection from the other reactors
            if id in self.inflow_sources.keys():
                self.create_internal_inlet(id)

            if isinstance(rxtor, PSR):
                # if the reactor is a PSR, get a better guessed solution
                this_inlet = self.internal_inflow.get(id, None)
                if this_inlet is not None:
                    rxtor.reset_estimate_composition(this_inlet.X, mode="mole")
            # run the individual reactor model
            reactor_run_status = rxtor.run()
            # check status
            if reactor_run_status != 0:
                msg = [
                    Color.RED,
                    "run failure: reactor",
                    self.get_reactor_label(id),
                    "\n",
                    Color.SPACEx6,
                    "error code =",
                    str(reactor_run_status),
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.error(this_msg)
                exit()
            else:
                # process the solution and add the solution stream to the outflow stream list
                self.reactor_solutions[id] = rxtor.process_solution()
            status += reactor_run_status

        # construct the outlet stream properties
        self.set_external_streams()

        return status

    def run_with_tearstream(self) -> int:
        """
        Run the individual reactors in the network one by one with
        tear stream iteration.

        Returns
        -------
            run_status: integer
                error code
        """
        # initialization
        status = 0
        self._tear_converged = False
        # storage of the last reactor network solution
        # { reactor index : Stream object}
        last_reactor_solutions: dict[int, Stream] = {}
        # tear loop counter
        loop_count = 0
        #
        while not self.tear_converged and self.check_iteration_count(loop_count):
            # start of the tear stream iteration to solve the reactor network
            print(f"<---- running tear loop # {loop_count} ---->")
            #
            for id, rxtor in self.reactor_objects.items():
                # loop over all reactors
                # check reactor inlet configuration
                if rxtor.number_external_inlets == 0:
                    # the reactor has no external inlet
                    # check if the reactor has any incoming connection from the other reactors
                    # if possible, creat the internal flow inlet
                    if id in self.inflow_sources.keys():
                        self.create_internal_inlet(id)
                    else:
                        # network configuration error
                        msg = [
                            Color.PURPLE,
                            "run failure: reactor",
                            self.get_reactor_label(id),
                            "\n",
                            Color.SPACEx6,
                            "might have faulty connection configuration",
                            Color.END,
                        ]
                        this_msg = Color.SPACE.join(msg)
                        logger.error(this_msg)
                        exit()
                # run the individual reactor model in order
                reactor_run_status = rxtor.run()
                # check status
                if reactor_run_status != 0:
                    msg = [
                        Color.RED,
                        "run failure: reactor",
                        self.get_reactor_label(id),
                        "\n",
                        Color.SPACEx6,
                        "error code =",
                        str(reactor_run_status),
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.error(this_msg)
                    exit()
                else:
                    # process the solution and add the solution stream to the outflow stream list
                    self.reactor_solutions[id] = rxtor.process_solution()

                # overall run status
                status += reactor_run_status

            # tear loop residual
            loop_residual = 0.0
            #
            for id, rxtor in self.reactor_objects.items():
                # construct the through flow inlet for the next iteration
                # Note: the first iteration does not include any recycling/tear stream
                # check if the reactor has any incoming connection from the other reactors
                if id in self.inflow_sources.keys():
                    # the reactor has incoming connection from the other reactors
                    self.create_internal_inlet(id)
                    if id in self.internal_inflow.keys():
                        # reset the inlet
                        rxtor.reset_inlet(self.internal_inflow[id])
                if isinstance(rxtor, PSR):
                    # if the reactor is a PSR, get a better guessed solution
                    this_inlet = self.internal_inflow.get(id, None)
                    if this_inlet is not None:
                        rxtor.reset_estimate_composition(this_inlet.X, mode="mole")
                        rxtor.set_estimate_conditions(
                            option="TP", guess_temp=rxtor.temperature
                        )

                # update the previous reactor solutions
                if id in last_reactor_solutions.keys():
                    # previous reactor solution already exists
                    # get residual for this tear point (reactor)
                    stream_new = self.reactor_solutions[id]
                    stream_old = last_reactor_solutions[id]
                    # check tear point residuals
                    if id in self.tearpoint:
                        # compare last and current reactor solutions
                        converged, residual = self.check_tearstream_convergence(
                            stream_old, stream_new
                        )
                        # update overall residual
                        loop_residual = max(loop_residual, residual)
                        # check flow rate change in the tear stream
                        flow_residual = abs(
                            self.internal_inflow[id].mass_flowrate
                            - stream_old.mass_flowrate
                        )
                        flow_residual /= stream_old.mass_flowrate
                        loop_residual = max(loop_residual, flow_residual)
                        converged = converged and loop_residual <= self.tolerance
                        # check overall loop convergence
                        self.tear_converged = self.tear_converged or converged
                    # update internal stream properties
                    updated_stream = self.update_tear_solution(stream_new, stream_old)
                    # clone the streams
                    clone_stream(updated_stream, self.reactor_solutions[id])
                    clone_stream(updated_stream, last_reactor_solutions[id])
                    del updated_stream
                else:
                    # there is no previous reactor solution
                    # create storage for the last reactor solution and save the current reactor solution there
                    last_reactor_solutions[id] = copy.deepcopy(
                        self.reactor_solutions[id]
                    )
                    self.tear_converged = False

            # increment the tear loop iteration count
            print(f">---- loop {loop_count}: max residual = {loop_residual} ----<")
            print()
            loop_count += 1

        # done tear iteration loop
        if self.tear_converged:
            # convergence information
            msg = [
                Color.YELLOW,
                "the reactor network is converged.\n",
                Color.SPACEx6,
                "number of tear iteration = ",
                str(loop_count),
                "\n",
                Color.SPACEx6,
                "max tear loop residual =",
                str(loop_residual),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.info(this_msg)
            # construct the outlet stream properties
            self.set_external_streams()
        else:
            # go over the mximum iteration count
            msg = [
                Color.PURPLE,
                "failure to solve the reactor network:\n",
                Color.SPACEx6,
                "max tear iteration count reached",
                str(self.max_tearloop_count),
                "\n",
                Color.SPACEx6,
                "max tear loop residual =",
                str(loop_residual),
                Color.END,
            ]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            status = 10
        # clear up
        del last_reactor_solutions

        return status

    # tear stream utilities
    def remove_tearpoint(self, reactor_name: str):
        """
        Remove the tear point from the list.

        Parameters
        ----------
            reactor_name: string
                reactor name/label
        """
        # check reactor has already joined the network
        #
        reactor_index = self.reactor_map.get(reactor_name, 0)
        if reactor_index > 0:
            # check if it is in the list
            if reactor_index in self.tearpoint:
                # delete the tear point
                self.tearpoint.remove(reactor_index)
                # decrease the tear point count
                self.numb_tearpoints -= 1
            else:
                msg = [Color.MAGENTA, "reactor is NOT a tear point.", Color.END]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                exit()
        else:
            # reactor does not exist in the network
            msg = [Color.MAGENTA, "reactor is NOT in the network.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()

    def add_tearingpoint(self, reactor_name: str):
        """
        Add a new tear point to the list.

        Parameters
        ----------
            reactor_name: string
                reactor name/label
        """
        # check reactor has already joined the network
        reactor_index = self.reactor_map.get(reactor_name, 0)
        if reactor_index > 0:
            # valid reactor name
            # check if the reactor is already declared as a tear point
            if reactor_index in self.tearpoint:
                # reactor already decalred as a tear point
                msg = [
                    Color.MAGENTA,
                    "reactor",
                    reactor_name,
                    "already declared as a tear point.",
                    Color.END,
                ]
                this_msg = Color.SPACE.join(msg)
                logger.warning(this_msg)
                exit()
            else:
                # add new tear point
                self.tearpoint.append(reactor_index)
                # increase the tear point count
                self.numb_tearpoints += 1
                #
                if verbose():
                    msg = [
                        Color.YELLOW,
                        "reactor",
                        reactor_name,
                        "is declared as tear point no.",
                        str(self.numb_tearpoints),
                        "of the network.",
                        Color.END,
                    ]
                    this_msg = Color.SPACE.join(msg)
                    logger.info(this_msg)
        else:
            # reactor does not exist in the network
            msg = [Color.MAGENTA, "reactor is NOT in the network.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.warning(this_msg)
            exit()

    def set_tear_tolerance(self, tol: float = 1.0e-6):
        """
        Set the relative tolerance to test the tear stream convergencce.

        Parameters
        ----------
            tol: double, default = 1.0e-6
                relative tolerance
        """
        if tol > 0.0:
            self.tolerance = tol
        else:
            msg = [Color.PURPLE, "tolerance must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def set_tear_iteration_limit(self, max_count: int):
        """
        Set the maximum number of tear loop iterations.

        Parameters
        ----------
            max_count: integer
                tear loop iteration limit
        """
        if max_count > 0:
            self.max_tearloop_count = max_count
        else:
            msg = [Color.PURPLE, "iteration limit must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def check_iteration_count(self, count: int) -> bool:
        """
        Check the iteration count for over the set limit.

        Parameters
        ----------
            count: integer
                current tear loop iteration count

        Returns
        -------
            status: boolean
                True=the current count is under the limit
                False=the current count is over the limit
        """
        if count <= self.max_tearloop_count:
            return True
        else:
            return False

    def set_relaxation_factor(self, relax: float):
        """
        Set the relaxation factor when updating the tear stream
        properties from their values of the previous iteration step.

        Parameters
        ----------
            relax: double
                iteration relaxation factor
        """
        if relax > 0.0:
            self.relaxation_factor = relax
        else:
            msg = [Color.PURPLE, "relaxation factor must > 0.", Color.END]
            this_msg = Color.SPACE.join(msg)
            logger.error(this_msg)
            exit()

    def check_tearstream_convergence(self, streamA, streamB) -> tuple[bool, float]:
        """
        Compare the last and the current reactor solution at the tear point.

        Parameters
        ----------
            streamA: Stream object
                reactor solution to be compared against
            streamB: Stream object
                reactor solution used for the comparison

        Returns
        -------
            converged: boolean
                True=the two streams are close within the relative tolerance
                False=the differences of the two streams are greater than the relative tolerance
            residual: double
                the maximum relative difference between the properties of the streams
        """
        # compare the streams
        converged, max_atol, max_rtol = compare_streams(
            streamA, streamB, atol=1.0e-8, rtol=self.tolerance
        )
        return converged, max_rtol

    def update_tear_solution(self, new_stream: Stream, old_stream: Stream) -> Stream:
        """
        Update the old/temporary tear point stream properties (with the used of a relaxation factor).

        Parameters
        ----------
            new_stream: Stream object
                tear point stream properties from the current iteration step
            old_stream: Stream object
                tear point stream properties from the last iteration step

        Returns
        -------
            updated_stream: Stream object
                updated tear point stream properties
        """
        #
        updated_stream = copy.deepcopy(new_stream)
        # relaxation factor
        fac = self.relaxation_factor
        # compute the new property values
        # Note: the new stream stays the same if the relaxation factor = 1.0
        diff = new_stream.temperature - old_stream.temperature
        updated_stream.temperature = old_stream.temperature + fac * diff
        #
        diff = new_stream.mass_flowrate - old_stream.mass_flowrate
        updated_stream.mass_flowrate = old_stream.mass_flowrate + fac * diff
        #
        old_fractions = copy.deepcopy(old_stream.Y)
        new_fractions = copy.deepcopy(new_stream.Y)
        for k in range(len(old_fractions)):
            diff = new_fractions[k] - old_fractions[k]
            new_fractions[k] = old_fractions[k] + fac * diff

        updated_stream.Y = new_fractions
        # clean up
        del old_fractions
        del new_fractions
        return updated_stream

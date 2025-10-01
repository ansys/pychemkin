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
.. _ref_SI_engine_knock:

============================
SI engine knocking analysis
============================

The 0-D SI engine model can be used to predict the occurrence of engine knock. There is no consensus on
the determination of the onset of engine knock. In this example, the engine knock is defined as when the end gas
(the gas mixture in the unburned zone) auto-ignites. And the knock intensity is associated to the peak value of the
total thermicity of the end gas, that is, the higher the peak total thermicity the more intense the engine knock.

There are many parameters that have impact on the occurrence of engine knock, for example, the fuel composition,
the start of combustion timing, the engine speed, the wall heat loss and etc. This tutorial focuses on the effect of
start of combustion timing (or spark timing) on the engine knock. It also shows the steps of setting up an SI engine
parameter study for knock analysis. The predicted knocking occurrence (represented by the end gas autoignition timing) are
presented as a function of the start of combustion (SOC) crank angle. Detailed description of the 0-D SI engine model can be
found in the "0-D engine model" sections of the "Chemkin Theory" manual, respectively. You can use the
``ansys.chemkin.manuals()`` method to get to the latest version of the Chemkin online manuals.

The parameter study is performed in the multi-process (or multi-core) mode by using ``ProcessPoolExecutor``
from the ``concurrent.futures`` package. If the computer has enough number of CPU cores, the parameter
study can be run in parallel with each case running on a designated CPU core.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_SI_engine_knock.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

from concurrent.futures import ProcessPoolExecutor
import os
import time

import ansys.chemkin as ck  # Chemkin

# chemkin spark ignition (SI) engine model (transient)
from ansys.chemkin.engines.SI import SIengine
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger
from ansys.chemkin.utilities import workingFolders
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching

# check working directory
current_dir = os.getcwd()
logger.debug("working directory: " + current_dir)
# set verbose mode
ck.set_verbose(True)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = True


#######################################
# Create an SI engine calculator class
# =====================================
# Create a local class that wraps around the actual ``SIengine`` class
# to make the setup of the SI engine knock analysis parameter study
# more convenient.
#
# .. note::
#   The auto-ignition detection must be turned on to obtain the correct
#   knock occurrence crank angle. Use the ``set_ignition_delay()`` method to
#   specify the ignition definition criterion and to turn on the ignition delay time
#   detection.
#
class SIengineCalculator:
    """
    SI engine calculator with fixed set up parameters
    """

    def __init__(self, fresh_mixture: Stream, index: int, start_combustion: float):
        """
        SI engine calculator that instantiates an SI engine object with
        the given fresh (unburnt) mixture condition and predefined engine
        parameters.

        Parameters
        ----------
            fresh_mixture: Mixture object
                the initial/fresh/unburnt condition
            index: integer
                run index of this SI engine calculator
            start_combustion: double
                start of combustion crank angle [CA]
        """
        # instantiate an SI engine object
        # set up the run and working directory name
        self.name = "SI_engine_" + str(index)
        # instantiate the SIengine object for this run
        self.SIcalculator = SIengine(fresh_mixture, label=self.name)
        # run status
        self.runstatus = -100
        # calculated engine knock occurrence crank angle [CA]
        self.knock_CA = -720.0
        # calculated maximum thermicity value in the unburned zone [1/sec]
        self.max_sigma = 0.0
        # peak end gas temperature [K]
        self.max_temp = 0.0
        # IMEP
        self.IMEP = 0.0
        # set the required premixed flame model parameters
        self.set_engine_parameters()
        # set burn profile
        # use fixed Wiebe function
        wiebe_b = 7.0
        wiebe_n = 4.0
        duration = 45.6  # degrees
        self.set_burn_profile(start_combustion, duration, wiebe_n, wiebe_b)

    def run(self):
        """
        Run the SI engine calculation in a separate working directory
        """
        # run the flame speed calculation
        self.runstatus = self.SIcalculator.run()
        # extract the laminar flame speed from the solution
        if self.runstatus == 0:
            # postprocess the unburned zone solution
            unburnedzone = 1
            # burnedzone = 2
            iErr = self.SIcalculator.process_engine_solution(zoneID=unburnedzone)
            if iErr == 0:
                # get the engine knock (end-gas auto-ignition) occurrence crank angle [CA]
                # if no engine knock, the value = -720 CA
                # because the memory is shared, it must be done as soon as the run is finished
                self.knock_CA = self.SIcalculator.check_engine_knock()
                temperature = self.SIcalculator.get_solution_variable_profile(
                    "temperature"
                )
                thermicity = self.SIcalculator.get_solution_variable_profile(
                    "thermicity"
                )
                self.max_temp = np.max(temperature)
                self.max_sigma = np.max(thermicity)
                self.IMEP = self.SIcalculator.get_engine_IMEP()
                del thermicity, temperature
            else:
                self.runstatus = iErr

    def set_engine_parameters(self):
        """
        Set up basic engine parameters
        """
        # Set up basic engine parameters
        # These engine parameters are used to describe the cylinder volume during the
        # simulation. The ``starting_CA`` argument should be the crank angle corresponding
        # to the cylinder IVC. The ``ending_CA`` argument is typically the EVO crank angle.
        # cylinder bore diameter [cm]
        self.SIcalculator.bore = 8.5
        # engine stroke [cm]
        self.SIcalculator.stroke = 10.82
        # connecting rod length [cm]
        self.SIcalculator.connecting_rod_length = 17.853
        # compression ratio [-]
        self.SIcalculator.compression_ratio = 12
        # engine speed [RPM]
        self.SIcalculator.RPM = 1200
        # simulation start CA [degree]
        self.SIcalculator.starting_CA = -120.2
        # simulation end CA [degree]
        self.SIcalculator.ending_CA = 139.8
        # wall heat transfer parameters
        # Set up engine wall heat transfer model
        heattransferparameters = [0.1, 0.8, 0.0]
        # set cylinder wall temperature [K]
        Twall = 434.0
        self.SIcalculator.set_wall_heat_transfer(
            "dimensionless", heattransferparameters, Twall
        )
        # in-cylinder gas velocity correlation parameter (Woschni)
        # [<C11> <C12> <C2> <swirl ratio>]
        GVparameters = [2.28, 0.318, 0.324, 0.0]
        self.SIcalculator.set_gas_velocity_correlation(GVparameters)
        # set piston head top surface area [cm2]
        self.SIcalculator.set_piston_head_area(area=56.75)
        # set cylinder clearance surface area [cm2]
        self.SIcalculator.set_cylinder_head_area(area=56.75)
        # other engine model parameters
        # set tolerances in tuple: (absolute tolerance, relative tolerance)
        self.SIcalculator.tolerances = (1.0e-10, 1.0e-6)
        # turn on the force non-negative solutions option in the solver
        self.SIcalculator.force_nonnegative = False
        # set minimum zonal mass [g]
        self.SIcalculator.set_minimum_zone_mass(minmass=1.0e-5)
        # set the maximum solver time step
        self.SIcalculator.max_CAstep = 0.1
        # set the number of crank angles between saving solution
        self.SIcalculator.CAstep_for_saving_solution = 0.5
        # set the number of crank angles between printing solution
        self.SIcalculator.CAstep_for_printing_solution = 10.0
        # turn ON adaptive solution saving
        self.SIcalculator.adaptive_solution_saving(mode=True, steps=20)
        # set to detect auto-ignition
        self.SIcalculator.set_ignition_delay(method="Thermicity_peak")
        # suppress text output
        self.SIcalculator.stop_output()

    def set_burn_profile(
        self,
        start_combustion: float,
        burn_duration: float,
        wiebe_n: float,
        wiebe_b: float,
    ):
        """
        Set up the fuel burn rate profile

        Parameters
        ----------
            start_combustion: double
                start of combustion (SOC) crank angle [CA]
            burn_duration: double
                burn duration in crank angle [CA]
            wiebe_n: double
                Wiebe function n parameter value
            wiebe_b: double
                Wiebe function b parameter value
        """
        # start of combustion CA
        self.SIcalculator.set_burn_timing(SOC=start_combustion, duration=burn_duration)
        self.SIcalculator.wiebe_parameters(n=wiebe_n, b=wiebe_b)


#############################################################################
# Set up the SI engine knocking analysis parameter study for multi-processing
# ===========================================================================
# Create ``SIengineCalculator`` object with different start of combustion (SOC) crank angle.
# Each object represents one parameter study case and will be run on an designated cpu core
# when the parameter study is executed.


def SIengine_run(case: tuple[int, float]) -> tuple[float, float, float, float, float]:
    """
    Set up the parameter study runs for multi-processing.

    Parameters
    ----------
        case: tuple (integer, double)
            SI engine calculation case condition: case index and start of combustion crank angle [CA]

    Returns
    -------
        startCA: double
            the parameter: start of combustion crank angle [degree CA]
        knock_CA: double
            crank angle when engine knock occurs [degree CA]
        max_sigma: double
            the peak total thermicity of the end gas [1/sec]
        max_temp: double
            peak end gas temperature [K]
        IMEP: double
            indicated mean effective pressure [bar]
    """
    ##########################################
    # Create an instance of the Chemistry Set
    # ========================================
    # For PRF, the encrypted 14-component gasoline mechanism, ``gasoline_14comp_WBencrypted.inp``,
    # is used. The chemistry set is named ``gasoline``.
    #
    # .. note::
    #   Because this gasoline mechanism does not come with any transport data, you do not need to provide
    #   a transport data file.
    #
    # case index
    index = case[0]
    # start of combustion crank angle
    startCA = case[1]
    # create and change to the working directory for this run
    name = "SI_engine_" + str(index)
    work_folder = workingFolders(name, current_dir)
    # set mechanism directory (the default Chemkin mechanism data directory)
    data_dir = os.path.join(ck.ansys_dir, "reaction", "data")
    mechanism_dir = data_dir
    # create a chemistry set based on the gasoline 14 components mechanism
    MyGasMech = ck.Chemistry(label="Gasoline")
    # set mechanism input files
    # including the full file path is recommended
    MyGasMech.chemfile = os.path.join(mechanism_dir, "gasoline_14comp_WBencrypt.inp")

    #######################################
    # Preprocess the gasoline chemistry set
    # =====================================

    # preprocess the mechanism files
    iError = MyGasMech.preprocess()

    ################################################
    # Set up the stoichiometric gasoline-air mixture
    # ==============================================
    # You must set up the stoichiometric gasoline-air mixture for the subsequent
    # SI engine calculations. Here the ``X_by_Equivalence_Ratio()`` method is used.
    # You create the ``fuel`` and the ``air`` mixtures first. You then define the
    # *complete combustion product species* and provide the *additives* composition if applicable.
    # Finally, you simply set ``equivalenceratio=1`` to create the stoichiometric
    # gasoline-air mixture.
    #
    # For PRF 90 gasoline, the recipe is ``[("ic8h18", 0.9), ("nc7h16", 0.1)]``.

    # create the fuel mixture
    fuelmixture = ck.Mixture(MyGasMech)
    # set fuel composition
    fuelmixture.X = [("ic8h18", 0.6), ("nc7h16", 0.0), ("mch", 0.1), ("c6h5c3h7", 0.3)]
    # setting pressure and temperature is not required in this case
    fuelmixture.pressure = 2.5 * ck.Patm
    fuelmixture.temperature = 473.0

    # create the oxidizer mixture: air
    air = ck.Mixture(MyGasMech)
    air.X = [("o2", 0.21), ("n2", 0.79)]
    # setting pressure and temperature is not required in this case
    air.pressure = fuelmixture.pressure
    air.temperature = fuelmixture.temperature

    # products from the complete combustion of the fuel mixture and air
    products = ["co2", "h2o", "n2"]
    # species mole fractions of added/inert mixture. You can also create an additives mixture here.
    add_frac = np.zeros(MyGasMech.KK, dtype=np.double)  # no additives: all zeros

    # create the unburned fuel-air mixture
    fresh = ck.Mixture(MyGasMech)

    # mean equivalence ratio
    equiv = 1.0
    iError = fresh.X_by_Equivalence_Ratio(
        MyGasMech, fuelmixture.X, air.X, add_frac, products, equivalenceratio=equiv
    )

    ##########################################################
    # Specify pressure and temperature of the fuel-air mixture
    # ========================================================
    # Since you are going to use ``fresh`` fuel-air mixture to instantiate
    # the engine object later, setting the mixture pressure and temperature
    # is equivalent to setting the initial temperature and pressure of the
    # engine cylinder.
    fresh.temperature = fuelmixture.temperature
    fresh.pressure = fuelmixture.pressure

    ###########################################
    # Add EGR to the fresh fuel-air mixture
    # =========================================
    # Many engines have the configuration for exhaust gas recirculation (EGR). Chemkin
    # engine models let you add the EGR mixture to the fresh fuel-air mixture entering
    # the cylinder. If the engine you are modeling has EGR, you should have the EGR ratio, which
    # is generally the volume ratio of the EGR mixture and the fresh fuel-air mixture.
    # However, because you know nothing about the composition of the exhaust gas, you cannot simply
    # combine these two mixtures. In this case, you use the ``get_EGR_mole_fraction()`` method to estimate
    # the major components of the exhaust gas from the combustion of the fresh fuel-air mixture. The
    # ``threshold=1.0e-8`` parameter tells the method to ignore any species with a mole fraction below
    # the threshold value. Once you have the EGR mixture composition, use the ``X_by_Equivalence_Ratio()``
    # method a second time to re-create the fuel-air mixture ``fresh`` with the original
    # ``fuelmixture`` and ``air`` mixtures, along with the EGR composition that you just got as the
    # *additives*.
    EGRratio = 0.25
    # compute the EGR stream composition in mole fractions
    add_frac = fresh.get_EGR_mole_fraction(EGRratio, threshold=1.0e-8)
    # recreate the initial mixture with EGR
    iError = fresh.X_by_Equivalence_Ratio(
        MyGasMech,
        fuelmixture.X,
        air.X,
        add_frac,
        products,
        equivalenceratio=equiv,
        threshold=1.0e-8,
    )
    if iError != 0:
        print("error...creating the initial fuel-oxidizer mixture.")
        exit()
    ##############################
    # Set up the SI engine knock case
    # ============================
    # create an SI engine calculation instance
    this_run = SIengineCalculator(fresh, index=index, start_combustion=startCA)
    # run the case
    this_run.run()
    # change back to the original top folder
    work_folder.done()
    return (
        startCA,
        this_run.knock_CA,
        this_run.max_sigma,
        this_run.max_temp,
        this_run.IMEP,
    )


#########################################
# Set up and start the multi-process runs
# =======================================
# Use the ``ProcessPoolExecutor()`` method to assign the SI engine runs.
# In this project, each SI engine run/case runs on an exclusive cpu core.
# The first parameter of ``map()`` method should be ``SIengine_run()``.
# Make the ``SIengine_run()`` method return the required parameter and results
# to make post-processing the results from the parameter study easier.
#
# .. note::
#   ``if __name__ == '__main__'`` is required when ``multiprocessing``
#   package is used.
#

if __name__ == "__main__":
    # number of available cpu cores
    numb_cores = max(os.cpu_count(), 1)
    # set the number of worker cpu cores to be used by the parameter study
    numb_workers = numb_cores - 2
    numb_workers = max(numb_workers, 1)
    # start the multi-process
    # starting value for the start of combustion crank angle
    start_combustion = -10.5
    startCA = start_combustion
    # increment
    dCA = 2.0
    # number of cases to run in the parameter study
    numb_cases = 5
    # create the cases
    SI_cases: list[tuple[int, float]] = []
    for i in range(numb_cases):
        SI_cases.append((i, startCA))
        startCA += dCA
    # set the start wall time
    start_time = time.time()
    #
    numb_workers = min(numb_workers, numb_cases)
    knock_CA: float = []
    knock_intensity: float = []
    soc_CA: float = []
    knock_CA_soc: float = []
    engine_IMEP: float = []
    peak_temp: float = []
    with ProcessPoolExecutor(max_workers=numb_workers) as e:
        for ret_value in e.map(SIengine_run, SI_cases):
            # results returned by the SI engine calculator
            # start of combustion CA
            soc_CA.append(ret_value[0])
            # knock CA (end gas auto-ignition CA)
            if ret_value[1] > -720.0:
                knock_CA.append(ret_value[1])
                knock_CA_soc.append(ret_value[0])
            # knock intensity as measured by the peak total thermicity
            # of the end gas in the unburned zone [1/sec]
            knock_intensity.append(ret_value[2])
            # peak end gas temperature [K]
            peak_temp.append(ret_value[3])
            # IMEP [bar]
            engine_IMEP.append(ret_value[4])
    # compute the total runtime
    runtime = time.time() - start_time
    print()
    print(f"total simulation duration: {runtime} [sec]")
    print()

    ###########################################
    # Plot the premixed flame solution profiles
    # =========================================
    # Plot the predicted flame speeds against the experimental data.
    plt.subplots(2, 2, sharex="col", figsize=(12, 6))
    plt.subplot(221)
    plt.plot(knock_CA_soc, knock_CA, linestyle="-", marker="^", color="blue")
    plt.ylabel("Knock Crank Angle [degree CA]")
    plt.subplot(222)
    plt.plot(soc_CA, knock_intensity, linestyle="-", marker="o", color="green")
    plt.ylabel("Knock Intensity [1/sec]")
    plt.subplot(223)
    plt.plot(soc_CA, engine_IMEP, linestyle="-", marker="s", color="magenta")
    plt.ylabel("IMEP [bar]")
    plt.xlabel("Start of Combustion [degree CA]")
    plt.subplot(224)
    plt.plot(soc_CA, peak_temp, linestyle="-", marker="v", color="red")
    plt.ylabel("Peak End-Gas Temperature [K]")
    plt.xlabel("Start of Combustion [degree CA]")
    # plot results
    if interactive:
        plt.show()
    else:
        plt.savefig("plot_SI_engine_knock.png", bbox_inches="tight")

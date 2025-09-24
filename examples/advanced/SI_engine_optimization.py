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
.. _ref_SI_engine_optimization:

===========================================================
Find the optimal operating parameters for a given SI engine
===========================================================

Finding the optimal solution is one of the prominent parts of the research and development activities.
When it comes to *PyChemkin*, reaction mechanism (rate parameters) optimization, fuel surrogate blend
optimization, and process/reactor parameters optimization are some of the common applications.

This example project describes the steps of integrating *PyChemkin* and *pymoo*, a well-known
Python multi-objectives optimization package, to determine the optimal combination of the "start of combustion (SOC)
timing" and the "exhaust gas recirculation (EGR) ratio" such that the engine would generate the highest power output
with minimal nitric oxide (NO) emission and no engine knock. The engine power output in this example
is measured by the indicated mean effective pressure (IMEP) between the intake valve close (IVC)
and the exhaust valve open (EVO). The NO emission is indicated by the cylinder-averaged peak NO mass fraction
during the same time period. The engine knock is defined as the occurrence of auto-ignition of the "endgas"
which is the gas mixture of the "unburned zone" of the SI engine model. The intensity of the engine knock in this
example is represented by the peak value of the total thermicity of the endgas. When the peak thermicity value
of a particular engine run goes above the preset threshold, the occurrence of engine knock is predicted. Consequently, the
operating condition and the solution associated with the run will be discarded because of the violation of the no engine knock
constraint enforced by the optimization process.

The multi-objective algorithm ``NGSA-II`` is employed to perform the optimization. The best pair of SOC and EGR ratio
that satisfies the no engine knock constraint will be "determined" by applying a "decision-making" method from
the *pymoo* package. The parameters and the solutions of cases around the Pareto frontier of the optimization process
will be plotted in the "design space" and the "objective space", respectively, with the best parameter/solution marked.

Details of the Chemkin 0-D SI engine model are available at the "Ansys Product Help" site. Information about
*pymoo* can be found online at pymoo.org.
"""

# sphinx_gallery_thumbnail_path = '_static/plot_SI_engine_optimization_objspace.png'

################################################
# Import PyChemkin packages and start the logger
# ==============================================

import datetime
import multiprocessing
import os

import ansys.chemkin as ck  # Chemkin
from ansys.chemkin.inlet import Stream  # external gaseous inlet
from ansys.chemkin.logger import logger
from ansys.chemkin.utilities import workingFolders

# chemkin spark ignition (SI) engine model (transient)
from ansys.chemkin.engines.SI import SIengine
import matplotlib.pyplot as plt  # plotting
import numpy as np  # number crunching
# pymoo: multi-objective optimization in Python
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.problem import StarmapParallelization
from pymoo.decomposition.asf import ASF
from pymoo.operators.sampling.lhs import LHS
from pymoo.optimize import minimize

# check working directory
current_dir = os.getcwd()
# create a separate working folder for the optimization runs
top_level = workingFolders("SI_optimize", current_dir)
top_dir = top_level.work
top_level.done()
#
logger.debug("working directory: " + os.getcwd())
# set verbose mode
ck.set_verbose(False)
# set interactive mode for plotting the results
# interactive = True: display plot
# interactive = False: save plot as a PNG file
global interactive
interactive = False


#######################################
# Create an SI engine calculator class
# =====================================
# Create a local class that wraps around the actual ``SIengine`` class
# to make the setup of the SI engine parameter study more convenient.
#
class SIengineCalculator:
    """
    SI engine calculator with fixed set up parameters
    """
    def __init__(self, fresh_mixture: Stream):
        """
        SI engine calculator that instantiates an SI engine object with
        the given fresh (unburnt) mixture condition and predefined engine
        parameters.

        Parameters
        ----------
            fresh_mixture: Mixture object
                the initial/fresh/unburnt condition
        """
        # instantiate an SI engine object
        # set up the run and working direcotry name
        self.name = "SI_engine"
        self.SIcalculator = SIengine(fresh_mixture, self.name)
        # run status
        self.runstatus = -100
        # IMEP (objective #1)
        self.IMEP = 0.0
        # peak NO mass fraction (objective #2)
        self.max_NO = 1.0
        # calculated maximum thermicity value in the unburned zone [1/sec] (objective #3)
        self.max_sigma = 1.0e5
        # peak endgas temperature [K]
        self.max_temp = 0.0
        # set the required premixed flame model parameters
        self.set_engine_parameters()
        # set burn profile (variables/parameters)
        # the burn profile will be set right before each run during the optimization process

    def run(self) -> int:
        """
        Run the SI engine calculation

        Returns
        -------
            status: integer
                run status code
        """
        # reset status
        self.runstatus = -100
        # run the flame speed calculation
        self.runstatus = self.SIcalculator.run()
        # extract the laminar flame speed from the solution
        if self.runstatus == 0:
            # postprocess the unburned zone solution
            unburnedzone = 1
            burnedzone = 2
            # because the memory is shared, it must be done as soon as the run is finished
            # get solution variables from the unburned zone
            iErr = self.SIcalculator.process_engine_solution(zoneID=unburnedzone)
            if iErr == 0:
                temperature = self.SIcalculator.get_solution_variable_profile("temperature")
                thermicity = self.SIcalculator.get_solution_variable_profile("thermicity")
                self.IMEP = self.SIcalculator.get_engine_IMEP()
                self.max_temp = np.max(temperature)
                self.max_sigma = np.max(thermicity)
                # get the cylinder-averaged solution variables
                iErr_a = self.SIcalculator.process_average_engine_solution()
                NO_mass_frac = self.SIcalculator.get_solution_variable_profile("no")
                self.max_NO = np.max(NO_mass_frac)
                del thermicity, temperature, NO_mass_frac
            else:
                self.runstatus = iErr
        return self.runstatus

    def set_engine_parameters(self):
        """
        Set up basic engine parameters (key geometries and operating conditions)
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
        self.SIcalculator.RPM = 900
        # simulation start CA [degree]
        self.SIcalculator.starting_CA = -120.2
        # simulation end CA [degree]
        self.SIcalculator.ending_CA = 139.8
        # wall heat transfer parameters
        # Set up engine wall heat transfer model
        heattransferparameters = [0.15, 0.8, 0.0]
        # set cylinder wall temperature [K]
        Twall = 434.0
        self.SIcalculator.set_wall_heat_transfer("dimensionless", heattransferparameters, Twall)
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
        self.SIcalculator.tolerances = (1.0e-12, 1.0e-6)
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
        # self.SIcalculator.set_ignition_delay(method="T_inflection")
        # set species lower bound for equilibrium calculation
        self.SIcalculator.set_burned_products_minimum_mole_fraction(1.0e-12)
        # suppress text output from the SI engine model
        self.SIcalculator.stop_output()

    def set_burn_profile(self, start_combustion: float, burn_duration: float, wiebe_n: float, wiebe_b: float):
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

    def get_solution_parameters(self) -> tuple[float, float, float]:
        """
        Get the SI engire solution varibles to be optimized

        Returns
        -------
            IMEP: double
                indicated mean effective pressure [bar]
            max_NO: double
                peak cylinder-averaged NO mass fraction [-]
            max_sigma: double
                peak total thermicity of the eng gas [1/sec]
        """
        return self.IMEP, self.max_NO, self.max_sigma



##########################################
# Create an instance of the Chemistry Set
# ========================================
# For PRF, the encrypted 14-component gasoline mechanism, ``gasoline_14comp_WBencrypted.inp``,
# is used. The chemistry set is named ``gasoline``.

def setup_fresh_mixture(fuel_compo: list[tuple[str, float]], phi: float, egr_rate: float) -> Stream:
    """
    Set up the chemistry set and create the fresh fuel-air mixture at the intake valve closing (IVC).

    Parameters
    ----------
        fuel: tuple (str, double)
            fuel species mole fractions
        phi: double
            equivalence ratio
        egr_rate: double
            EGR rate ratio

    Returns
    -------
        fresh: Mixture object
            fresh unburned mixture (initial charge)
    """
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
#
# .. note::
#   The "surrogate blend utility" from the Chemkin "reaction workbench" can be
#   used to create a surrogate that matches the properties (heating value, distillation curve,
#   density, viscosity, etc.) of an actual fuel.
#

    # create the fuel mixture
    fuelmixture = ck.Mixture(MyGasMech)
    # set fuel composition
    fuelmixture.X = fuel_compo
    # setting pressure and temperature is not required in this case
    fuelmixture.pressure = 1.15 * ck.Patm
    # Update the fresh mixture temperature according to the EGR ratio
    # intake gas temperature [K]
    temp_intake = 310.0
    # engine exhaust gas temperature [K]
    temp_egr = 950.0
    # gas charge temperature [K]
    fuelmixture.temperature = (1.0 - egr_rate) * temp_intake + egr_rate * temp_egr

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
    equiv = phi
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
    fresh.pressure = fuelmixture.pressure
    fresh.temperature = fuelmixture.temperature

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
    EGRratio = egr_rate
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

    return fresh


##################################
# Set up the optimization problem
# ================================
# Define the optimization problem. Set up the optimization parameters, variables, objectives, constraints, and
# the bounds of the variables.
#
# Since the ``_evaluate()`` method perform one SI engine simulation at a time, the optimization problem is defined as an
# ``ElementwiseProblem``. By implementing the SI engine model as an element-wise problem, the parallelization available in
# *pymoo* can be directly utilized.
#
# There are two design parameters/variables in this project. The variable SOC is used to set up the mass burned fraction profile
# of the SI engine model. And the EGR ratio is applied to create the initial gas charge inside the engine cylinder. The objective and
# the constraint values come from the result of the SI engine simulation. The first objective is the negative value of IMEP,
# the second objective is the peak NO mass fraction in the cylinder, and the constraint value is the knock intensity
# as indicated by the total thermicity of the endgas.
# 
# In *pymoo*, the objectives are "minimized" so the objective of maximizing the IMEP must be reformulated into minimizing the negative value
# of IMEP, and the constraint must be written in the form of "g <= 0.0".
#
# .. note::
#   See the *pymoo* documentation for details about how the optimization problem is set up and evaluated.
#
class SI_engine_opt(ElementwiseProblem):
    def __init__(self, numb_vars: int, numb_objs: int, fuel: list[tuple[str, float]], phi: float, **kwargs):
        """
        Set up the SI engine optimization.
        The objective is to find the optimal fuel mass burn profile to
        maximize the engine indicated mean effective pressure (IMEP) while
        avoiding engine knock and minimizing the NO emission.
        The fuel mass burn profile is described by the anchor points:
        crank angles at 10%, 50%, and 90% mass burned.

        Parameters
        ----------
            numb_vars: integer
                number of variables
            numb_objs: integer
                number of objectives
            fuel: tuple (str, double)
                fuel species mole fractions
            phi: double
                equivalence ratio of the initial gas charge
        """
        # number of inequality constraints
        numb_constrs = 1
        # lower bounds for the variables
        lower_bounds = np.zeros (numb_vars, dtype=np.double)
        upper_bounds = np.zeros (numb_vars, dtype=np.double)
        lower_bounds[0] = -15.0
        lower_bounds[1] = 0.0
        # upper bounds for the variables
        upper_bounds[0] = 10.0
        upper_bounds[1] = 0.60
        # initializa the "Probelm" object
        super().__init__(n_var=numb_vars,
                         n_obj=numb_objs,
                         n_ieq_constr=numb_constrs,
                         n_eq_constr=0,
                         xl=lower_bounds,
                         xu=upper_bounds,
                         )
        # number of evaluation calls
        self.call_count = 0
        self.fail_count = 0
        # initial gas charge parameters
        # fuel composition
        self.fuel = fuel
        # equivalence ratio
        self.phi = phi

    def _evaluate(self, x, out, *args, **kwargs):
        """
        Evaluate the objective function values with the given
        set of variable values
    
        Parameters
        ----------
            x: list[float]
                values of the pymoo variables to be optimized.
            out: dict[str: list[float]]
               pymoo OUT dictionary containing evaluated function values and constrain values
        """
        # initialize "bad" return values to drive the search away from this state point
        f1 = 0.0
        f2 = 1.0
        # 
        g1 = 5.0e4
        # variables
        # x[0]: start of combustion crank angle
        # x[1]: exhaust gas recirculation (EGR) rate
        # create and change to the working directory for this run
        name = "SI_engine_" + str(self.call_count)
        sub_work_folder = workingFolders(name, top_dir)
        # set up the unburned mixture of the SI engine at IVC
        init_mixture = setup_fresh_mixture(self.fuel, self.phi, x[1])
        # instatiate the SI engine
        engine_case = SIengineCalculator(init_mixture)
        # set the burned mass fraction profile
        # use two sets of Wiebe function parameters to describe the burned mass fraction profile
        # in the SI engine
        # profile #1: fast burn
        wiebe_n_egr = 1.8
        wiebe_b_egr = 8.5
        # profile #2: regular burn
        wiebe_n_fresh = 2.0
        wiebe_b_fresh = 8.0
        # EGR ratio
        r = x[1]
        r1 = 1.0 - r
        # burned mass profile as a function of EGR ratio
        wiebe_b = r * wiebe_b_egr + r1 * wiebe_b_fresh
        wiebe_n = r * wiebe_n_egr + r1 * wiebe_n_fresh
        burn_duration = 45.5  # degrees
        # set up the mass burned profile in the SI engine
        # the burn duration is fixed
        engine_case.set_burn_profile(x[0], burn_duration, wiebe_n, wiebe_b)
        # objects
        # f1: IMEP [bar] (to be maximized)
        # f2: peak cylinder-averaged NO mass fraction [-] (to be minimized)
        runstatus = engine_case.run()
        if runstatus == 0:
            # get the results only when the simulation is successful
            f1, f2, max_sigma = engine_case.get_solution_parameters()
            # compute the constraints
            # g1: minimal engine knock
            # max_sigma: peak endgas total thermicity (knock indicator) [1/sec]
            g1 = max_sigma - 5.0e4
        else:
            # tally the number of failed runs
            self.fail_count += 1
        # clean up
        del engine_case
        sub_work_folder.done()
        #
        self.call_count += 1
        # return solution values to the optimizer
        # return negative solution value for maximization objective
        out["F"] = [-f1, f2]
        # return constraints
        out["G"] = [g1]


######################
# Run the optimization
# ====================
# Select the optimization algorithm and the decision making method. Since this project
# has two objectives and one constraint, the multi-objective algorithm ``NSGA-II`` is selected
# as the optimization algorithm. A relatively small population size of 50 is used for this simple
# optimization problem.
#
# The ``StarmapParallelization()`` method from *pymoo* is used to parallelize the
# SI engine simulation runs during the optimization process. 
#
def run_optimization(numb_vars: int, numb_objs: int, fuel: list[tuple[str, float]], phi: float):
    """
    Runing the optimization project

    Parameters
     ----------
        numb_vars: integer
            number of variables
        numb_objs: integer
            number of objectives
        fuel: tuple (str, double)
            fuel species mole fractions
        phi: double
            equivalence ratio of the initial gas charge

    Returns
    -------
        res: pymoo Result object
            pymoo minimization results and history
    """
    # create the "problem" to be optimized
    # set up the optimization algorithm
    # NSGA: Non-dominated Sorting Genetic Algorithm
    # LHS: Latin Hypercube Sampling
    algorithm = NSGA2(
                    pop_size=50,
                    sampling=LHS(),
                    eliminate_duplications=True,
                )
    # set up the thermination condition
    # run the optimization
    run_parallel = True
    if run_parallel:
        # parallel mode
        # initialize the process pool and create the runner
        # number of available cpu cores
        numb_cores = max(os.cpu_count(), 1)
        numb_workers = numb_cores - 2
        numb_workers = 4  # max(numb_workers, 1)
        pool = multiprocessing.Pool(numb_workers)
        runner = StarmapParallelization(pool.starmap)
        # define the problem by passing the starmap interface of the thread pool
        problem = SI_engine_opt(numb_vars=numb_vars, numb_objs=numb_objs, fuel=fuel, phi=phi, elementwise_runner=runner)
    else:
        # series mode
        problem = SI_engine_opt(numb_vars=numb_vars, numb_objs=numb_objs, fuel=fuel, phi=phi)

    res = minimize(problem, algorithm, termination=("n_gen", 1), seed=1)
    # run summary
    # execution time [sec]
    print(f"Total wall time: {datetime.timedelta(seconds=res.exec_time)}")
    # run statistics
    print(f"Failed runs: {problem.fail_count}/{problem.call_count}")
############################
# Post-process the solutions
# ==========================
# Parse the solutions from the optimization process and display them with the best solution
# marked by a green cross mark. The variable values are stored in the ``Result.X`` array and the
# corresponding objective values are in the ``Result.F`` array. Because the objective of maximizing
# the "IMEP" has to be converted to a minimization objective in *pymoo*, an alternative objective
# of minimizing "-IMEP" is adopted. Consequently, the raw -IMEP objective solution must be converted back
# to IMEP to make sense of the optimization and the subsequent "decision-making" results.
#
# In this project, the decomposition function ``ASF()`` is applied to retrieve the best solution. Weight factors
# are used to make the IMEP objective slightly more important than the NO emission objective.

    # process the results
    X = res.X
    F = res.F
    # convert the objective result back to the original IMEP
    IMEP = -1.0 * F[:, 0]
    #
    approx_ideal = F.min(axis=0)
    approx_nadir = F.max(axis=0)
    # finding the optimal result
    # normalize the results
    diff = approx_ideal - approx_nadir
    if np.isclose(diff[0], 0.0, atol=1.0e-9) or np.isclose(diff[1], 0.0, atol=1.0e-9):
        nF = F - approx_ideal
    else:
        nF = (F - approx_ideal) / (approx_ideal - approx_nadir)
    # set the weight factors
    weight = np.zeros(n_objs, dtype=np.double)
    weight[0] = 1.0 / 0.6
    weight[1] = 1.0 / 0.4
    # find the best solution
    # use the Augmented Scalarization Function (ASF) method
    decomp = ASF()
    i = decomp.do(nF, weight).argmin()
    print("Best outcome according to ASF:")
    print(f"run # = {i}")
    print(f"Predicted IMEP = {IMEP[i]} [bar]")
    print(f"Predicted peak NO mass fraction = {F[i, 1]} [-]")
    print(f"Start of combustion = {X[i, 0]} [degree CA ATDC]")
    print(f"EGR ratio = {X[i, 1] * 100.0} % vol.")

    # plot the design space
    xl, xu = problem.bounds()
    plt.figure(figsize=(7, 5))
    plt.title("Design Space")
    plt.scatter(X[:, 0], X[:, 1], s=30, facecolors='none', edgecolors='b')
    plt.scatter(X[i, 0], X[i, 1], color='green', marker="x", s=75, label="Best")
    plt.xlim(xl[0], xu[0])
    plt.ylim(xl[1], xu[1])
    plt.ylabel("EGR ratio [-]")
    plt.xlabel("Start of Combustion [degree CA]")
    plt.legend()
    if interactive:
        plt.show()
    else:
        plt.savefig("plot_SI_engine_optimization_designspace.png",bbox_inches="tight")

    # plot the objective space to show the Pareto frontier
    plt.figure(figsize=(7, 5))
    plt.scatter(IMEP, F[:, 1], s=30, facecolors='none', edgecolors='b')
    # note that because IMEP = -F[:, 0], approx_ideal[0] and approx_nadir[0] must be modified accordingly
    ideal_IMEP = -1.0 * approx_ideal[0]
    nadir_IMEP = -1.0 * approx_nadir[0]
    plt.scatter(ideal_IMEP, approx_ideal[1], facecolors='none', edgecolors='red', marker="*", s=100, label="Ideal Point (Approx)")
    plt.scatter(nadir_IMEP, approx_nadir[1], facecolors='none', edgecolors='black', marker="p", s=100, label="Nadir Point (Approx)")
    plt.scatter(IMEP[i], F[i, 1], color='green', marker="x", s=75, label="Best")
    plt.title("Objective Space - Power & Emission")
    plt.xlabel("IMEP [bar]")
    plt.ylabel("Peak NO mass fraction [-]")
    plt.legend()
    if interactive:
        plt.show()
    else:
        plt.savefig("plot_SI_engine_optimization_objspace.png",bbox_inches="tight")
    #
    del IMEP, nF, xl, xu
    #
    return res

####################################################
# Set up the main driver of the optimization process
# ==================================================
# A main method for the optimization process is necessary especially with the intension to run the
# design point cases in parallel.
# 
# You can further post-process the results from the optimization process. Here the correlations between
# the design variables (SOC and EGR ratio) and the objectives (IMEP and NO emission)
# are explored. You can consult *pymoo* online API document to find out additional information offered by
# the ``Result`` object. 
#
if __name__ == "__main__":
    # properties of the initial gas charge
    # fuel composition
    # the "surrogate blend utility" from the Chemkin "reaction workbench" can be
    # used to create a surrogate that matches the properties (heating value, distillation curve,
    # density, viscocity, etc.) of an actual fuel
    fuel = [("ic8h18", 0.6), ("nc7h16", 0.0), ("mch", 0.1), ("c6h5c3h7", 0.3)]
    # equivalence ratio
    phi = 1.0
    # perform optimization
    # number of variables
    n_vars = 2
    # number of objectives
    n_objs = 2
    opt_results = run_optimization(n_vars, n_objs, fuel, phi)
    # process the optimization results
    X = opt_results.X
    F = opt_results.F
    IMEP = -1.0 * F[:, 0]
    # plot correlations between the variable/parameter and the objectives
    plt.figure(figsize=(7, 5))
    plt.scatter(X[:, 0], IMEP, s=30, facecolors='none', edgecolors='r')
    plt.title("Correlation - SOC & Power")
    plt.ylabel("IMEP [bar]")
    plt.xlabel("Start of combustion [CA]")
    if interactive:
        plt.show()
    else:
        plt.savefig("plot_SI_engine_optimization_correlation1.png",bbox_inches="tight")
    #
    plt.figure(figsize=(7, 5))
    plt.scatter(X[:, 1], IMEP, s=30, facecolors='none', edgecolors='r')
    plt.title("Correlation - EGR & Power")
    plt.ylabel("IMEP [bar]")
    plt.xlabel("EGR ratio [-]")
    if interactive:
        plt.show()
    else:
        plt.savefig("plot_SI_engine_optimization_correlation2.png",bbox_inches="tight")

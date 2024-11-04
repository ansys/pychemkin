from ctypes import c_double, c_int

import numpy as np

from chemkin.color import Color as Color
from chemkin.inlet import Inlet
from chemkin.mixture import Mixture
from chemkin.reactormodel import Keyword, ReactorModel


class openreactor(ReactorModel):
    """
    Generic open reactor unit
    """

    def __init__(self, guessedmixture, label=None):
        """
        Create a steady-state flow reactor
        :param guessedmixture: guessed/estimate reactor condition (Mixture object)
        :param label: inlet name (string)
        """
        # check Inlet
        if isinstance(guessedmixture, Mixture):
            # initialization
            super().__init__(guessedmixture, label)
        else:
            # wrong argument type
            print(
                Color.RED + "** the first argument must be a Mixture object",
                end=Color.END,
            )
            raise TypeError
        if label is None:
            self.label = "Reactor"
        else:
            self.label = label
        # use API mode for steady-state open reactor/flame simulations
        Keyword.noFullKeyword = True
        # FORTRAN file unit of the text output file
        self._myLOUT = c_int(158)
        # inlet inforamtion
        # number of external inlets
        self.numbexternalinlets = 0
        # dict of external inlet objects {inlet label: inlet object}
        self.externalinlets = {}
        # total mass flow rate into this reactor [g/sec]
        self.totalmassflowrate = 0.0
        #
        self.SolverTypes = {"Transient": 1, "SteadyState": 2}
        self.EnergyTypes = {"ENERGY": 1, "GivenT": 2}
        # specify "reactor residence time" = "CONP" or
        # specify "reactor volume" = "CONV"
        self.ProblemTypes = {"CONP": 1, "CONV": 2}
        #
        # initialize the steady-state solver control parameter with the detault setting
        # >>> steady-state search algorithm:
        # absolute tolerance for the steady-state solution
        self.SSabsolutetolerance = 1.0e-9
        # relative tolerance for the steady-state solution
        self.SSrelativetolerance = 1.0e-4
        # max number of iterations per steady state search
        self.SSmaxiteration = 100
        # number of steady-state searches before evaluating new Jacobian matrix
        self.SSJacobianage = 20
        # max number of calls to pseudo transient algorithm
        self.maxpseudotransient = 100
        # number of pseudo transient "steps" before calling the steady-state search algorithm
        self.numbinitialpseudosteps = 0
        # upper bound of the temperature value during iteration
        self.maxTbound = 5000.0  # [K]
        # floor value (lower bound) of the gas species mass fraction during iteration
        self.speciesfloor = -1.0e-14
        # reset negative gas species fraction to the given value in intermediate solution
        self.speciespositive = None
        # use legacy steady-state solver algorithm
        self.uselegacytechnique = None
        # use damping in search: 0 = OFF; 1 = ON
        self.SSdamping = 1
        # absolute perturbation for Jacobian evaluation
        self.absoluteperturbation = None
        # relative perturbation for Jacobian evaluation
        self.relativeperturbation = None
        # >>> pseudo trasient (time stepping) algorithm:
        # absolute tolerance for the time stepping solution
        self.TRabsolutetolerance = 1.0e-9
        # relative tolerance for the time stepping solution
        self.TRrelativetolerance = 1.0e-4
        # max number of iterations per pseudo time step before cutting the time step size
        self.TRmaxiteration = 25
        # max number of pseudo time steps before increasing the time step size
        self.timestepsizeage = 25
        # minimum time step size allowed
        self.TRminstepsize = 1.0e-10  # [sec]
        # maximum time step size allowed
        self.TRmaxstepsize = 1.0e-2  # [sec]
        # time step size increasing factor
        self.TRupfactor = 2.0
        # time step size decreasing factor
        self.TRdownfactor = 2.2
        # number of pseudo time steps before evaluating new Jacobian matrix
        self.TRJacobianage = 20
        # initial stride and number of steps per pseudo time stepping call
        # for fixed-temperature solution
        self.TRstride_fixT = 1.0e-6  # [sec]
        self.TRnumbsteps_fixT = 100
        # for energy equation solution
        self.TRstride_ENRG = 1.0e-6  # [sec]
        self.TRnumbsteps_ENRG = 100
        # solver message output level: 0 ~ 2
        self.printlevel = 1
        #
        # steady-state solver keywords
        self.SSsolverkeywords = {}

    def setinlet(self, extinlet):
        """
        Add an external inlet to the reactor
        :param extinlet: external inlet (Inlet object)
        """
        # check Inlet
        if not isinstance(extinlet, Inlet):
            # wrong argument type
            print(
                Color.RED + "** the argument must be an Inlet object",
                end=Color.END,
            )
            raise TypeError
        # current external inlet count
        count = self.numbexternalinlets + 1
        if extinlet.label is None:
            inletname = self.label + "_inlet_" + str(count)
        else:
            # inlet has label
            inletname = self.label + "_" + extinlet.label
        # check inlet name uniqueness
        if inletname in self.externalinlets:
            print(Color.YELLOW + f"** inlet {inletname} already connected")
            print("   will append '_dup' to the original name", end=Color.END)
            inletname += "_dup"
        # check inlet flow rate
        if extinlet._flowratemode < 0:
            # no given in the inlet
            print(Color.PURPLE + "** inlet flow rate is not set")
            print("   specify flow rate of the 'Inlet' object", end=Color.END)
            exit()
        else:
            flowrate = extinlet.massflowrate
            if flowrate <= 0.0:
                print(Color.PURPLE + "** inlet flow rate is not set correctly")
                print("   specify flow rate of the 'Inlet' object", end=Color.END)
                exit()
        # add the inlet object to the inlet dict of the reactor
        self.externalinlets[inletname] = extinlet
        self.numbexternalinlets = count
        self.totalmassflowrate += flowrate
        print(Color.YELLOW + f"** new inlet {inletname} is added", end=Color.END)

    def removeinlet(self, name):
        """
        Delete an existing external inlet from the reactor by the inlet name
        :param name: external inlet name (string)
        """
        # check input
        if not isinstance(name, str):
            print(Color.PURPLE + "** the argument must be a string", end=Color.END)
            exit()
        # delete the named inlet from the externalinlets dict
        missed = False
        if name in self.externalinlets:
            # existing inlet
            extinlet = self.externalinlets.pop(name, None)
            if extinlet is None:
                missed = True
            elif not isinstance(extinlet, Inlet):
                # some internal messed up
                missed = True
                print(Color.RED + "** failed to found the inlet data", end=Color.END)
                raise TypeError
            else:
                # decrease the external inlet count by 1
                self.numbexternalinlets -= 1
                # take out the mass flow rate contribution
                self.totalmassflowrate -= extinlet.massflowrate
                print(Color.YELLOW
                      + f"** inlet {name} is removed from reactor {self.label}",
                      end=Color.END)
        else:
            # not in the external inlet dict
            missed = True
        if missed:
            print(Color.PURPLE + f"** inlet {name} not found", end=Color.END)

    #
    # solver parameters keyword processing
    @property
    def absolutetolerance(self):
        """
        absolute tolerance for the steady-state solution
        """
        return self.SSabsolutetolerance

    @property
    def relativetolerance(self):
        """
        relative tolerance for the steady-state solution
        """
        return self.SSrelativetolerance

    def settolerances(self, atol, rtol):
        """
        set the absolute and the relative tolerances for the steady-state solution
        :param atol: absolutie tolerance (double scalar)
        "param rtol: relative tolerance (double scalar)
        """
        iErr = 0
        if atol > 0.0:
            self.SSsolverkeywords["ATOL"] = atol
            self.SSabsolutetolerance = atol
        else:
            iErr = 1

        if rtol > 0.0:
            self.SSsolverkeywords["RTOL"] = rtol
            self.SSrelativetolerance = rtol
        else:
            iErr = 1

        if iErr > 0:
            print(Color.PURPLE + "** tolerance must > 0", end=Color.END)

    def settimesteppingtolerances(self, atol, rtol):
        """
        set the absolute and the relative tolerances for the pseudo time stepping solution
        :param atol: absolutie tolerance (double scalar)
        "param rtol: relative tolerance (double scalar)
        """
        iErr = 0
        if atol > 0.0:
            self.SSsolverkeywords["ATIM"] = atol
            self.TRabsolutetolerance = atol
        else:
            iErr = 1

        if rtol > 0.0:
            self.SSsolverkeywords["RTIM"] = rtol
            self.TRrelativetolerance = rtol
        else:
            iErr = 1

        if iErr > 0:
            print(Color.PURPLE + "** tolerance must > 0", end=Color.END)

    def setmaxpseudotransientcall(self, maxtime):
        """
        set the maximum number of call to the pseudo transient algorithm
        in an attempt to find the steady-state solution
        :param maxtime: max number of pseudo transient calls/attempts (integer scalar)
        """
        if maxtime >= 1:
            self.SSsolverkeywords["MAXTIME"] = maxtime
            self.maxpseudotransient = maxtime
        else:
            print(Color.PURPLE + "** parameter must > 0", end=Color.END)

    def setmaxnumbtimestepiteration(self, maxiteration):
        """
        set the maximum number of iterations per time step when performing the pseudo transient algorithm
        :param maxtime: max number of iterations per pseudo time step (integer scalar)
        """
        if maxiteration >= 1:
            self.SSsolverkeywords["TRMAXITER"] = maxiteration
            self.TRmaxiteration = maxiteration
        else:
            print(Color.PURPLE + "** parameter must > 0", end=Color.END)

    def setmaxnumbsearchiteration(self, maxiteration):
        """
        set the maximum number of iterations per search when performing the steady-state search algorithm
        :param maxtime: max number of iterations per steady-state search (integer scalar)
        """
        if maxiteration >= 1:
            self.SSsolverkeywords["SSMAXITER"] = maxiteration
            self.SSmaxiteration = maxiteration
        else:
            print(Color.PURPLE + "** parameter must > 0", end=Color.END)

    def setinitialtimesteps(self, initsteps):
        """
        set the number of pseudo time steps to be performed to establish a "better"
        set of guessed solution before start the actual steady-state solution search
        :param initsteps: number of initial pseudo time steps (integer scalar)
        """
        if initsteps >= 1:
            self.SSsolverkeywords["ISTP"] = initsteps
            self.numbinitialpseudosteps = initsteps
        else:
            print(Color.PURPLE + "** parameter must > 0", end=Color.END)

    def setspeciesfloor(self, floorvalue):
        """
        set the minimum species fraction value allowed during steady-state solution search
        :param floorvalue: minimum species fraction value (double scalar)
        """
        if np.abs(floorvalue) < 1.0:
            self.SSsolverkeywords["SFLR"] = floorvalue
            self.speciesfloor = floorvalue
        else:
            print(Color.PURPLE + "** species floor value must < 1", end=Color.END)

    def settemperatureceiling(self, ceilingvalue):
        """
        set the maximum temperature value allowed during steady-state solution search
        :param ceilingvalue: maximum temperature value (double scalar)
        """
        if ceilingvalue > 300.0:
            self.SSsolverkeywords["TBND"] = ceilingvalue
            self.maxTbound = ceilingvalue
        else:
            print(Color.PURPLE + "** temperature value must > 300", end=Color.END)

    def setspeciesresetvalue(self, resetvalue):
        """
        set the positive value to reset any negative species fraction in
        intermediate solutions during iterations
        :param resetvalue: positive value to reset negative species fraction (double scalar)
        """
        if resetvalue >= 0.0:
            self.SSsolverkeywords["SPOS"] = resetvalue
            self.speciespositive = resetvalue
        else:
            print(Color.PURPLE + "** species fraction value must >= 0", end=Color.END)

    def setmaxpseudotimestepsize(self, dtmax):
        """
        set the maximum time step sizes allowed by the pseudo time stepping solution
        :param dtmax: maximum time step size allowed (double scalar)
        """
        if dtmax > 0.0:
            self.SSsolverkeywords["DTMX"] = dtmax
            self.TRmaxstepsize = dtmax
        else:
            print(Color.PURPLE + "** time step size must > 0", end=Color.END)

    def setminpseudotimestepsize(self, dtmin):
        """
        set the minimum time step size allowed by the pseudo time stepping solution
        "param dtmin: minimum time step size allowed (double scalar)
        """
        if dtmin > 0.0:
            self.SSsolverkeywords["DTMN"] = dtmin
            self.TRminstepsize = dtmin
        else:
            print(Color.PURPLE + "** time step size must > 0", end=Color.END)

    def setpseudotimestepage(self, age):
        """
        set the minimum number of time steps taken before allowing time step size increase
        "param age: min age of the pseudo time step size (integer scalar)
        """
        if age > 0:
            self.SSsolverkeywords["IRET"] = age
            self.timestepsizeage = age
        else:
            print(Color.PURPLE + "** number of time step must > 0", end=Color.END)

    def setJacobianage(self, age):
        """
        set the number of steady-state searches before re-evaluate the Jacobian matrix
        "param age: age of the steady-state Jacobian matrix (integer scalar)
        """
        if age > 0:
            self.SSsolverkeywords["NJAC"] = age
            self.SSJacobianage = age
        else:
            print(Color.PURPLE + "** number of time step must > 0", end=Color.END)

    def setpseudoJacobianage(self, age):
        """
        set the number of time steps taken before re-evaluate the Jacobian matrix
        "param age: age of the pseudo time step Jacobian matrix (integer scalar)
        """
        if age > 0:
            self.SSsolverkeywords["TJAC"] = age
            self.TRJacobianage = age
        else:
            print(Color.PURPLE + "** number of time step must > 0", end=Color.END)

    def setdampingoption(self, ON):
        """
        turn ON (True) or OFF (False) the damping option of the steady-state solver
        :param ON: turn On the damping option (boolean scalar)
        """
        if type(ON) == bool:
            if ON:
                self.SSdamping = 1
            else:
                self.SSdamping = 0
            self.SSsolverkeywords["TWOPNT_DAMPING_OPTIN"] = self.SSdamping
        else:
            print(Color.PURPLE + "** parameter must be True or False", end=Color.END)

    def setlegacyoption(self, ON):
        """
        turn ON (True) or OFF (False) the legacy steady-state solver
        :param ON: turn On the legacy solver (boolean scalar)
        """
        if type(ON) == bool:
            self.uselegacytechnique = ON
            if ON:
                self.SSsolverkeywords["USE_LEGACY_TECHNIQUE"] = "4X"
        else:
            print(Color.PURPLE + "** parameter must be True or False", end=Color.END)

    def setprintlevel(self, level):
        """
        set the level of information to be provided by the steady-state solver
        to the text output
        :param level: solver message details level (0 ~ 2) (integer scalar)
        """
        if level in [0, 1, 2]:
            self.SSsolverkeywords["PRNT"] = level
            self.printlevel = level
        else:
            print(Color.PURPLE + "** print level must be 0, 1, or 2", end=Color.END)

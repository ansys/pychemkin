# Steady-State Solver Control Parameters
#
import numpy as np

from .color import Color


class steadystatesolver:
    """
    Common steady-state solver controlling parameters
    """

    def __init__(self):
        # steady-state solver control parameter class
        # mostly just keyword processing
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

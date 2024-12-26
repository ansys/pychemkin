# The "Inlet" class is an extension of the "Mixture" class
#
import copy

from .chemistry import Patm
from .color import Color
from .mixture import Mixture


class Inlet(Mixture):
    """
    define an inlet based on the gas species in the given chemistry set for open reactor models
    """

    def __init__(self, chem):
        """
        set up an inlet with a given chemistry set for open reactor models
        :param chem: Chemistry object
        """
        super().__init__(chem)
        # 0=mass flow rate/1=volumetric flow rate/2=velocity/3=SCCM
        # flag
        self._flowratemode = -1  # not given
        self._inletflowrate = [0.0] * 4
        # types of flow rate allowed
        self._massflowrate = 0.0  # mass flow rate FLRT [g/sec]
        self._volflowrate = 0.0  # volumetric flow rate VDOT [cm3/sec]
        self._velocity = 0.0  # gas velocity VEL [cm/sec] for plug flow reactor model
        self._SCCM = 0.0  # standard (198.15K, 1atm) cubic centimeters per minute SCCM [standard cm3/min]
        # inlet velocity gradient [1/sec] (for premixed, oppdif, amd spin)
        self._velgrad = 0.0
        # flow area (for velocity in plug flow reactor model)
        self._haveflowarea = False
        # cross-sectional flow area [cm2]
        self._flowarea = 1.0

    def converttomassflowrate(self):
        """
        convert different types of flow rate value to mass flow rate
        :return: mass flow rate [g/sec] (double scalar)
        """
        #
        if self._flowratemode == 1:
            # volumetric flow rate
            # get inlet gas mixture density
            mrate = self.RHO * self._volflowrate
            return mrate
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                mrate = self.RHO * self._flowarea * self._velocity
                return mrate
            else:
                # no flow area
                print(
                    Color.PURPLE + "** flow area value is not given for this inlet",
                    end=Color.END,
                )
                return 0.0

        elif self._flowratemode == 3:
            # SCCM
            chemID = self._chemset_index.value
            # set standard condition
            p = Patm  # [atm]
            t = 298.15  # [K]
            # set mass fractions
            frac = copy.deepcopy(self.Y)
            # molecular masses
            wt = self._WT
            # get gas density at the standard condition
            standard_den = Mixture.density(chemID, p, t, frac, wt, mode="mass")
            mrate = standard_den * self._SCCM / 60.0
            del frac
            return mrate
        else:
            print(Color.PURPLE + "** unknown flow rate units", end=Color.END)
            return 0.0

    def converttovolflowrate(self):
        """
        convert different types of flow rate value to volumetric flow rate
        :return: volmetric flow rate [cm3/sec] (double scalar)
        """
        #
        if self._flowratemode == 0:
            # mass flow rate
            # get inlet gas mixture density
            vrate = self._massflowrate / self.RHO
            return vrate
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                vrate = self._flowarea * self._velocity
                return vrate
            else:
                # no flow area
                print(
                    Color.PURPLE + "** flow area value is not given for this inlet",
                    end=Color.END,
                )
                return 0.0

        elif self._flowratemode == 3:
            # SCCM
            chemID = self._chemset_index.value
            # set standard condition
            p = Patm  # [atm]
            t = 298.15  # [K]
            # set mass fractions
            frac = copy.deepcopy(self.Y)
            # molecular masses
            wt = self._WT
            # get gas density at the standard condition
            standard_den = Mixture.density(chemID, p, t, frac, wt, mode="mass")
            mrate = standard_den * self._SCCM / 60.0
            vrate = mrate / self.RHO
            del frac
            return vrate
        else:
            print(Color.PURPLE + "** unknown flow rate units", end=Color.END)
            return 0.0

    def converttoSCCM(self):
        """
        convert different types of flow rate value to SCCM
        :return: SCCM [standard cm3/min] (double scalar)
        """
        #
        chemID = self._chemset_index.value
        # set standard condition
        p = Patm  # [atm]
        t = 298.15  # [K]
        # set mass fractions
        frac = copy.deepcopy(self.Y)
        # molecular masses
        wt = self._WT
        # get gas density at the standard condition
        standard_den = Mixture.density(chemID, p, t, frac, wt, mode="mass")
        del frac
        #
        if self._flowratemode == 0:
            # mass flow rate
            sccm = self._massflowrate / standard_den * 60.0
            return sccm
        elif self._flowratemode == 1:
            # volumetric flow rate
            # get inlet gas mixture density
            mrate = self.RHO * self._volflowrate
            sccm = mrate / standard_den * 60.0
            return sccm
        elif self._flowratemode == 2:
            # velocity
            if self._haveflowarea:
                mrate = self.RHO * self._flowarea * self._velocity
                sccm = mrate / standard_den * 60.0
                return sccm
            else:
                # no flow area
                print(
                    Color.PURPLE + "** flow area value is not given for this inlet",
                    end=Color.END,
                )
                return 0.0
        else:
            print(Color.PURPLE + "** unknown flow rate units", end=Color.END)
            return 0.0

    @property
    def flowarea(self):
        """
        Get inlet flow area
        :return: flow cross-sectional area [cm2] (double scalar)
        """
        if self._haveflowarea:
            return self._flowarea
        else:
            print(
                Color.PURPLE + "** flow area value is not given for this inlet",
                end=Color.END,
            )
            return self._flowarea

    @flowarea.setter
    def flowarea(self, farea):
        """
        Set inlet cross-sectional flow area
        :param farea: flow area [cm2] (double scalar)
        :return: None
        """
        if farea <= 0.0:
            print(Color.RED + "** invalid flow area value", end=Color.END)
            return
        self._haveflowarea = True
        self._flowarea = farea

    @property
    def massflowrate(self):
        """
        Get inlet mass flow rate
        :return: mass flow rate [g/sec] (double scalar)
        """
        if self._flowratemode == 0:
            return self._massflowrate
        else:
            return self.converttomassflowrate()

    @massflowrate.setter
    def massflowrate(self, mflowrate):
        """
        Set inlet mass flow rate
        :param mflowrate: mass flow rate [g/sec] (double scalar)
        """
        if mflowrate <= 0.0:
            print(Color.RED + "** invalid mass flow rate value", end=Color.END)
            return
        # reset the flow rates
        self._volflowrate = 0.0
        self._velocity = 0.0
        self._SCCM = 0.0
        # set flow rate mode to mass flow rate
        self._flowratemode = 0
        self._inletflowrate[self._flowratemode] = mflowrate
        self._massflowrate = mflowrate

    @property
    def volflowrate(self):
        """
        Get inlet volumetric flow rate
        :return: mass flow rate [cm3/sec] (double scalar)
        """
        if self._flowratemode == 1:
            return self._volflowrate
        else:
            return self.converttovolflowrate()

    @volflowrate.setter
    def volflowrate(self, vflowrate):
        """
        Set inlet volumetric flow rate
        :param vflowrate: volumetric flow rate [cm3/sec] (double scalar)
        """
        if vflowrate <= 0.0:
            print(Color.RED + "** invalid volumetric flow rate value", end=Color.END)
            return
        # reset the flow rates
        self._massflowrate = 0.0
        self._velocity = 0.0
        self._SCCM = 0.0
        # set flow rate mode to volumetric flow rate
        self._flowratemode = 1
        self._inletflowrate[self._flowratemode] = vflowrate
        self._volflowrate = vflowrate

    @property
    def sccm(self):
        """
        Get inlet SCCM volumetric flow rate
        :return: SCCM volumetric flow rate [standard cm3/min] (double scalar)
        """
        if self._flowratemode == 3:
            return self._SCCM
        else:
            return self.converttoSCCM()

    @sccm.setter
    def sccm(self, vflowrate):
        """
        Set inlet volumetric flow rate in SCCM
        :param vflowrate: SCCM volumetric flow rate [standard cm3/min] (double scalar)
        :return: None
        """
        if vflowrate <= 0.0:
            print(
                Color.RED + "** invalid SCCM volumetric flow rate value", end=Color.END
            )
            return
        # reset the flow rates
        self._massflowrate = 0.0
        self._volflowrate = 0.0
        self._velocity = 0.0
        # set flow rate mode to volumetric flow rate
        self._flowratemode = 3
        self._inletflowrate[self._flowratemode] = vflowrate
        self._SCCM = vflowrate

    @property
    def velocity(self):
        """
        Get inlet gas velocity
        :return: velocity [cm/sec] (double scalar)
        """
        if self._flowratemode == 2:
            return self._velocity
        else:
            if self._haveflowarea:
                # have flow area
                if self._flowratemode == 1:
                    vrate = self._volflowrate
                else:
                    vrate = self.converttovolflowrate()
                # convert volumetric flow rate to velocity
                return vrate / self._flowarea
            else:
                # flow area not defined
                print(
                    Color.PURPLE + "** flow area value is not given for this inlet",
                    end=Color.END,
                )
                return 0.0

    @velocity.setter
    def velocity(self, vel):
        """
        Set inlet velocity
        :param vel: velocity [cm/sec] (double scalar)
        :return: None
        """
        if vel <= 0.0:
            print(Color.RED + "** invalid inlet velocity value", end=Color.END)
            return
        # reset the flow rates
        self._massflowrate = 0.0
        self._volflowrate = 0.0
        self._SCCM = 0.0
        # set flow rate mode to velocity
        self._flowratemode = 2
        self._inletflowrate[self._flowratemode] = vel
        self._velocity = vel

    @property
    def velocitygradient(self):
        """
        Get inlet gas axial velocity gradient (for premixed, oppdif, and spin)
        or radial velocity spreading rate (v_r/r) at the inlet.
        :return: velocity gradient [1/sec] (double scalar)
        """
        return self._velgrad

    @velocitygradient.setter
    def velocitygradient(self, velgrad):
        """
        Set inlet axial velocity gradient
        :param velgrad: axial velocity gradient [1/sec] (double scalar)
        :return: None
        """
        if velgrad <= 0.0:
            print(
                Color.RED + "** invalid inlet radial velocity spreading rate value",
                end=Color.END,
            )
            return
        # set velocity gradient
        self._velgrad = velgrad

from cy23 import cythermWF as th  # Import thermodynamic workflow module
import numpy as np  # Import NumPy for numerical calculations

# Class representing a single turbine S1 plane (cascade plane)
class S1T:
    FAR = 0.0  # Fuel-to-air ratio (global class variable)

    def __init__(self, radii, omega):
        """
        Initialize the turbine S1 plane with blade station radii and rotational speed.
        """
        # Create blade stations for Stator and Rotor (leading and trailing edges)
        self.SL = BladeStation('SL', radii[0])  # Stator Leading
        self.ST = BladeStation('ST', radii[1])  # Stator Trailing
        self.RL = BladeStation('RL', radii[2])  # Rotor Leading
        self.RT = BladeStation('RT', radii[3])  # Rotor Trailing

        # Assign rotational speed to all blade stations
        self.SL.omega = self.ST.omega = omega
        self.RL.omega = self.RT.omega = omega

        # Initialize required variables
        self.dht = 0.0  # Specific work [J/kg]
        self.loss_R = 0.0  # Rotor loss coefficient
        self.loss_S = 0.0  # Stator loss coefficient
        self.alpha1 = 0.0  # Inlet flow angle [degrees]
        self.alpha2 = 0.0  # Outlet flow angle [degrees]
        self.PtIn = 101325.0  # Total pressure at inlet [Pa]
        self.TtIn = 288.15  # Total temperature at inlet [K]
        self.DoR = 0.0  # Degree of Reaction

    def __call__(self):
        """
        Perform calculations for the turbine S1 plane.
        """
        # Set inlet conditions for the Stator leading edge
        self.SL.setTtPt(self.TtIn, self.PtIn)
        self.SL.FAR = S1T.FAR  # Assign fuel-to-air ratio
        self.SL.alpha = self.alpha1  # Set inlet flow angle

        # Copy Stator leading edge properties to Stator trailing edge
        self.ST.copy(self.SL)
        self.ST.setTtPt(self.SL.Tt, self.SL.Pt - self.loss_S * (self.SL.Pt - self.SL.Ps))  # Account for stator losses
        self.ST.alpha = self.alpha2  # Set outlet flow angle

        # Copy Stator trailing edge properties to rotor leading edge
        self.RL.copy(self.ST)
        self.RL.Vu = self.ST.Vu * self.ST.r / self.RL.r  # Adjust tangential velocity based on radius

        # Copy rotor leading edge properties to rotor trailing edge
        self.RT.copy(self.RL)
        self.RT.Vu = (self.RL.U * self.RL.Vu - self.dht) / self.RT.U  # Adjust tangential velocity for work extraction
        self.RT.setHtPtrel(self.RL.ht - self.dht, self.RL.rel.Pt - self.loss_R * (self.RL.rel.Pt - self.RL.Ps))  # Account for rotor losses

        # Calculate the degree of reaction
        self.DoR = (self.RL.hs - self.RT.hs) / (self.SL.hs - self.RT.hs)

    def data(self):
        """
        Generate a dictionary containing turbine stage data.
        """
        data_dict = {
            "SL": self.SL.data(),  # Data for Stator leading edge
            "ST": self.ST.data(),  # Data for Stator trailing edge
            "RL": self.RL.data(),  # Data for rotor leading edge
            "RT": self.RT.data(),  # Data for rotor trailing edge
            "DoR": self.DoR,  # Degree of Reaction
            "dht [kJ/kg]": self.dht / 1000.0,  # Specific work in kJ/kg
            "DTt [K]": self.RT.Tt - self.SL.Tt,  # Total temperature difference
            "loss_R": self.loss_R,  # Rotor loss coefficient
            "loss_S": self.loss_S,  # Stator loss coefficient
        }

        # Flatten nested dictionaries for easier access
        return flatten_dict(data_dict)


# Class representing a blade station in the turbine
class BladeStation(th.flowStation):
    
    def __init__(self, ID, r=None):
        """
        Initialize a blade station with an ID and radius.
        """
        super().__init__(ID)  # Initialize parent class
        self.rel = th.flowStation(ID)  # Relative flow station
        self._r = r  # Radius
        self._omega = None  # Rotational speed
        self._Vax = None  # Axial velocity
        self._U = None  # Blade speed
        self._alpha = 0.0  # Absolute flow angle
        self._beta = 0.0  # Relative flow angle
        self._V = None  # Absolute velocity
        self._Vrel = None  # Relative velocity
        self._Vu = 0.0  # Tangential velocity
        self._Vurel = 0.0  # Relative tangential velocity

    # Property for blade speed
    @property
    def U(self):
        return self._U

    # Property for axial velocity
    @property
    def Vax(self):
        return self._Vax

    @Vax.setter
    def Vax(self, value):
        self._Vax = value
        if self._omega is not None:
            self._update("alpha")  # Update dependent properties

    # Property for radius
    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value
        if self._Vax is not None:
            self._update("alpha")  # Update dependent properties

    # Property for rotational speed
    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, value):
        self._omega = value
        self._U = self._omega * self._r  # Calculate blade speed
        if self._Vax is not None:
            self._update("alpha")  # Update dependent properties

    # Property for absolute flow angle
    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value
        self._update("alpha")  # Update dependent properties

    # Property for relative flow angle
    @property
    def beta(self):
        return self._beta

    @beta.setter
    def beta(self, value):
        self._beta = value
        self._update("beta")  # Update dependent properties

    # Property for absolute velocity
    @property
    def V(self):
        return self._V

    @V.setter
    def V(self, value):
        self._V = value
        self._update("V")  # Update dependent properties

    # Property for relative velocity
    @property
    def Vrel(self):
        return self._Vrel

    @Vrel.setter
    def Vrel(self, value):
        self._Vrel = value
        self._update("Vrel")  # Update dependent properties

    # Property for tangential velocity
    @property
    def Vu(self):
        return self._Vu

    @Vu.setter
    def Vu(self, value):
        self._Vu = value
        self._update("Vu")  # Update dependent properties

    # Property for relative tangential velocity
    @property
    def Vurel(self):
        return self._Vurel

    @Vurel.setter
    def Vurel(self, value):
        self._Vurel = value
        self._update("Vurel")  # Update dependent properties

    def _update(self, property_name):
        """
        Update dependent properties based on the changed property.
        """
        # Update logic for each property
        if property_name == "alpha":
            self._V = self._Vax / np.cos(np.deg2rad(self._alpha))
            self._Vu = np.sqrt(self._V**2 - self._Vax**2)
            self.Vurel = self._Vu - self._U
        elif property_name == "beta":
            self._Vrel = self._Vax / np.cos(np.deg2rad(self._beta))
            self._Vurel = np.sqrt(self._Vrel**2 - self._Vax**2)
            self.Vu = self._Vurel + self._U
        elif property_name == "V":
            self._alpha = np.arctan(self._Vu / self._Vax) * 180 / np.pi
            self._Vu = np.sqrt(self._V**2 - self._Vax**2)
            self.Vurel = self._Vu - self._U
        elif property_name == "Vu":
            self._V = np.sqrt(self._Vax**2 + self._Vu**2)
            self._alpha = np.arctan(self._Vu / self._Vax) * 180 / np.pi
            self._Vurel = self._Vu - self._U
            self._Vrel = np.sqrt(self._Vax**2 + self._Vurel**2)
            self._beta = np.arctan(self._Vurel / self._Vax) * 180 / np.pi
        elif property_name == "Vurel":
            self._Vrel = np.sqrt(self._Vax**2 + self._Vurel**2)
            self._beta = np.arctan(self._Vurel / self._Vax) * 180 / np.pi
            self._Vu = self._Vurel + self._U
            self._V = np.sqrt(self._Vax**2 + self._Vu**2)
            self._alpha = np.arctan(self._Vu / self._Vax) * 180 / np.pi

        # Update static properties
        self.setStaticsV(self._V)
        self.rhos = self.Ps / (self.R * self.Ts)  # Density
        htrel = self.hs + self._Vrel**2 / 2  # Relative total enthalpy
        self.rel.setHtS(htrel, self.s)  # Update relative station properties
        self.rel.setStaticsV(self._Vrel)

    def setHtPtrel(self, ht, Ptrel):
        """
        Set total enthalpy and relative total pressure.
        """
        htrel = ht - (0.5 * self._V**2) + (0.5 * self._Vrel**2)
        self.rel.setHtPt(htrel, Ptrel)
        self.setHtS(ht, self.rel.s)
        self.setStaticsV(self._V)

    def data(self):
        """
        Generate a dictionary of station data using variable names as keys.
        """
        data_dict = {
            "Pt [kPa]": float(self.Pt) / 1000.0,
            "Ps [kPa]": float(self.Ps) / 1000.0,
            "Tt [K]": float(self.Tt),
            "Ts [K]": float(self.Ts),
            "ht [kJ/kg]": float(self.ht) / 1000.0,
            "hs [kJ/kg]": float(self.hs) / 1000.0,
            "Pt rel [kPa]": float(self.rel.Pt) / 1000.0,
            "Tt rel [K]": float(self.rel.Tt),
            "Vax [m/s]": float(self.Vax),
            "Alpha [deg]": float(self.alpha),
            "Beta [deg]": float(self.beta),
            "V [m/s]": float(self.V),
            "Vrel [m/s]": float(self.Vrel),
            "Vu [m/s]": float(self.Vu),
            "Vurel [m/s]": float(self.Vurel),
            "U [m/s]": float(self.U),
        }

        return data_dict


"""
Auxiliary functions.
"""

def flatten_dict(d, parent_key='', sep='_'):
    """
    Flatten a nested dictionary into a single-level dictionary.
    Keys from nested dictionaries are concatenated with the parent key.
    """
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)
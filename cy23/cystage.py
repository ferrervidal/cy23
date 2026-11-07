from cy23.cythermWF import flowStation, gamma
import numpy as np

class VelStation(flowStation):

    def __init__(self, name=None):
        super().__init__(name)
        self._rel=flowStation() # Aux relative flow station for calcs.
        self.U=None # Blade speed.
        self.Vax=None # Axial velocity.
        self._alpha=None # Angle between axial and stationary frame velocity.
        self._beta=None # Angle between axial and relative frame velocity.
        self._Vr=None #  Relative frame velocity.
        self._Vu=None #  Stationary frame tangential velocity.
        self._Vru=None #  Relative frame tangential velocity.
        self._Ttr=None # Relative frame total temperaure.
        self._Ptr=None # Relative frame total pressure.
        self._htr=None # Relative frame total enthalpy.
    
    def setThermo(self):
        self.setStaticsV(self.V)
        self._htr=self._hs+self._Vr**2/2
        self._rel.copy(self)
        self._rel.setHtS(self.htr, self.s)
        self._Ttr=self._rel.Tt
        self._Ptr=self._rel.Pt
    
    @property
    def Vr(self):
        return self._Vr
    
    @property
    def Vu(self):
        return self._Vu
    
    @property
    def Vru(self):
        return self._Vru
    
    @property
    def Ptr(self):
        return self._Ptr
    
    @property
    def Ttr(self):
        return self._Ttr
    
    @property
    def htr(self):
        return self._htr
    
    @property
    def alpha(self):
        return self._alpha
    
    @property
    def beta(self):
        return self._beta
    
    @alpha.setter
    def alpha(self, value):
        self._alpha=value
        self._V=self.Vax/(np.cos(np.radians(self._alpha)))
        self._Vu=self.Vax*np.tan(np.radians(self._alpha))
        self._Vru=self._Vu-self.U
        self._Vr=np.sqrt(self.Vax**2+self._Vru**2)
        self._beta=np.degrees(np.arctan(self._Vru/self.Vax))
        self.setThermo()

    @beta.setter
    def beta(self, value):
        self._beta=value
        self._Vr=self.Vax/(np.cos(np.radians(self._beta)))
        self._Vru=self.Vax*np.tan(np.radians(self._beta))
        self._Vu=self._Vru+self.U
        self.alpha=np.degrees(np.arctan(self._Vu/self.Vax))

    @Vu.setter
    def Vu(self, value):
        self._Vu=value
        self.alpha=np.degrees(np.arctan(self._Vu/self.Vax))
    
    @Vru.setter
    def Vru(self, value):
        self._Vru=value
        self.Vu=self._Vru+self.U

    def setPtr(self, Ptr):
        self._Ptr=Ptr
        self._htr=self._ht-self._V**2/2+self._Vr**2/2
        self._rel.copy(self)
        self._rel.setHtPt(self._htr, self._Ptr)
        self._rel.setStaticsV(self.Vr)
        Ps=self._rel._Ps
        Ts=self._rel._Ts
        adc=gamma(Ts, self.FAR)
        Pt=Ps*(self.Tt/Ts)**(adc/(adc-1))
        self.setTtPt(self.Tt, Pt)


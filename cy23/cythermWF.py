"""
WFgasModel.py

An implementation of the Walsh & Fletcher thermodynamic 
gas model for air and kerosene combustion products.

Ref: 
    
    Walsh, P.P, Fletcher, P. (2004), "Gas Turbine Performance",
    Second Edition, Wiley (Blackwell Science). 
    Chapter 3, Pg 116-117

Author: Ferrer-Vidal, 2023
"""

import numpy as np
import scipy as sp

"""
Polynomial coefficients.
"""
# Air
A=np.array([
    0.992313, 
    0.236688,
    -1.852148,
    6.083152,
    -8.893933,
    7.097112,
    -3.234725,
    0.794571,
    -0.081873,
    0.422178,
    0.001053
])
# Combustion Products
B=np.array([
    -0.718874,
    8.747481,
    -15.863157,
    17.254096,
    -10.233795,
    3.081778,
    -0.361112,
    -0.003919,
    0.055593,
    -0.0016079
])

"""
Thermodynamic property functions.
"""
def cp(T, FAR): 
    """
    Specific heat capacity in [J/(kg-K)].
    """
    Tz=T/1000
    air_cp=A[0]+A[1]*Tz+A[2]*Tz**2+A[3]*Tz**3+A[4]*Tz**4+A[5]*Tz**5+A[6]*Tz**6+A[7]*Tz**7+A[8]*Tz**8 # cp of pure air.
    comb_cp= FAR/(1+FAR)*(B[0]+B[1]*Tz+B[2]*Tz**2+B[3]*Tz**3+B[4]*Tz**4+B[5]*Tz**5+B[6]*Tz**6+B[7]*Tz**7) # contribution of combustion products.
    return (air_cp+comb_cp)*1E3

def h(T, FAR):
    """
    Specific enthalpy in [J/kg].
    """
    Tz=T/1000
    air_h=A[0]*Tz+(A[1]*Tz**2)/2+(A[2]*Tz**3)/3+(A[3]*Tz**4)/4+(A[4]*Tz**5)/5+(A[5]*Tz**6)/6+(A[6]*Tz**7)/7+(A[7]*Tz**8)/8+(A[8]*Tz**9)/9+A[9] # enthalpy of pure air.
    comb_h=FAR/(1+FAR)*(B[0]*Tz+(B[1]*Tz**2)/2+(B[2]*Tz**3)/3+(B[3]*Tz**4)/4+(B[4]*Tz**5)/5+(B[5]*Tz**6)/6+(B[6]*Tz**7)/7+(B[7]*Tz**8)/8+B[8]) # contribution of combustion products.
    return (air_h+comb_h)*1E6


def intCpQT(T1, T2, FAR): 
    """
    F2-F1 in [J/(kg-K)]. This is int(cp/T)*dT between T1 and T2.
    """
    Tz1=max(T1/1000,1e-4)
    Tz2=max(T2/1000,1e-4)

    air_F2=A[0]*np.log(Tz2)+A[1]*Tz2+A[2]/2*Tz2**2+A[3]/3*Tz2**3+A[4]/4*Tz2**4+A[5]/5*Tz2**5+A[6]/6*Tz2**6+A[7]/7*Tz2**7+A[8]/8*Tz2**8+A[10] #  pure air.
    comb_F2= FAR/(1+FAR)*(B[0]*np.log(Tz2)+B[1]*Tz2+B[2]/2*Tz2**2+B[3]/3*Tz2**3+B[4]/4*Tz2**4+B[5]/5*Tz2**5+B[6]/6*Tz2*6+B[7]/7*Tz2**7+B[9]) # contribution of combustion products.
    
    air_F1=A[0]*np.log(Tz1)+A[1]*Tz1+A[2]/2*Tz1**2+A[3]/3*Tz1**3+A[4]/4*Tz1**4+A[5]/5*Tz1**5+A[6]/6*Tz1**6+A[7]/7*Tz1**7+A[8]/8*Tz1**8+A[10] #  pure air.
    comb_F1= FAR/(1+FAR)*(B[0]*np.log(Tz1)+B[1]*Tz1+B[2]/2*Tz1**2+B[3]/3*Tz1**3+B[4]/4*Tz1**4+B[5]/5*Tz1**5+B[6]/6*Tz1*6+B[7]/7*Tz1**7+B[9]) # contribution of combustion products.
    
    return (air_F2+comb_F2-(air_F1+comb_F1))*1E3

def R(FAR=0.0):
    """
    Specific gas constant in [J/(kg-K)].
    """
    return (287.05-0.00990*FAR+1E-07*FAR**2)

def dS(T1,P1,T2,P2,FAR=0.0):
    """
    Specific entropy in [J/(kg-K)].
    """
    P2=max(P2,1)
    P1=max(P1,1)
    return (intCpQT(T1,T2,FAR)-R(FAR)*np.log(P2/P1))

def gamma(T, FAR=0):
    """
    Specific heat capacity ratio.
    """
    cpp=cp(T, FAR)
    RR=R(FAR)
    return cpp/(cpp-RR)



class flowStation:
    """
    A class implementing above WF functions as 
    flow station objects that can be used to build
    simple cycle models. 
    """

    # Reference T&P for corrected mass flow. 
    Tiso=288.15
    Piso=101325
    # Reference T&P for ht and s. 
    Tref=1.0
    Pref=101325

    def __init__(self, name=None):
        self.ID=name
        self._W=0.0
        self._Tt=288.15
        self._Pt=101325
        self._FAR=0.0
        self.__sett__()
        # Statics
        self._Ts=None
        self._hs=None
        self._Ps=None
        self._M=None
        self._a=None
        self._V=None
        self._A=None

    def __sett__(self):
        self._ht=h(self._Tt, self._FAR)-h(flowStation.Tref, self._FAR)
        self._cp=cp(self._Tt, self._FAR)
        self._s=dS(flowStation.Tref,flowStation.Pref, self._Tt, self._Pt, self._FAR)
        self._R=R(self._FAR)
        self._gamma=gamma(self._Tt, self._FAR)
        self._Wc=self._W*np.sqrt(self._Tt/flowStation.Tiso)/(self._Pt/flowStation.Piso)
        self._mu=(1.458E-6)*self._Tt**(3/2)/(self._Tt+110.4) # Simple Sutherland's equation for air.
        self._RNI=self._Pt/flowStation.Pref*np.sqrt((287*flowStation.Tref)/(self._R*self._Tt))*((1.716E-5)/self._mu) # Reynolds number index with air as reference.

    @property
    def W(self):
        return self._W
    
    @property
    def Wc(self):
        return self._Wc

    @property
    def Tt(self):
        return self._Tt
    
    @property
    def Pt(self):
        return self._Pt
    
    @property
    def FAR(self):
        return self._FAR
    
    @property
    def ht(self):
        return self._ht
    
    @property
    def hs(self):
        return self._hs
    
    @property
    def cp(self):
        return self._cp
    
    @property
    def s(self):
        return self._s
    
    @property
    def R(self):
        return self._R
    
    @property
    def gamma(self):
        return self._gamma
    
    @property
    def mu(self):
        return self._mu

    @property
    def RNI(self):
        return self._RNI

    @property
    def Ts(self):
        return self._Ts
    
    @property
    def Ps(self):
        return self._Ps
    
    @property
    def M(self):
        return self._M
    
    @property
    def V(self):
        return self._V
    
    @property
    def A(self):
        return self._A
    
    @property
    def a(self):
        return self._a

    @W.setter
    def W(self, W):
        self._W=W
        self._Wc=self._W*np.sqrt(self._Tt/flowStation.Tiso)/(self._Pt/flowStation.Piso)


    @FAR.setter
    def FAR(self, far):
        self._FAR=far
        self.__sett__()

    def setTtPt(self, Tt=None, Pt=None):
        if Tt:
            self._Tt=Tt
        if Pt:
            self._Pt=Pt
        self.__sett__()
        
    def setHtPt(self, ht, Pt=None):
        def obj(T):
            return ((h(T, self.FAR)-h(flowStation.Tref, self.FAR))-ht)/max(abs(ht),1)*100
        # Tt=sp.optimize.fsolve(obj, self._Tt)[0]
        Tt=sp.optimize.newton(obj, self._Tt)
        self.setTtPt(Tt, Pt)

    def setHtS(self, ht, s=None):
        def obj(T):
            return ((h(T, self.FAR)-h(flowStation.Tref, self.FAR))-ht)/max(abs(ht),1)*100
        Tt=sp.optimize.fsolve(obj, self._Tt)[0]
        self.setTtS(Tt, s)

    def setTtS(self,Tt, s):
        def obj(P):
            return (dS(flowStation.Tref,flowStation.Pref, Tt, P, self._FAR)-s)/max(abs(s),1)*100.0
        Pt=sp.optimize.fsolve(obj, self._Pt)[0]
        self.setTtPt(Tt, Pt)

    def setPtS(self,Pt, s):
        def obj(T):
            return (dS(flowStation.Tref,flowStation.Pref, T, Pt, self._FAR)-s)/max(abs(s),1)*100.0
        # Tt=sp.optimize.fsolve(obj, self._Tt*1.1)[0]
        Tt=sp.optimize.newton(obj, self._Tt)
        self.setTtPt(Tt, Pt)

    def copy(self, toCopy):
        self._W=toCopy.W
        self._Tt=toCopy.Tt
        self._Pt=toCopy.Pt
        self._FAR=toCopy.FAR
        self.__sett__()

    def mix(self, toMix):
        W=self.W+toMix.W
        PtIn=self.Pt
        ht=(self.W*self.ht+toMix.W*toMix.ht)/W
        FAR=(self.W*self.FAR+toMix.W*toMix.FAR)/W
        # s=(self.s*self.W+toMix.s*toMix.W)/W
        self.W=W
        self.FAR=FAR
        self.setHtPt(ht, self.Pt)
        # self.setTtS(self.Tt, s)
        self.setTtPt(self.Tt, PtIn)

    def procPolyEff(self, init, PR, eff):
        PR=max(0, PR)
        self.copy(init)
        Pt=init.Pt*PR
        if PR>1:
            PtS=init.Pt*(PR)**(1/eff)
        else:
            PtS=init.Pt*(PR)**(eff)
        self.setPtS(PtS, init.s)
        self.setHtPt(self.ht,Pt)

    def procIsenEff(self, init, PR, eff):
        self.copy(init)
        Pt=init.Pt*PR
        self.setPtS(Pt, init.s)

        if PR>1.0:
            dh_isen=self.ht-init.ht
            dh=dh_isen/eff
            htOut=init.ht+dh
        else:
            dh_isen=init.ht-self.ht
            dh=dh_isen*eff
            htOut=init.ht-dh

        self.setHtPt(htOut,Pt)
    
    def setStaticsM(self, M):
        Tt=self.Tt
        Pt=self.Pt
        adc=self.gamma
        Ts=Tt
        TsOld=1.0
        while abs(Ts-TsOld)/TsOld>0.001:
            TsOld=Ts
            Ts=Tt/(1+((adc-1)/2)*M**2)
            adc=gamma(Ts, self.FAR)
        Ps= Pt*(Ts/Tt)**(adc/(adc-1))
        a=np.sqrt(adc*self.R*Ts)
        V=a*M
        rho=Ps/(self.R*Ts)
        self._Ts=Ts
        self._Ps=Ps
        self._M=M
        self._a=a
        self._V=V
        if abs(V)>0:
            self._A=abs(self.W)/(rho*V)
        else:
            self._A=float('inf')
        self._hs=self.ht-V**2/2

    def setStaticsPs(self, Ps):
        adc=self.gamma
        Pt=self.Pt
        Tt=self.Tt
        Ts=Tt
        TsOld=1.0
        while abs(Ts-TsOld)/TsOld>0.001:
            TsOld=Ts
            Ts=Tt*(Ps/Pt)**((adc-1)/adc)
            adc=gamma(Ts, self.FAR)
        a=np.sqrt(adc*self.R*Ts) 
        rho=Ps/(self.R*Ts)
        hdyn=abs(h(Tt, self.FAR)-h(Ts, self.FAR))
        V=np.sqrt(2*hdyn)
        M=V/a
        self._Ts=Ts
        self._Ps=Ps
        self._M=M
        self._a=a
        self._V=V
        if abs(V)>0:
            self._A=abs(self.W)/(rho*V)
        else:
            self._A=float('inf')
        self._hs=self.ht-V**2/2
        
    def setStaticsA(self, A, supers=False):
        def obj(Mn):
            self.setStaticsM(Mn)
            return (self.A-A)/A
        if supers:
            sol=sp.optimize.root_scalar(obj, bracket=[1.0,1000], method='brentq')
        else:
            sol=sp.optimize.root_scalar(obj, bracket=[0,1.0], method='brentq')
        Mn=sol.root
        # # Set initial guess for Mach.
        # if self.M is None:
        #     if supers:
        #         MnGuess=3.0
        #     else:
        #         MnGuess=0.1
        # else:
        #     MnGuess=self.M
        # Mn=sp.optimize.fsolve(obj, MnGuess)[0]
        # if Mn>1.0 and not supers:
        #     Mn=sp.optimize.fsolve(obj, 0.1)[0]
        # elif Mn<1.0 and supers:
        #     Mn=sp.optimize.fsolve(obj, 3.0)[0] 
 
    def setStaticsV(self, V):
        Tt=self.Tt
        Pt=self.Pt
        cpp=self.cp
        Ts=Tt
        TsOld=1.0
        while abs(Ts-TsOld)/TsOld>0.001:
            TsOld=Ts
            cpp=cp(Ts, self.FAR)
            Ts=Tt-V**2/(2*cpp)
        adc=gamma(Ts, self.FAR)
        Ps= Pt*(Ts/Tt)**(adc/(adc-1))
        a=np.sqrt(adc*self.R*Ts)
        M=V/a
        rho=Ps/(self.R*Ts)
        self._Ts=Ts
        self._Ps=Ps
        self._M=M
        self._a=a
        self._V=V
        if abs(V)>0:
            self._A=abs(self.W)/(rho*V)
        else:
            self._A=float('inf')
        self._hs=self.ht-V**2/2

    def setPsA(self, Ps, Atarg):

        PtMax=self.Pt*1.1
        def obj(PtGuess): 
            self.setTtPt(self.Tt, PtGuess) 
            self.setStaticsPs(Ps) 
            return (self.A-Atarg)/Atarg 
        sol=sp.optimize.root_scalar(obj, bracket=[Ps,PtMax], method='brentq')
        Pt=sol.root
        self.setTtPt(self.Tt, Pt)
        self.setStaticsPs(Ps)
        


        



# F01=WF_flowStation('01')

# F01.setTtPt(300.0)
# print(F01.Tt)
# print(F01.s)
# print(F01.cp)
# print(F01.gamma)

# F01.setPtS(300.0, F01.s)

# print('---')
# print(F01.Tt)
# print(F01.s)
# print(F01.cp)
# print(F01.gamma)


# F01.setPtS(100.0, F01.s)
# print('---')
# print(F01.Tt)
# print(F01.s)
# print(F01.cp)
# print(F01.gamma)
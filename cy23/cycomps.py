"""
CY23.CYCOMPS

This module defines a set of classes for simulating the components of a thermodynamic cycle.

Each component is represented by a class that can simulate its operation based on input parameters such as mass flow rate, 
temperature, pressure, and specific thermodynamic properties. Components can be connected in series (upstream/downstream relationships), 
and the flow properties are updated dynamically as the system operates. 

Classes:
- component: A base class for components that stores common attributes and provides methods for updating variables.
- ambient: Represents ambient conditions, including temperature, pressure, and Mach number.
- burner: Models a combustion chamber, including fuel-to-air ratio, efficiency, and power generation.
- compressor: Models a compressor stage, including pressure ratio, efficiency, and off-design scaling.
- duct: Represents a duct in the flow path, accounting for pressure drops.
- turbine: Models a turbine, including expansion ratio, efficiency, and off-design performance.
- bleedSource: Represents a source of bleed air, drawing mass from the flow.
- bleedSink: A component that receives bleed air from a source and mixes it into the flow.
- mixer: Combines two flow streams into a single output, adjusting pressures and areas as needed.
- nozzle: Models a convergent or convergent-divergent nozzle, including thrust calculation and choked flow conditions.
- splitter: Splits a flow into two streams, adjusting the mass flow split ratio.
- diffuser: A component that slows down a flow and increases its pressure, with the option for off-design operation.
- flowEnd: Represents the end of the flow path, where additional post-processing can be applied.
- shaft: A mechanical element connecting various components, used to track and tally torques within the system.

Each component is initialized with default parameters but can be customized with specific properties such as 
efficiency, mass flow, and design constraints. The system can simulate both design and off-design conditions, 
allowing for performance analysis across a range of operating points.
"""


from cy23.cythermWF import flowStation
import numpy as np

class component:

    def __init__(self):
        """
        Callables are parameters that must be reocomputed if the model state changes. 
        A good example are bleed schedules, which may depend on some model parameter that
        has changed. When each component is executed, all parameters in the callable list
        are executed. Many callables may be cyutil.Table objects, which must be looked up
        anytime the model is executed. 
        """
        self.callables=[]

    def __updatevars__(self):
        if self.callables: 
            for name in self.callables:
                var = getattr(self, name)
                if callable(var):
                    val=var()
                else:
                    val=var
                setattr(self, name, val)

class ambient(component):

    def __init__(self, name=None):
        self.station=None
        self.outlet=flowStation()
        self.T=288.15
        self.P=101325
        self.M=0.0
        self.W=0.0
        self.inlet=self.outlet
        self.ID=name
        self._FS=flowStation()
        self.callables=[]

    def __call__(self):
        self.__updatevars__()
        if self.M==0:
            self.outlet.W=self.W
            self.outlet.setTtPt(self.T, self.P)
            self.V=0.0
        else:
            self._FS.setTtPt(self.T, self.P)
            Ts=self.T
            Tt=Ts*1.2
            Ps=self.P
            TtOld=1.0
            while abs((Tt-TtOld))/TtOld>0.001:
                TtOld=Tt
                gamma=self._FS.gamma
                Tt=Ts*(1+(gamma-1)/2*self.M**2)
                Pt= self.P*(Tt/Ts)**(gamma/(gamma-1))
                self._FS.setTtPt((Tt+Ts)/2, (Pt+Ps)/2)
            self.outlet.setTtPt(Tt, Pt)
            self.outlet.W=self.W
            self.outlet.setStaticsM(self.M)
            self.V=self.outlet.V
        self.inlet=self.outlet
            
class burner(component):
    
    def __init__(self, name=None, upstream=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.outlet=flowStation()
        self.eff=1.0
        self._FS=flowStation()
        self.upstream=None
        self.FAR=0.0
        self.dPqP=0.0
        self.callables=['dPqP', 'eff']
        self.Wfuel_in=None
        self.OFF=False
        
        # Lower heating value of kerosene.
        self.LHV=43.0e6 # [J/kg]

        self.COT=None
        self.ERR=1.0
    def __call__(self):
        self.__updatevars__()
        PtOut=self.inlet.Pt*(1-self.dPqP)
        if self.OFF:
            self.FAR=0.0
            self.ERR=0.0
            self.Wf=0.0
            self.outlet.copy(self.inlet)
            self.outlet.setTtPt(self.inlet.Tt, PtOut)
        else:
            wair=self.inlet.W*(1-self.inlet.FAR)
            if self.Wfuel_in:
                self.FAR=self.Wfuel_in/wair
                
            wfuel=self.FAR*wair
            dH=wfuel*self.LHV*self.eff

            self._FS.copy(self.inlet)
            self._FS.FAR=(self.inlet.W*self.inlet.FAR+wfuel)/wair # if air already had FAR, we need to add to it.
            self._FS.W=self.inlet.W+wfuel
            htOut=(self.inlet.ht*self.inlet.W+dH)/self._FS.W

            self._FS.setHtPt(htOut, PtOut)
            self.outlet.copy(self._FS)
            self.Wf=wfuel
            if self.COT:
                self.ERR=(self.outlet.Tt-self.COT)/self.COT

    def __postSolve__(self):
        self.Wfuel_in=None
        
     


class compressor(component):

    def __init__(self, name=None, upstream=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.outlet=flowStation()
        self.PR=1.0
        self.eff=1.0
        self.poly=True
        self.Nmech=None
        self.power=None
        self.callables=[]

        # OFFDESIGN VARIABLES
        self.OFFDESIGN=False
        self.SCALE=False
        self.beta=10.0
        self.Nmap=1.0
        self.Nc=None
        self.WcMatch=None
        self.MAP=None
        self.ERR=None
        self.shaft=None
        
        # REYNOLDS CORRECTION PARAMETERS
        self.effcorrFlag=True
        self.micronsRough=1.5   # Blade roughness in microns. 
        self.HT=0.5             # Hub to tip ratio of first stage.
        self.HC=3               # Height to chord ratio of first stage.
        self.MnIn=0.5           # Assumed inlet Mach number for rough geometry.
        self.rr=0               # Relative roughness.
        self.fref=None           # Reference Darcy friction factor. 
        self.fsAux=flowStation() 
        
        # MAP SCALING VARIABLES
        self.DPscaleNc=1.0
        self.DPscalePR=1.0
        self.DPscaleWc=1.0
        self.DPaddEff=0.0

        self.AuditScalePR=1.0
        self.AuditScaleWc=1.0
        self.AuditAddEff=0.0
        self.AuditEffScalar=1.0

    def __call__(self):
        if self.AuditEffScalar is not None:
            self.AuditAddEff=(self.AuditEffScalar-1)*0.1
        self.Nmech=self.shaft.Nmech
        self.Nc=self.Nmech/(np.sqrt(self.inlet.Tt/288.15))
        self.__updatevars__()
        if self.OFFDESIGN:
            self.__offdes__()
        elif self.MAP is not None:
            self.__scale__()

        self.outlet.copy(self.inlet)

        if self.poly:
            self.outlet.procPolyEff(self.inlet, self.PR, self.eff)
        else:
            self.outlet.procIsenEff(self.inlet, self.PR, self.eff)
        
        self.dh=self.outlet.ht-self.inlet.ht
        self.power=-self.dh*self.inlet.W

    def __offdes__(self):
        self.ERR=1.0
        self.Nmap=self.Nc/self.DPscaleNc
        mapWc, mapPR, mapEff = self.MAP(self.Nmap, self.beta)
        self.PR=(mapPR-1)*self.DPscalePR*self.AuditScalePR+1
        self.eff=self.effcorr(mapEff+self.DPaddEff+self.AuditAddEff)
        self.WcMatch=mapWc*self.DPscaleWc*self.AuditScaleWc
        self.ERR=(-self.inlet.Wc+self.WcMatch)/self.WcMatch

    def __scale__(self):
        self.DPscaleNc=self.Nc/self.Nmap
        mapWc, mapPR, mapEff = self.MAP(self.Nmap, self.beta)
        self.DPscalePR=(self.PR-1)/(mapPR-1)
        self.DPscaleWc=self.inlet.Wc/(mapWc)
        self.DPaddEff=self.eff-mapEff 

        # Calculate a reference Re and relative roughness.
        self.fsAux.copy(self.inlet)
        self.fsAux.setStaticsM(self.MnIn)
        rhos=self.fsAux.Ps/(self.fsAux.R*self.fsAux.Ts)
        rt=np.sqrt(self.fsAux.A/(np.pi*(1-self.HT**2)))
        rh=rt*self.HT
        Lref=(rt-rh)/self.HC
        ReRef=rhos*self.fsAux.V*Lref/self.fsAux.mu
        self.RNIscale=ReRef/self.inlet.RNI
        self.rr=self.micronsRough*1e-6/Lref
        self.fref= (1/(-1.8*np.log(6.9/ReRef+(self.rr/3.7)**1.11)))**2
    
    def effcorr(self, eff):
        # Reynolds correction based on Haaland equation.
        # Darcy friction factor calculated from Haaland eq.
        # As Re decreases, friction factor increases.
        # Reynolds effect are assumed to account for up to 70% of efficiency
        # defecit, hence 0.7 factor.
        # Reference PTG TN-2409A0.
        if self.effcorrFlag:
            Re=self.inlet.RNI*self.RNIscale
            f = (1/(-1.8*np.log(6.9/Re+(self.rr/3.7)**1.11)))**2
            ReFc=1+0.7*(f/self.fref-1)
        else:
            ReFc=1.0
        return (1-(1-eff)*ReFc)
        

class duct(component):

    def __init__(self, name=None, upstream=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.outlet=flowStation()
        self.dPqP=0.0
        self.dPqPscale=True
        self.callables=['dPqP']

        #OFFDESIGN VARIABLES
        self.OFFDESIGN=False
        self.WcDes=None

        self.Tset=None
        self.Pset=None
        self.ERR=1.0

    def __call__(self):
        self.__updatevars__()
        self.outlet.copy(self.inlet)
        if not self.OFFDESIGN:
            self.WcDes=self.inlet.Wc
        if self.dPqPscale:
            dPqP=self.dPqP/self.WcDes*self.inlet.Wc
        else:
            dPqP=self.dPqP

        self.outlet.setTtPt(self.inlet.Tt, self.inlet.Pt*(1-dPqP))

        if self.Tset:
            self.ERR=(self.inlet.Tt-self.Tset)/self.Tset
        elif self.Pset:
            self.ERR=(self.inlet.Pt-self.Pset)/self.Pset
        else:
            self.ERR=1.0

class turbine(component):

    def __init__(self, name=None, upstream=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.outlet=flowStation()
        self.ER=2.0
        self.eff=1.0
        self.poly=True
        self._FS=flowStation()
        self.capFS=None # used for cooling floe capacity calcs.
        self.power=None
        self.callables=[]
        self.Nmech=None
        self.shaft=None

        # OFFDESIGN VARIABLES
        self.Nmap=1.0
        self.ERmap=3.0
        self.Nc=None
        self.WcMatch=None
        self.MAP = None
        self.OFFDESIGN=False
        self.ERR=None
        self.swirl=None

        # MAP SCALING VARIABLES
        self.DPscaleNc=1.0
        self.DPscaleER=1.0
        self.DPscaleWc=1.0
        self.DPaddEff=0.0

        self.AuditScalePR=1.0
        self.AuditScaleWc=1.0
        self.AuditAddEff=0.0
        self.AuditEffScalar=1.0

    def __call__(self):
        if self.AuditEffScalar is not None:
            self.AuditAddEff=(self.AuditEffScalar-1)*0.1
        self.Nmech=self.shaft.Nmech
        self.Nc=self.Nmech/(np.sqrt(self.inlet.Tt/288.15))
        self.__updatevars__()
        if self.OFFDESIGN:
            self.__offdes__()
        elif self.MAP is not None:
            self.__scale__()

        self.outlet.copy(self.inlet)
        self._FS.copy(self.inlet)
       

        if self.poly:
            self.outlet.procPolyEff(self.inlet, 1/self.ER, self.eff)
        else:
            self.outlet.procIsenEff(self.inlet, 1/self.ER, self.eff)
        
        self.dh=self.outlet.ht-self.inlet.ht
        self.power=-self.dh*self.inlet.W

    def __offdes__(self):
        self.ERR=1.0
        self.Nmap=self.Nc/self.DPscaleNc
        self.ERmap=(self.ER-1)/self.DPscaleER+1
        mapWc, mapEff, mapSwirl = self.MAP(self.Nmap, self.ERmap)
        self.swirl=mapSwirl
        self.eff=mapEff+self.DPaddEff+self.AuditAddEff
        self.WcMatch=mapWc*self.DPscaleWc*self.AuditScaleWc
        self.ERR=(-self.inlet.Wc+self.WcMatch)/self.WcMatch

    def __scale__(self):
        self.DPscaleNc=self.Nc/self.Nmap
        self.DPscaleER=(self.ER-1)/(self.ERmap-1)
        mapWc, mapEff, mapSwirl = self.MAP(self.Nmap, self.ERmap)
        self.swirl=mapSwirl
        self.DPscaleWc=self.inlet.Wc/mapWc
        self.DPaddEff=self.eff-mapEff


class bleedSource(component):
    
    def __init__(self, name=None, upstream=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.outlet=flowStation()
        self.fracW=0.0
        self.ref=None
        self.offtake=flowStation()
        self.callables=['fracW']
    
    def __call__(self):
        self.__updatevars__()
        self.outlet.copy(self.inlet)
        self.offtake.copy(self.inlet)
        if not self.ref:
            Wref=self.inlet.W
        else:
            if hasattr(self.ref, 'inlet'):
                Wref=self.ref.inlet.W
            else:
                Wref=self.ref
        self.outlet.W=self.inlet.W-Wref*self.fracW
        self.offtake.W=Wref*self.fracW

class bleedSink(component):

    def __init__(self, name=None, upstream=None, bleedIn=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.outlet=flowStation()
        self.bleedIn=bleedIn
        self.callables=[]
    
    def __call__(self):
        self.__updatevars__()

        bleed_list = self.bleedIn
        if not isinstance(bleed_list, (list, tuple)):
            bleed_list = [bleed_list]

        self.outlet.copy(self.inlet)

        for bldsrc in bleed_list:
            if hasattr(bldsrc, 'offtake'):
                bld = bldsrc.offtake
            else:
                bld = bldsrc
            self.outlet.mix(bld)

class mixer(component):

    def __init__(self, name=None, upstream1=None, upstream2=None):
        self.station=None
        self.ID=name
        self.inlet=upstream1
        self.outlet=flowStation()
        self.inlet2=upstream2
        self.callables=[]

        self.A1=None
        self.A2=None
        self.M1=None
        self.M2=None
        self.Ps1=None
        self.Ps2=None

        self.fsAux=flowStation()

        # OFFDESIGN VARIABLES
        self.OFFDESIGN=False
        self.ERR=None
    
    def __call__(self):
        self.__updatevars__()
        self.outlet.copy(self.inlet)
        self.outlet.mix(self.inlet2)

        
        if not self.OFFDESIGN: # Design
            if None not in (self.M1, self.M2):
                self.fsAux.copy(self.inlet)
                self.fsAux.setStaticsM(self.M1)
                self.A1=self.fsAux.A
                self.Ps1=self.fsAux.Ps

                self.fsAux.copy(self.inlet2)
                self.fsAux.setStaticsM(self.M2)
                self.A2=self.fsAux.A
                self.Ps2=self.fsAux.Ps
            else: # If statics not specified, use total pressures. 
                self.A1=self.A2=None
                self.Ps1=self.inlet.Pt
                self.Ps2=self.inlet2.Pt
            self.ERR=(self.Ps2-self.Ps1)/self.Ps1

        else: # Off-design
            if None not in (self.A1, self.A2):
                self.fsAux.copy(self.inlet)
                self.fsAux.setStaticsA(self.A1)
                self.Ps1=self.fsAux.Ps

                self.fsAux.copy(self.inlet2)
                self.fsAux.setStaticsA(self.A2)
                self.Ps2=self.fsAux.Ps
                self.ERR=(self.Ps2-self.Ps1)/self.Ps1
            else:
                self.ERR=(self.inlet2.Pt-self.inlet.Pt)/self.inlet.Pt

class nozzle(component):

    def __init__(self, name=None, upstream=None, ambient=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.throat=flowStation()
        self.outlet=flowStation()
        self.Ath=None
        self.Aex=None
        self.Cd=1.0
        self.Cv=1.0
        self.ambient=ambient
        self.Pback=None
        self.CONDI=False
        self.Fg=0.0
        self.callables=['Cd','Cv']
        self.power=None
    
        # OFFDESIGN VARIABLES
        self.OFFDESIGN=False
        self.ERR=None
        self.VARI=False


    def __call__(self):
        self.__updatevars__()
        if self.ambient is not None and hasattr(self.ambient, 'P'):
            Pback = self.ambient.P
        elif self.Pback is not None:
            Pback = self.Pback
        else:
            raise ValueError("No ambient or back pressure specified for nozzle.")
        self.throat.copy(self.inlet)
        self.throat.setStaticsM(1.0)
        if not self.CONDI: 
            # convergent nozzle.
            if self.throat.Ps>Pback: 
                # choked
                AthCalc=self.throat.A/self.Cd
                self.outlet.copy(self.throat)
                self.outlet.setStaticsM(1.0)
            else:
                # not choked.
                self.throat.setStaticsPs(Pback) 
                AthCalc=self.throat.A/self.Cd
                self.outlet.copy(self.throat)
                self.outlet.setStaticsPs(Pback)
            
            if self.OFFDESIGN and not self.VARI:
                self.ERR=(AthCalc-self.Ath)/(self.Ath)
            else:
                self.Ath=AthCalc
                self.Aex=self.Ath
                if self.VARI:
                    self.ERR=0.0
        else:  # NB: Convergent-divergent nozzle has not been tested.
            # Con-di nozzle
            if not self.OFFDESIGN:
                if self.throat.Ps>Pback: 
                    # choked
                    self.Ath=self.throat.A/self.Cd
                    self.outlet.copy(self.throat)
                    self.outlet.setStaticsPs(Pback)
                    self.Aex=self.outlet.A
                else:
                    raise BadDesign('A convergent-divergent nozzle has been designed to run unchoked.')
            else: 
                if self.throat.Ps>Pback: 
                    # choked
                    AthCalc=self.throat.A/self.Cd
                    self.outlet.copy(self.throat)
                    self.outlet.setStaticsA(self.Aex, super=True)
                    self.ERR=(self.outlet.Ps-Pback)/Pback
                else:
                    # unchoked
                    self.throat.setStaticsA(self.Ath)
                    self.outlet.setStaticsA(self.Aex)
                    self.ERR=(self.outlet.Ps-Pback)/Pback

        # Calculate gross thrust
        self.Fg=self.outlet.W*self.outlet.V*self.Cv+self.outlet.A*(self.outlet.Ps-Pback)
        
        # Calculate stream power (kinetic energy flux) at nozzle exit. 
        self.power=self.outlet.W*self.outlet.V**2/2
        

class splitter(component):

    def __init__(self, name=None, upstream=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.outlet1=flowStation()
        self.outlet2=flowStation()
        self.BPR=1.0
        self.callables=['BPR']

    def __call__(self):
        self.__updatevars__()
        self.outlet1.copy(self.inlet)
        self.outlet2.copy(self.inlet)

        self.outlet1.W=self.inlet.W/(self.BPR+1)
        self.outlet2.W=self.BPR*self.inlet.W/(self.BPR+1)

class diffuser(component):

    def __init__(self, name=None, upstream=None, ambient=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.cpr=1.0
        self.ambient=ambient
        self.Pback=101325.0
        self.outlet=flowStation()
        self.fsAux=flowStation()
        self.Mdes=None
        self.AR=1
        self.callables=['cpr']

        # OFFDESIGN VARIABLES
        self.OFFDESIGN=False
        self.ERR=None

    def __call__(self):
        self.ERR=1.0
        self.__updatevars__()
        if self.ambient is not None and hasattr(self.ambient, 'P'):
            Pback = self.ambient.P
        else:
            Pback = self.Pback
        self.fsAux.copy(self.inlet)
        self.outlet.copy(self.inlet)
        if self.OFFDESIGN:
            if self.Ain is None or self.AR is None:
                raise ValueError("Diffuser Ain and/or AR must be set before OFFDESIGN evaluation.")
            else:
                self.fsAux.setStaticsA(self.Ain)
        else:
            ARideal=np.sqrt(1/(1-self.cpr))
            self.AR=max(ARideal, self.AR)

            self.fsAux.setStaticsM(self.Mdes)
            self.Ain=self.fsAux.A
        
        PsOut=self.fsAux.Ps+self.cpr*(self.fsAux.Pt-self.fsAux.Ps)
        
        # Diffuser match error calculated the same way for off and on-design. 
        self.ERR=(PsOut-Pback)/Pback
        
        self.outlet.setPsA(PsOut, self.Ain*self.AR)
       
                    
        
class flowEnd:

    def __init__(self, name=None, upstream=None):
        self.station=None
        self.ID=name
        self.inlet=upstream
        self.post=None
    
    def __call__(self):
        if callable(self.post):
            self.post()

class shaft:
    # Tallies all torques.

    def __init__(self, name=None, cnxn=[]):
        self.station=None
        self.ID=name
        self.Nmech=1.0
        self.power=0.0
        self.torque=0.0
        self.mechEff=1.0
        self.refTrq=1.0
        self.omega=self.Nmech*np.pi/30.0
        self.pwrND=None
        self.callables=['mechEff']
        self.cnxn=cnxn
        self.Ptrq=0.0
        self.Ntrq=0.0
        self.ref=None
        self.STEADY=True
        for cnx in self.cnxn:
            cnx.shaft=self

    def __call__(self):
        self.ERR=1.0
        self.omega=self.Nmech*np.pi/30.0
        self.power=0
        self.Ptrq=0
        self.Ntrq=0
        self.torque=0
        for cnx in self.cnxn:
            pwr=cnx.power
            if pwr>0:
                self.Ptrq+=pwr*self.mechEff/self.omega
            else:
                self.Ntrq+=pwr/self.omega
        
        self.torque=self.Ptrq+self.Ntrq
        self.power=self.torque*self.omega
        
        if self.ref is None :
            if self.Ptrq>0:
                self.ref=self.Ptrq
            elif abs(self.Ntrq)>1:
                self.ref=abs(self.Ntrq)
            else:
                self.ref=1.0

        if self.STEADY:
            self.ERR=(self.torque/self.ref)

    
class load:

    def __init__(self, name=None):
        self.station=None
        self.ID=name
        self.Nmech=None
        self.kWe=0.0 # Electrical power of the load in kW.
        self.eff=1.0 # Electrical efficiency of the load.
        self.power=0.0
        self.shaft=None
        
    def __call__(self):
        # Sign of shaft load power is opposite of electrical power. 
        # If the load is negative, it is a starter and efficiency def is reciprocal.
        if self.kWe < 0.0: 
            self.power=-self.kWe*1000.0*self.eff
        else: 
            self.power=-self.kWe*1000.0/self.eff
        self.Nmech=self.shaft.Nmech

class BadDesign(Exception):
    pass
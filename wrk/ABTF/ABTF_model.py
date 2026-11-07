import cy23 as cy

"""
DECLARE MODEL
""" 
# Create model and name it.
eng=cy.Model('ABTF')

"""
ADD COMPONENTS
""" 
# Add all components, name them and declare connections.
amb=eng.add(cy.ambient, 'amb')
intk=eng.add(cy.duct,'intk', amb.outlet)
spl=eng.add(cy.splitter,'spl', intk.outlet)

lpc=eng.add(cy.compressor,'lpc',spl.outlet1)
lpbw=eng.add(cy.bleedSource,'lpbw', lpc.outlet)
lpbn=eng.add(cy.bleedSource,'lpbn', lpbw.outlet)
ess=eng.add(cy.duct, 'ess', lpbn.outlet)
swa=eng.add(cy.duct,'swa',ess.outlet)
hpc=eng.add(cy.compressor,'hpc',swa.outlet)
hpbw=eng.add(cy.bleedSource,'hpbw', hpc.outlet)
hpbn=eng.add(cy.bleedSource,'hpbn', hpbw.outlet)
comb=eng.add(cy.burner,'comb', hpbn.outlet)
hclw=eng.add(cy.bleedSink,'hclw',comb.outlet,hpbw)
hpt=eng.add(cy.turbine,'hpt',hclw.outlet)
hcln=eng.add(cy.bleedSink,'hcln',hpt.outlet,hpbn)
lclw=eng.add(cy.bleedSink,'lclw',hcln.outlet,lpbw)
lpt=eng.add(cy.turbine,'lpt',lclw.outlet)
lcln=eng.add(cy.bleedSink,'lcln',lpt.outlet,lpbn)

fan=eng.add(cy.compressor, 'fan',spl.outlet2)
bpd=eng.add(cy.duct,'bpd',fan.outlet)

mxr=eng.add(cy.mixer, 'mxr', lcln.outlet, bpd.outlet)
abrnr=eng.add(cy.burner, 'abrnr', mxr.outlet)
noz=eng.add(cy.nozzle,'noz',abrnr.outlet)
fe=eng.add(cy.flowEnd,'fe',noz.outlet)

lpshft=eng.add(cy.shaft, 'lpshft', [fan, lpc, lpt])
hpshft=eng.add(cy.shaft, 'hpshft', [hpc, hpt])

"""
INITIALIZE VARIABLES
"""
# Initial variables are needed to get solution started. Always good to declare defaults.
amb.W=10.0
amb.P=24.54*1E3
amb.T=220.0
amb.M=0.85
fan.eff=0.87
lpc.eff=0.87

hpc.eff=0.87
hpt.eff=0.90

lpt.eff=0.90
lpshft.mechEff=hpshft.mechEff=0.99
comb.dPqP=0.05
abrnr.dPqP=0.07
swa.dPqP=0.02
comb.eff=0.99
abrnr.eff=0.99
noz.Pback=amb.P
eng.params.FN_lbf=0
# mxr.M1=mxr.M2=0.3 # Can be used for mixer matching calc. Otherwise total pressures used.


noz.CONDI=False # When sizing the nozzle, do not allow con-di nozzle.


lpbn.ref=lpbw.ref=hpbn.ref=hpbw.ref=lpc # All references for bleed flow fractional flow definitions are the same: lpc

lpbw.fracW=0.03
lpbn.fracW=0.01

hpbw.fracW=0.10
hpbn.fracW=0.01

fan.Nmap=lpc.Nmap=hpc.Nmap=0.95

lpt.Nmap=hpt.Nmap=0.875
lpt.ERmap=hpt.ERmap=9.0

comb.FAR=0.02
abrnr.FAR=0.02
"""
LOAD MAPS
"""

# Load pickled maps. These were pickled from excel data using a separate script.
fan.MAP=lpc.MAP=hpc.MAP=cy.loadmap(r'.\genericCmap.pkl')
lpt.MAP=hpt.MAP=cy.loadmap(r'.\genericTmap.pkl')

"""
STATION NUMBERS (SAE ARP 755A style)
Only declaring those of interest. cy23 assumes the component inlet as the station number.
"""
amb.station    = 0   # ambient "component" represents freestream
intk.station  = 1   # connects to amb outlet (freestream)
lpc.station   = 2   # LPC inlet / core inlet
hpc.station   = 26   # LPC inlet / core inlet
comb.station  = 4   # combustor inlet
hpt.station   = 41  # HPT inlet
lpt.station   = 44  # LPT inlet
fan.station   = 12  # fan inlet (bypass stream side)
bpd.station   = 16  # bypass duct inlet
mxr.station   = 5   # core stream inlet to mixer (LPT exit)
# Afterburner
abrnr.station = 6   # afterburner inlet
# Nozzle and jet
noz.station   = 7  # nozzle inlet
fe.station    = 9   # jet plane used by flight effects

# Now that the model is declared, it can be saved in a pickle for loading elsewhere, rather than 
# using module imports. 
cy.saveModel(eng, 'ABTF.pkl')

import cy23 as cy

"""
DECLARE MODEL
""" 
# Create model and name it.
eng=cy.Model('SSGG')

"""
ADD COMPONENTS
""" 
# Add all components, name them and declare connections.
amb=eng.add(cy.ambient, 'amb')
intk=eng.add(cy.duct,'intk', amb.outlet)

# Compressor Group 1 
c1=eng.add(cy.compressor,'c1',intk.outlet)
bld1W=eng.add(cy.bleedSource,'bld1W', c1.outlet)
bld1=eng.add(cy.bleedSource,'bld1', bld1W.outlet)
hdl1=eng.add(cy.bleedSource,'hdl1', bld1.outlet)

# Compressor Group 2 
c2=eng.add(cy.compressor,'c2',hdl1.outlet)
bld2W=eng.add(cy.bleedSource,'bld2W', c2.outlet)
bld2=eng.add(cy.bleedSource,'bld2', bld2W.outlet)
hdl2=eng.add(cy.bleedSource,'hdl2', bld2.outlet)

# Compressor Group 3 
c3=eng.add(cy.compressor,'c3',hdl2.outlet)
bld3W=eng.add(cy.bleedSource,'bld3W', c3.outlet)
bld3N=eng.add(cy.bleedSource,'bld3N', bld3W.outlet)
bld3=eng.add(cy.bleedSource,'bld3', bld3N.outlet)
hdl3=eng.add(cy.bleedSource,'hdl3', bld3.outlet)

comb=eng.add(cy.burner,'comb', hdl3.outlet)

preNGV=eng.add(cy.bleedSink,'preNGV',comb.outlet,bld3N)
postNGV=eng.add(cy.bleedSink,'postNGV',preNGV.outlet,[bld1W, bld2W, bld3W])
turb=eng.add(cy.turbine,'turb',postNGV.outlet)
postTurb=eng.add(cy.bleedSink,'postTurb',turb.outlet,[bld1, bld2, bld3])

exh=eng.add(cy.duct, 'exh',postTurb.outlet)
diff=eng.add(cy.diffuser, 'diff',exh.outlet, amb)
fe=eng.add(cy.flowEnd,'fe',diff.outlet)

gen=eng.add(cy.load, 'gen')

shft=eng.add(cy.shaft, 'shft', [c1, c2, c3, turb, gen])


"""
INITIALIZE VARIABLES
"""
# Initial variables are needed to get solution started. Always good to declare defaults.
amb.W=100.0
amb.P=101325
amb.T=288.15
amb.M=0.0

gen.kWe=50E6

c1.eff=0.87
c2.eff=0.87
c3.eff=0.87

turb.eff=0.90
shft.mechEff=0.99
comb.dPqP=0.05
exh.dPqP=0.02
comb.eff=0.99
diff.cpr=0.6


# All references for bleed flow fractional flow definitions are the same: c1
bld1.ref=bld1W.ref=hdl1.ref=c1
bld2.ref=bld2W.ref=hdl2.ref=c1
bld3.ref=bld3W.ref=bld3N.ref=hdl3.ref=c1

bld1.fracW=0.00
bld1W.fracW=0.0

bld2.fracW=0.0
bld2W.fracW=0.0

bld3.fracW=0.0
bld3W.fracW=0.0
bld3N.fracW=0.0

c1.Nmap=c2.Nmap=c3.Nmap=0.95

turb.Nmap=0.875
turb.ERmap=9.0

comb.FAR=0.02

"""
LOAD MAPS
"""

# Load pickled maps. These were pickled from excel data using a separate script.
c1.MAP=c2.MAP=c3.MAP=cy.loadmap(r'.\genericCmap.pkl')
turb.MAP=cy.loadmap(r'.\genericTmap.pkl')

# Now that the model is declared, it can be saved in a pickle for loading elsewhere, rather than 
# using module imports. 
cy.saveModel(eng, 'SSGG.pkl')

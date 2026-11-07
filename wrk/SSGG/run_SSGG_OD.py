import cy23 as cy

"""
LOAD EXISTING ENGINE MODEL
"""

eng = cy.loadModel("SSGG.pkl")

"""
DECLARE DEPENDENTS & INDEPENDENTS
"""
# Shaft torque match error
shftERR=cy.Dependent(eng.shft, 'ERR', 0.0)

# Turbine flow match error
turbERR=cy.Dependent(eng.turb, 'ERR', 0.0)

# Compressor flow match errors
c1ERR=cy.Dependent(eng.c1, 'ERR', 0)
c2ERR=cy.Dependent(eng.c2, 'ERR', 0)
c3ERR=cy.Dependent(eng.c3, 'ERR', 0)

# TET match error 
tetTarg=cy.Dependent(eng.comb, 'ERR', 0.0)

# Diffuser back pressure match error.
diffERR=cy.Dependent(eng.diff, 'ERR', 0)

# Independents
turbER=cy.Independent(eng.turb, 'ER')
flow=cy.Independent(eng.amb, 'W')
c1beta=cy.Independent(eng.c1, 'beta')
c2beta=cy.Independent(eng.c2, 'beta')
c3beta=cy.Independent(eng.c3, 'beta')
nmech=cy.Independent(eng.shft, 'Nmech')
far=cy.Independent(eng.comb, 'FAR')
pwr=cy.Independent(eng.gen, 'kWe')

"""
DEFINE SOLVERS
"""

#fuelMatch: Match the fuel rate directly.
fuelMatch=cy.Solver(eng) # Create a solver for model eng.
fuelMatch.offdesign=True # This is an off-design solver, so it will use maps and all areas are fixed.
fuelMatch.indeps=[flow, c1beta, c2beta, c3beta, pwr, turbER]
fuelMatch.deps=[diffERR, c1ERR, c2ERR, c3ERR, shftERR, turbERR]

#Power Match: match desired power at fixed RPM, TET is an output.
pwrmatch=cy.Solver(eng) # Create a solver for model eng.
pwrmatch.offdesign=True # This is an off-design solver, so it will use maps and all areas are fixed.
pwrmatch.indeps=[flow, c1beta, c2beta, c3beta, turbER, far]
pwrmatch.deps=[diffERR, c1ERR, c2ERR, c3ERR, shftERR, turbERR]


"""
SET SCALARS AND BLEEDS
"""

eng.bld1.fracW=0.01
eng.bld1W.fracW=0.005

eng.bld2.fracW=0.02
eng.bld2W.fracW=0.02

eng.bld3.fracW=0.06
eng.bld3W.fracW=0.07
eng.bld3N.fracW=0.02

eng.c1.AuditScaleWc=0.98
eng.c1.AuditAddEff=-0.01

eng.c2.AuditScaleWc=1.03
eng.c2.AuditAddEff=-0.026

eng.c3.AuditScaleWc=0.99
eng.c3.AuditAddEff=0.005

eng.turb.AuditScaleWc=0.975
eng.turb.AuditAddEff=0.003

"""
RUN OFF DESIGN
"""
# Use batch run to run multiple off-design cases. 
# First we define the outputs desired.
# The BatchRun instance is declared specifying the solver, inputs and outputs.
# The batch run is run by calling the spreadsheet and tab name where the input values are
# declared. All output values are written to the right on the same tabe, including a
# timestamp and convergence status. The inputs are also written as outputs for confirmation.

outputs=["shft.Nmech", "gen.kWe", "c1.outlet.Pt","c1.outlet.Tt", "c2.outlet.Pt","c2.outlet.Tt", "c3.outlet.Pt","c3.outlet.Tt","exh.inlet.Tt"]

runTETs=cy.BatchRun(fuelMatch, ["amb.T", "amb.P", "comb.Wfuel_in"], outputs)

runTETs("calibModel.xlsx", "TET1750_wfactors")
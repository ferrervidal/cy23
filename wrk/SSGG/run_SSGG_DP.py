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

# Independent
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

# Deisgn point solver. For a given generator load at a given TET target, this calculates the required engine size. 
# Compressor PR's must be given. 
sizeEngine=cy.Solver(eng)
sizeEngine.offdesign=False # This is a design-point solver, so it calculates map scalars and areas are free to change.
sizeEngine.indeps=[turbER, far, flow]
sizeEngine.deps=[shftERR, tetTarg, diffERR]

# Diagnostics file to assess solver convergence.
sizeEngine.diag=True 

"""
RUN DESIGN POINT SIZING
"""

eng.comb.COT=1750.0 # Target TET.
eng.c1.PR=3 # Desired PR.
eng.c2.PR=2.5 # Desired PR.
eng.c3.PR=2.5 # Desired PR.
eng.gen.kWe=50E3 # Desired power output
eng.amb.W=140.0 # Guess initial mass flow
eng.turbER=4.0 # initial guess
eng.amb.P=101325
eng.amb.T=288.15
eng.diff.Mdes=0.6
eng.diff.cpr=0.6
eng.diff.AR=2.0

sizeEngine() # Run the solver.
cy.write(eng, 'SSGG') # Write variables from the "eng" model.

#Save model again after sizing, so it can be run in off-design with other scripts. 
cy.saveModel(eng, 'SSGG.pkl')

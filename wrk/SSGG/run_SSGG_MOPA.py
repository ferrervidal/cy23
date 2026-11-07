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

#Create a casewriter to save data.
crw=cy.CaseWriter(eng)


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

eng.c1.AuditScaleWc=1.0
eng.c1.AuditAddEff=0

eng.c2.AuditScaleWc=1
eng.c2.AuditAddEff=0

eng.c3.AuditScaleWc=1
eng.c3.AuditAddEff=0

eng.turb.AuditScaleWc=1
eng.turb.AuditAddEff=0

eng.c1.AuditEffScalar=1.0
eng.c2.AuditEffScalar=1.0
eng.c3.AuditEffScalar=1.0
eng.turb.AuditEffScalar=1.0

"""
RUN OFF DESIGN
"""

# Create a multiple-operating point gas path analysis calibrator. 
# The solver to use is specified along with measurements and the calibration factors
# (map scalars) to use. 
# MOPA uses linear least squares to determine the best set of calibration factors for a given data set. 
# MOPA load excel is used to load the test data used for calibration. 

mopaCal=cy.MOPA(
    solver=fuelMatch,
    inputs=["amb.T", "amb.P", "comb.Wfuel_in"],
    measurements=["gen.kWe","c1.outlet.Pt","c1.outlet.Tt","c2.outlet.Pt","c2.outlet.Tt","c3.outlet.Pt","c3.outlet.Tt","exh.inlet.Tt"],
    calFactors=[
        "c1.AuditScaleWc",
        "c1.AuditEffScalar",
        "c2.AuditScaleWc",
        "c2.AuditEffScalar",
        "c3.AuditScaleWc",
        "c3.AuditEffScalar",
        "turb.AuditScaleWc",
        "turb.AuditEffScalar",
        ]
)

mopaCal.load_excel("MOPAinput.xlsx", sheet_name=0)

mopaCal.calibrate(savename="calibrated_SSGG.pkl")


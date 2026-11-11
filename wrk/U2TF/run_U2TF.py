import cy23 as cy

"""
LOAD EXISTING ENGINE MODEL
"""

eng = cy.loadModel("U2TF.pkl")

"""
DECLARE DEPENDENTS & INDEPENDENTS
"""
# Dependents and indenpendents are created here.
# We need to know the name of the dep/indep variables within each component.
# For dependent values, the last number is the target value.
# Note that we may be using dependents to match a desired target, see
# for example "tetTarg" or "fnTarg" below, which target a specific TET or thrust.
# We can set the target value to a random value (say 0.0) and later updated via
# fnTarg.target etc.
# In the case of fnTarg, the FN_lbf variable was previously created in the model, but gets
# updated here in the performance fucntion.

lpshftERR=cy.Dependent(eng.lpshft, 'ERR', 0.0)
hpshftERR=cy.Dependent(eng.hpshft, 'ERR', 0.0)

lptERR=cy.Dependent(eng.lpt, 'ERR', 0.0)
hptERR=cy.Dependent(eng.hpt, 'ERR', 0.0)

tetTarg=cy.Dependent(eng.comb, 'ERR', 0.0)
fnTarg=cy.Dependent(eng.params, 'FN_lbf', 0.0)

nozERR=cy.Dependent(eng.noz, 'ERR', 0)
bpnERR=cy.Dependent(eng.bpn, 'ERR', 0)
fanERR=cy.Dependent(eng.fan, 'ERR', 0)
lpcERR=cy.Dependent(eng.lpc, 'ERR', 0)
hpcERR=cy.Dependent(eng.hpc, 'ERR', 0)

lptER=cy.Independent(eng.lpt, 'ER')
hptER=cy.Independent(eng.hpt, 'ER')

flow=cy.Independent(eng.amb, 'W')
bprind=cy.Independent(eng.spl, 'BPR')
fprind=cy.Independent(eng.fan, 'PR')
fanbeta=cy.Independent(eng.fan, 'beta')
lpcbeta=cy.Independent(eng.lpc, 'beta')
hpcbeta=cy.Independent(eng.hpc, 'beta')
n1=cy.Independent(eng.lpshft, 'Nmech')
n3=cy.Independent(eng.hpshft, 'Nmech')
far=cy.Independent(eng.comb, 'FAR')

"""
DEFINE SOLVERS
"""

tetmatch=cy.Solver(eng) # Create a solver for model eng.
tetmatch.offdesign=True # This is an off-design solver, so it will use maps and all areas are fixed.
tetmatch.indeps=[flow, fanbeta, lpcbeta, hpcbeta,n1, n3, lptER, hptER, far, bprind]
tetmatch.deps=[nozERR, fanERR, lpcERR, hpcERR, lptERR, hptERR, lpshftERR, hpshftERR, tetTarg, bpnERR]

thrustmatch=cy.Solver(eng) # Create a solver for model eng.
thrustmatch.offdesign=True # This is an off-design solver, so it will use maps and all areas are fixed.
thrustmatch.indeps=[flow, fanbeta, lpcbeta, hpcbeta,n1, n3, lptER, hptER, far, bprind]
thrustmatch.deps=[nozERR, fanERR, lpcERR, hpcERR, lptERR, hptERR, lpshftERR, hpshftERR, fnTarg, bpnERR]

sizeEngine=cy.Solver(eng)
sizeEngine.offdesign=False # This is a design-point solver, so it calculates map scalars and areas are free to change.
sizeEngine.indeps=[lptER, hptER, far, flow]
sizeEngine.deps=[lpshftERR, hpshftERR, tetTarg, fnTarg]

"""
PERFORMANCE CALCULATIONS
"""
# Model performance function.
eng.post=cy.Performance(eng,"jet")

"""
RUN MODEL
"""

#Create a casewriter to save data.
crw=cy.CaseWriter(eng)
prw=cy.ParamWriter(eng, ["amb.W", "params.FN_lbf", "noz.Ath", "params.TSFC_US", "comb.outlet.Tt", "lpc.PR", "hpc.PR", "spl.BPR"])

sizeEngine.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1750 # Target TET.
eng.lpc.PR=eng.fan.PR=1.5 # Desired PR.
eng.hpc.PR=26 # Desired PR.
eng.spl.BPR=10 # Initial guess BPR, to be varied by solver.
fnTarg.target=80000 # Target thrust, in [lbf] as the dependent checks against eng.params.FN_lbf.
sizeEngine() # Run the solver.
prw.update()
# Pickle can be used to save the state of the model. This includes all component variables. So now
# that the engine has been sized using the sizeEngine solver, we would load it in in a different file and
# run off-design. Also, if there is a particularly hard case to converge and we must proceed in small steps,
# we can save the last converged state using pickle as below.
# Note that the below only saves the model and its variables but it does not save the solvers or the performance
# function. A user may choose to define all solvers and performance function in a separate "solver definitio"
# module and simply import that, or write them all as part of the run file as done here.
cy.saveModel(eng, 'U2TF_80klbf.pkl')


# Now we run an off-design case at the same designed TET.
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1750 # Target TET, in [K].
tetmatch() # Run the solver.
prw.update()

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1500 # Target TET, in [K].
tetmatch() # Run the solver.
prw.update()

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1400 # Target TET, in [K].
tetmatch() # Run the solver.
prw.update()

# Now we match a thrust setting
thrustmatch.diag=True # Diagnostics file to assess solver convergence.
fnTarg.target=80000*0.75
thrustmatch() # Run the solver.
prw.update()


prw.write('U2TF_paramOutput')




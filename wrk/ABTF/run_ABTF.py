import cy23 as cy

"""
LOAD EXISTING ENGINE MODEL
"""

eng = cy.loadModel("ABTF.pkl")

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
abTarg=cy.Dependent(eng.abrnr, 'ERR', 0.0)
fnTarg=cy.Dependent(eng.params, 'FN_lbf', 0.0)

nozERR=cy.Dependent(eng.noz, 'ERR', 0)
mxrERR=cy.Dependent(eng.mxr, 'ERR', 0)
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
farAB=cy.Independent(eng.abrnr, 'FAR')

"""
DEFINE SOLVERS
"""

tetmatch=cy.Solver(eng) # Create a solver for model eng.
tetmatch.offdesign=True # This is an off-design solver, so it will use maps and all areas are fixed.
tetmatch.indeps=[flow, fanbeta, lpcbeta, hpcbeta,n1, n3, lptER, hptER, far, farAB, bprind]
tetmatch.deps=[nozERR, fanERR, lpcERR, hpcERR, lptERR, hptERR, lpshftERR, hpshftERR,tetTarg, abTarg, mxrERR]


tetmatchABoff=cy.Solver(eng)
tetmatchABoff.offdesign=True
tetmatchABoff.indeps=[flow, fanbeta, lpcbeta, hpcbeta,n1, n3, lptER, hptER, far, bprind]
tetmatchABoff.deps=[nozERR, fanERR, lpcERR, hpcERR, lptERR, hptERR, lpshftERR, hpshftERR,tetTarg, mxrERR]

sizeEngine=cy.Solver(eng)
sizeEngine.offdesign=False # This is a design-point solver, so it calculates map scalars and areas are free to change.
sizeEngine.indeps=[lptER, hptER, far, farAB, flow, bprind]
sizeEngine.deps=[lpshftERR, hpshftERR, tetTarg, abTarg, fnTarg, mxrERR]

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
prw=cy.ParamWriter(eng, ["amb.W", "params.FN_lbf", "noz.Ath", "params.TSFC_US"])

sizeEngine.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1750 # Target TET.
eng.abrnr.COT=2000 # Target after-burner firing temperature.
eng.lpc.PR=eng.fan.PR=3 # Desired PR.
eng.hpc.PR=20 # Desired PR.
eng.spl.BPR=0.5 # Initial guess BPR, to be varied by solver.
fnTarg.target=10000 # Target thrust, in [lbf] as the dependent checks against eng.params.FN_lbf.
sizeEngine() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model to an excel file called "test".
# Pickle can be used to save the state of the model. This includes all component variables. So now
# that the engine has been sized using the sizeEngine solver, we would load it in in a different file and
# run off-design. Also, if there is a particularly hard case to converge and we must proceed in small steps,
# we can save the last converged state using pickle as below.
# Note that the below only saves the model and its variables but it does not save the solvers or the performance
# function. A user may choose to define all solvers and performance function in a separate "solver definitio"
# module and simply import that, or write them all as part of the run file as done here.
cy.saveModel(eng, 'ABTF_10klbf.pkl')


# Now we run an off-design case at the same designed TET.
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1750 # Target TET, in [K].
tetmatch() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1745# Target TET, in [K].
tetmatch() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1700# Target TET, in [K].
eng.abrnr.COT=1990 # Target after-burner firing temperature.
tetmatch() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1650 # Target TET, in [K].
eng.abrnr.COT=1700 # Target after-burner firing temperature.
tetmatch() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1650 # Target TET, in [K].
eng.abrnr.COT=1600 # Target after-burner firing temperature.
tetmatch() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1650 # Target TET, in [K].
eng.abrnr.COT=1200 # Target after-burner firing temperature.
tetmatch() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

# Now we run an off-design case at lower TET
tetmatch.diag=True # Diagnostics file to assess solver convergence.
eng.comb.COT=1600 # Target TET, in [K].
eng.abrnr.COT=1000 # Target after-burner firing temperature.
tetmatch() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

# Now we run an off-design case at the same designed TET but no AB.
tetmatchABoff.diag=True # Diagnostics file to assess solver convergence.
eng.abrnr.OFF=True 
eng.noz.Ath=eng.noz.Ath*0.4822
eng.comb.COT=1500 # Target TET, in [K].
tetmatchABoff() # Run the solver.
prw.update()
#crw.update()
#write(eng, 'test') # Write variables from the "eng" model.

#crw.write('exampleOutput')
prw.write('ABTF_paramOutput')




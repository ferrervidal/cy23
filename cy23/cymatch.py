"""
CY23.CYMATCH

This module provides classes and methods to enable solving system models. 

Key Components:
---------------
1. **Model**: A class to represent a system model composed of various components. It includes methods to add components, 
   set off-design operation mode, and execute the model.
2. **Independent** and **Dependent**: Classes to define the independent and dependent variables that are used by the solver 
   for system optimization. These classes ensure that the solver computes the system state by comparing dependent variables 
   against target values.
3. **Solver**: A solver class that uses non-linear equation solvers to find solutions for a given model's variables. 
   The solver computes independent and dependent variables, checks for convergence, and writes diagnostics.

Features:
---------
- **Modeling**: Users can define a system by assembling various components, where each component may have its own attributes 
  and parameters.
- **Solver Setup**: The solver can handle independent and dependent variables, calculate errors, and determine convergence.
- **Diagnostics**: The solver can write detailed diagnostic information to a file during the optimization process.

Usage:
------
This module can be used to solve models composed of cy23.cycomps components. 
"""

# Importing necessary libraries and modules
from scipy.optimize import fsolve  # For solving non-linear equations
from scipy.optimize import least_squares
from scipy.interpolate import interp1d  # For 1D interpolation
import pandas as pd  # For data manipulation and file reading
import numpy as np  # For numerical operations
from datetime import datetime  # For date and time handling
import os  # For file operations
import pickle  # For object serialization and deserialization
from cy23.cycomps import *  # Importing components from a custom module
import cy23.cymaps  # Importing additional module (assumed custom)
from typing import Callable, List, Optional, Tuple, Union  # For type annotations
from tabulate import tabulate # For MOPA output.

# Main Model class representing a system model composed of components
class Model:
    """
    A class to represent a system model composed of various components.
    
    Attributes:
    -----------
    name : str
        Name of the model.
    comps : List
        List of components in the model.
    _OFFDESIGN : bool
        Flag to indicate off-design mode.
    """

    def __init__(self, name: str = 'CYCLE'):
        """
        Initializes the model with a name and empty component list.
        
        Parameters:
        -----------
        name : str
            Name of the model (default is 'CYCLE').
        """
        self.name = name
        self.comps: List[object] = []  # List to store components
        self._OFFDESIGN: bool = False  # Flag for off-design mode
        self.status = 'Not run'  # Status of the model run
        self.params = Container('params')  # Container for model parameters
        self.post=[] # List of functions or function

    @property
    def OFFDESIGN(self) -> bool:
        """Property to get the off-design status."""
        return self._OFFDESIGN

    @OFFDESIGN.setter
    def OFFDESIGN(self, mode: bool):
        """
        Sets the off-design mode for the model and its components.
        
        Parameters:
        -----------
        mode : bool
            Off-design mode flag.
        """
        self._OFFDESIGN = bool(mode)
        # Propagate the mode to all components in the model
        for comp in self.comps:
            if hasattr(comp, 'OFFDESIGN'):
                comp.OFFDESIGN = self._OFFDESIGN

    def __call__(self):
        """
        Executes all components in the model.
        """
        for comp in self.comps:
            comp()  # Call each component's callable method.
        if self.post is not None:
            if callable(self.post):
                self.post()
            elif isinstance(self.post, list):
                for postFunc in self.post:
                    if callable(postFunc):
                        postFunc()

    def postSolve(self):
        """
        Executes all post-solve component functions if they exist.
        These are component-level functions that should be executed only after the model is converged.
        """
        for comp in self.comps:
            compPostSolve=getattr(comp, "__postSolve__", None)
            if callable(compPostSolve):
                compPostSolve()

    def add(self, comp: Callable, name: str, *args, **kwargs) -> object:
        """
        Adds a component to the model with a unique name.

        Parameters:
        -----------
        comp : Callable
            The component class to be added.
        name : str
            The name of the component.
        *args, **kwargs : 
            Additional arguments to initialize the component.

        Returns:
        --------
        obj : object
            The instantiated component object.
        """
        original_name = name
        i = 1
        # Ensure unique component names within the model
        while hasattr(self, name):
            name = f"{original_name}_{i}"
            i += 1
        obj = comp(name, *args, **kwargs)  # Create the component instance
        setattr(self, name, obj)  # Set the component as an attribute of the model
        self.comps.append(obj)  # Add the component to the list
        return obj
    
    def post(self):
        """Placeholder post function."""
        pass


# Class to define independent variables for a solver
class Independent:
    """
    Class to define the independent variables for solver computations.
    """
    def __init__(self, obj, attr):
        self.obj = obj  # Object containing the attribute
        self.attr = attr  # Attribute name as a string

        if not isinstance(self.attr, str):
            raise SolverSetupError(f'ERROR: Attribute should be a string. Got {type(self.attr)} for object {self.obj}')


# Class to define dependent variables and match targets for a solver
class Dependent:
    """ 
    Class to define dependent variables and their matching targets. 
    """
    def __init__(self, obj, attr, targ=None):
        self.obj = obj  # Object containing the attribute
        self.attr = attr  # Attribute name as a string
        self.target = targ  # Target value or callable for comparison
        self.ERR = 1.0  # Default error value
        
        if not isinstance(self.attr, str):
            raise SolverSetupError('ERROR: For matching targets, the attribute should be defined as a string.')
    
    def __call__(self):
        # Get the actual value of the attribute
        actVal = getattr(self.obj, self.attr)
        
        # Determine the target value based on its type
        if isinstance(self.target, tuple):
            targetVal = getattr(self.target[0], self.target[1])
        elif callable(self.target):
            targetVal = self.target()
        else:
            targetVal = self.target
            
        if callable(targetVal):
            targetVal = targetVal()
        
        # Calculate the relative error
        self.ERR = ((actVal - targetVal) / max(abs(targetVal), 1.0)) * 100


class Solver:
    """
    A solver for a model that computes solutions for independent and dependent variables,
    and checks convergence against a set tolerance.

    Attributes:
    -----------
    model : object
        The model to be solved.
    offdesign : bool
        Flag indicating if the system is off-design.
    indeps : list
        List of independent variables.
    deps : list
        List of dependent variables.
    msg : str
        Message about the solver status.
    initState : list
        Initial state of the system.
    diag : bool
        Flag to indicate if diagnostics should be written.
    state : np.ndarray
        Current state of the system.
    goodState : np.ndarray
        Best state after convergence.
    converged : bool
        Flag indicating if the solver has converged.
    tol : float
        Tolerance for convergence.
    maxfev : int
        Maximum number of function evaluations for the solver.
    perturb : float
        Perturbation value for numerical derivatives.
    rewriteDiag : bool
        Flag to indicate if diagnostics should be rewritten.
    ERRRMS : float
        Root mean square error.
    MAXERR : float
        Maximum error.
    post : Callable, optional
        Optional post-solution function to be called after solving.
    """

    def __init__(self, model):
        """
        Initializes the Solver with a model and sets default attributes.

        Parameters:
        -----------
        model : object
            The model to be solved.
        """
        self.model = model
        self.offdesign = False
        self.indeps = []
        self.deps = []
        self.msg = 'Solver not yet run.'
        self.initState = []
        self.diag = False
        self.state = None
        self.goodState = None
        self.converged = False
        self.tol = 1.5E-8
        self.maxfev = 500
        self.perturb = self.tol * 100
        self.rewriteDiag = True
        self.ERRRMS = 1.0
        self.MAXERR = 1.0
    
   

    def objFunc(self, sol):
        """
        Objective function to evaluate the error between the model's dependent variables and their target values.

        Parameters:
        -----------
        sol : np.ndarray
            The solution vector representing the state of independent variables.

        Returns:
        --------
        np.ndarray
            The error values for the dependent variables.
        """
        self.state = np.multiply(sol, self.initState)  # Apply solution scaling
        for i, indep in enumerate(self.indeps):
            setattr(indep.obj, indep.attr, self.state[i])  # Set the value of independent variables in the model
        self.model()  # Update model with new state
        ERRS = self.getERRS()  # Calculate the errors
        self.ERRRMS = np.sqrt(np.mean(ERRS**2))  # RMS error
        self.MAXERR = max(abs(ERRS))  # Maximum error
        if self.diag:
            self.writeDiag(ERRS)  # Write diagnostics if enabled

        return ERRS

    def getState(self):
        """
        Retrieves the current state of the independent variables.

        Returns:
        --------
        np.ndarray
            The current state values for the independent variables.
        """
        if self.converged or self.goodState is None or len(self.goodState) != len(self.indeps):
            state = np.zeros((len(self.indeps)))
            for i, indep in enumerate(self.indeps):
                state[i] = getattr(indep.obj, indep.attr)  # Get state from the model
        else:
            state = self.goodState  # Use stored good state if available
        return state

    def getERRS(self):
        """
        Calculates the errors between the model's predicted and actual dependent variable values.

        Returns:
        --------
        np.ndarray
            The error values for each dependent variable.
        """
        ERRS = np.ones((len(self.deps)))  # Initialize error array
        for i, dep in enumerate(self.deps):
            dep()  # Call the dependent variable's function to update the error
            ERRS[i] = dep.ERR  # Store the error
        return ERRS

    def initialize(self):
        """
        Initializes the solver by setting up the model, calculating initial states, and preparing for the solution.
        """
      
        self.model.OFFDESIGN = self.offdesign  # Set offdesign flag in model
        self.model()  # Update model
        self.state = self.getState()  # Get initial state
        self.initState = self.state  # Store initial state
        for i, init in enumerate(self.initState):
            if abs(init) < 10 * self.tol:
                self.initState[i] = 0.01  # Avoid zero initialization for small values

    def __call__(self, verbose=True):
        """
        Runs the solver to find the solution by calling the optimization solver.
        """
        if len(self.indeps) != len(self.deps):
            raise SolverSetupError('Number of independents and dependents is not the same.')  # Check for setup consistency
        self.initialize()  # Initialize solver
        guess = np.ones_like(self.initState) * 1.00  # Initial guess for independent variables
        sol, _, ier, self.msg = fsolve(self.objFunc, guess, full_output=True, epsfcn=self.perturb, maxfev=self.maxfev, factor=0.1, xtol=self.tol)  # Solve the system using fsolve
        self.state = np.multiply(sol, self.initState)  # Scale solution back to state
        self.objFunc(sol)  # Call objective function with solution
        if verbose:
            print(self.msg)  # Print solver status message
        if ier==1:
            self.model.status = "Converged" # Update model status
            self.model.postSolve() # Run the model post-solve function. 
        else:
            self.model.status = "Not converged" # Update model status
        self.converged = ier == 1  # Check if the solver converged (ier=1 means success)
        if not np.isnan(self.state).any():
            self.goodState = self.state  # Store good state if solution is valid
        if self.diag:
            with open("diag.txt", 'a') as f:  # Write diagnostics to a file
                f.write(' \n')
                f.write('######----->' + self.msg + '\n')
                f.write(' \n')
                f.write('VVVVVV --- NEXT POINT --- VVVVV\n')

    def writeDiag(self, ERRS):
        """
        Writes diagnostic information to a file.

        Parameters:
        -----------
        ERRS : np.ndarray
            The error values for the dependent variables.
        """
        inds = []  # List to store independent variable info
        deps = []  # List to store dependent variable info

        for indep in self.indeps:
            if indep.obj.ID is not None:
                ID = indep.obj.ID 
            else:
                ID ="No ID"
            indstr = ID + '|' + indep.attr
            inds.append(indstr)
        for dep in self.deps:
            if dep.obj.ID is not None:
                ID= dep.obj.ID
            else:
                ID = "No ID"
            depstr = ID + '|' + dep.attr
            deps.append(depstr)

        if not os.path.exists("diag.txt") or self.rewriteDiag:
            append_write = 'w'  # Overwrite if file doesn't exist or should be rewritten
        else:
            append_write = 'a'  # Append if file already exists
        with open("diag.txt", append_write) as f:
            f.write('----------------------\n')
            f.write('----------------------\n')
            for n in range(len(inds)):
                f.write("{: <10}\t {: <20f}\t {: <10}\t {: <20.4e}\n".format(inds[n], self.state[n], deps[n], ERRS[n]))
            strn = "ERROR RMS: %1.3E \n" % self.ERRRMS
            f.write('----------------------\n')
            f.write(strn)

        self.rewriteDiag = False  # Disable rewriting after the first write

class MOPA:
    """ MOPA Least-Squares Solver"""

    def __init__(self, solver, inputs, measurements, calFactors):
        self.solver=solver
        self.inputs=inputs # list like in batchsolve
        self.meas=measurements # list of strings with variable names 
        self.meas_start_col=len(self.inputs)
        self.calFactors=calFactors # The calibration factors or health parameters we want to determine.
        self.finalCalValues=None
        self.result=None
     
    def load_excel(self, filename, sheet_name=0):
        """
        Load operating-point data from Excel.

        The first len(self.inputs) columns are model inputs.
        The remaining columns are measured parameters.
        """

        df = pd.read_excel(filename, sheet_name=sheet_name)

        # Remove completely blank rows, if any.
        df = df.dropna(how="all")
        # Remove completely blank columns, if any.
        df = df.dropna(axis=1, how="all")

        expected_columns = len(self.inputs) + len(self.meas)

        if df.shape[1] != expected_columns:
            raise ValueError(
                f"Excel file contains {df.shape[1]} columns, but "
                f"{expected_columns} are required: "
                f"{len(self.inputs)} inputs and {len(self.meas)} measurements."
            )

        try:
            # Convert the DataFrame to a floating-point NumPy array.
            self.data = df.apply(pd.to_numeric, errors="raise").to_numpy(
                dtype=float
            )
        except ValueError as error:
            raise ValueError(
                "The Excel data contains nonnumeric values in a required column."
            ) from error

        self.points = self.data.shape[0]
        self.res = np.zeros((self.points, len(self.meas)))

        # Keeping the DataFrame can be useful for checking headers later.
        self.dataframe = df

        return self.data

    def _get_by_path(self, path: str):
        """
        Resolve a dotted path. 
        """
        obj = self.solver.model
        for name in path.split("."):
            if hasattr(obj, name):
                obj = getattr(obj, name)
            else:
                raise AttributeError(f"Could not resolve path {path!r}")
                break
        return obj

    def _set_by_path(self, path: str, value):
        """
        Resolve a dotted path and assign value. 
        """
        parent_path, _, attr = path.rpartition(".")
        obj = self._get_by_path(parent_path)
        if hasattr(obj, attr):
            setattr(obj,attr,value)
        else:
            print(f"WARNING:Could not resolve path {path!r}")

    def residuals(self,calValues):
        print("MOPA Pass...")
        residuals=np.zeros_like(self.res)
        for n, calValue in enumerate(calValues):
            self._set_by_path(self.calFactors[n], calValue)
        for point in range(self.points):
            # For each operating points in data array.
            for  col, inputt in enumerate(self.inputs):
                # For each input value, set model valuel
                value=self.data[point,col]
                self._set_by_path(inputt,value)

            # After all inputs set, run model solver
            self.solver(verbose=False)
            
            for col, target in enumerate(self.meas):
                # For each measurement, calculate model error, normalized as required.
                measValue=self.data[point,col+self.meas_start_col]
                modelValue=self._get_by_path(target)
                residuals[point,col]=100*(modelValue-measValue)/modelValue

            self.res=residuals

        return residuals.ravel()
    
    def solve(self):
        x0 = np.ones(len(self.calFactors))
        for n,calFactor in enumerate(self.calFactors):
            x0[n]=self._get_by_path(calFactor)

        result = least_squares(self.residuals, x0, diff_step=1e-4)
        self.finalCalValues=result.x.copy()
        if result.success:
            table_data=list(zip(self.calFactors, self.finalCalValues))
            headers=["Factor", "Value"]
            print("----------")
            print("MOPA Successful.")
            print("Cost: ", result.cost)
            print("----------")
            print(tabulate(table_data, headers=headers, tablefmt="fancy_grid"))
        else:
            print("----------")
            print("MOPA FAILED.")
            print("----------")
            print(result)
        return result.success

    def calibrate(self, save=True, savename="calibrated_model.pkl"):

        if self.finalCalValues is None:
            print("MOPA solve not yet run. Running MOPA Solve...")
            ok=self.solve()
        if not ok:
            print("MOPA FAILED. Did not calibrate.")
        else:
            print("MOPA successful. Calibrating model.")
            for n, calValue in enumerate(self.finalCalValues):
                self._set_by_path(
                    self.calFactors[n],
                    calValue,
                )
            if save:
                print("Saving calibrated model...")
                with open(savename, "wb") as f:
                    pickle.dump(self.solver.model, f)
                    print("Calibrated model saved as: ",savename)
        
    def calcIdent(self):

        calValues = np.ones([len(self.calFactors)])
        for n, val in enumerate(self.calFactors):
            calValues[n]=self._get_by_path(val)

        result = least_squares(
            self.residuals,
            calValues,
            jac="2-point",
            max_nfev=1,
        )

        self.J=result.jac
        self.svd= np.linalg.svd(result.jac, full_matrices=False, compute_uv=False)
        self.condition=np.linalg.cond(result.jac)
        self.norms=np.linalg.norm(result.jac, axis=0)


    

class SolverSetupError(Exception):
    """Custom exception for solver setup errors."""
    pass

class Container:
    """
    A simple container class.
    
    Attributes:
    -----------
    ID : str
        The identifier for the container.
    """
    
    def __init__(self, name):
        """
        Initializes the container with a given name.

        Parameters:
        -----------
        name : str
            The ID to assign to the container.
        """
        self.ID = name

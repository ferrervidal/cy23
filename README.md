# cy23
---
A toolset for thermodynamic cycle modeling.
Contact: luis.ferrer-vidal@erau.edu
---

### Description

cy23 is a set of tools useful for thermodynamic cycle modeling. It is a development research tool and thus not presently geared towards production use. It combines cycle matching capability alongside other general thermodynamic tools. cy23 continues to expand as need arises and thus grows organically in support of other projects. Therefore, cy23 might not include all functionality typically included in production cycle codes. 

Key features of the toolset: 
- Thermally perfect gas model for air and combustion products.
- Tools for the creation of callable turbomachinery performance maps from loaded spreadsheet data. 
- Completely customizable engine architectures, bleed arrangements, and solver setups for steady-state engine matching.
- Ability to save model states.
- Data output to formatted spreadsheets.
- Classes to enable automatic Multi-Operating Point Gas Path Analysis (MOPA).

### License
cy23 is distributed under a BSD 3-Clause License. For details refer to LICENSE.

### Source
cy23 is currently developed within the Performance & Turbomachinery group at Embry-Riddle Aeronautical University's Gas Tubine Laboratory. 

### Installation in developer mode

From the project root, install cy23 in editable mode:

```bash
python -m pip install -e .
```

This installs the required dependencies and allows changes to the local source code to take effect immediately.


### Usage & Examples
The cy23 folder contains the core modules.

To create a performance model, a .py file is created detailing the engine architecture and solver set-up. Maps should be available as .pkl files, which can be generated from properly formatted .csv data using the processCmap and processTmap functions within the cymaps modules. A separate .py "run file" is typically used to run models. Please refer to examples and the User guide in the wrk folder to best learn how to use this package.

Please refer to the "wrk" folder to find example models. These examples illustrate how to create and run models using most of the utilities available. The example models provided can be copied and modified to generate any desired architecture. 

The models included are:

-ABTF: After-burning turbofan. 
    This model illustrates the basic workflow as well as writing output using the basic, case and param writers.

-SSGG: Single-shaft gas generator.
    This model illustrates the basic workflow including how to model compressor interstage bleeds by splitting the compressor in groups. It also illustrates use of the BatchRun utility as well as the MOPA tool, to calibrate models to test data via multiple-operating point gas path analysis. 


General workflow to create a model:

- Use the makemaps.py script to import csv map data (see wkr folder for example) and create a pickle of the maps.
- Write a script to instantiate a "Model" object. Instantiate the model components, link them, and load the pickled maps. Finally, pickle the Model object.
- Write a run file where the pickled model object is loaded, create dependents, independents, and solvers. You can create auxiliary functions here too. 
- In the run file, write the run commands as desired. Note that you can save the model state at any time by pickling the model.
- Output data to a spreadsheet using one of the write utilities. 

##########################

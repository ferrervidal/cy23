"""
Top-level API for the cy23 package.

Usage:

    import cy23 as cy

    eng = cy.Model("MyEngine")
    amb = cy.ambient("amb")
    fan = cy.compressor("fan", amb.outlet)

    cy.write(eng, "exampleOutput")

Submodules (cymatch, cycomps, cymaps, cyutil) are still accessible as
cy23.cymatch, cy23.cycomps, etc.
"""

# Core modeling and solver
from .cymatch import (
    Model,
    Solver,
    Independent,
    Dependent,
    MOPA,
)

# Components
from .cycomps import (
    component,
    ambient,
    burner,
    compressor,
    duct,
    turbine,
    bleedSource,
    bleedSink,
    mixer,
    nozzle,
    splitter,
    diffuser,
    flowEnd,
    shaft,
    load,
)

# Map utilities
from .cymaps import (
    processCmap,
    processTmap,
    loadmap,
)

# Generic utilities and I/O
from .cyutil import (
    write,
    Performance,
    Table,
    DataSet,
    CaseWriter,
    ParamWriter,
    BatchRun,
    loadModel,
    saveModel,
)

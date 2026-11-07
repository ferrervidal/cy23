"""
cyutil

Utility classes and functions for the cy23 gas turbine performance code.

Contents
--------
Table
    One dimensional lookup table based on CSV data for maps and schedules.

DataSet
    Convenience wrapper for loading numeric data from CSV or Excel files.

loadmap
    Helper for loading pickled performance maps.

write
    Single point writer that exports a full model state to a formatted Excel sheet.

CaseWriter
    Case row writer that flattens the full model state over multiple runs
    into a single wide table (one row per case).

ParamWriter
    Lightweight case row writer that records only selected scalar parameters
    over multiple runs.

saveModel
    Save model and state to pickle file.

loadModel
    Load model and state.
"""

import numpy as np
import pandas as pd
from datetime import datetime
from typing import Optional, Any
from pathlib import Path
import pickle
from openpyxl import load_workbook
from cy23.cycomps import *  # Importing components from cycomps

def saveModel(model, filename: str):
    with open(filename, "wb") as f:
        pickle.dump(model, f)

def loadModel(filename: str):
    with open(filename, "rb") as f:
        return pickle.load(f)

class Table:
    """
    This class represents a lookup table that interpolates data from a CSV file.
    
    Attributes:
    -----------
    argument : tuple, optional
        A tuple containing an object and an attribute that is used for table lookup.
    tabInt : scipy.interpolate.interp1d
        Interpolation function created from the CSV data.
        
    Methods:
    --------
    __call__:
        Returns the interpolated value for the specified argument attribute.
    """
    
    def __init__(self, csvData, interpType='linear'):
        """
        Initializes the Table object by reading a CSV file and creating an interpolation function.

        Parameters:
        -----------
        csvData : str
            The path to the CSV file containing the data.
        interpType : str, optional
            The interpolation method ('linear', 'nearest', 'zero', etc.). Default is 'linear'.
        """
        self.argument = None  # (obj, attr) tuple for interpolation argument
        # Read the CSV data into a pandas DataFrame
        with open(csvData) as data:
            tableData = pd.read_csv(data, delimiter=',', index_col=0)

        tableData = tableData.dropna()  # Drop rows with NaN values
        xs = np.array(tableData.index.astype(float))  # Independent variable (index of the table)
        ys = tableData.values[:, 0]  # Dependent variable (first column of the table)
        
        # Create an interpolation function
        self.tabInt = interp1d(xs, ys, kind=interpType, axis=0, bounds_error=False, fill_value="extrapolate")
    
    def __call__(self):
        """
        Interpolates the value for the given argument attribute.

        Returns:
        --------
        float
            The interpolated value for the argument attribute.

        Raises:
        -------
        UndefinedTableArgument
            If the argument has not been properly defined.
        """
        if self.argument[0] is not None:
            # Interpolate using the argument attribute value
            out = float(self.tabInt(getattr(self.argument[0], self.argument[1])))
            return out
        else:
            raise UndefinedTableArgument("The argument for this table has not been properly defined.")

class DataSet:
    """
    DataSet

    Utility for loading measurement / tuning data from CSV or Excel files.

    Assumptions:
    -----------
    - First row contains column headers.
    - Each column header is the desired variable name.
    - All columns are numeric (or at least convertible to numeric).
    """

    def __init__(
        self,
        filename: str,
        sheet_name: int | str | None = 0,
        index_col: int | str | None = None,
        dropna: bool = True,
    ):
        """
        Parameters
        ----------
        filename : str
            Path to CSV or Excel file.
        sheet_name : int | str | None, optional
            Sheet name or index for Excel files. Ignored for CSV.
        index_col : int | str | None, optional
            Column to use as index (for example a test point ID).
            If None, no explicit index is set.
        dropna : bool, optional
            If True, drop any rows that contain NaNs in any column.
        """
        self.filename = filename

        lower = filename.lower()
        if lower.endswith((".xls", ".xlsx", ".xlsm", ".xlsb", ".ods")):
            df = pd.read_excel(filename, sheet_name=sheet_name)
        else:
            # Default to CSV
            df = pd.read_csv(filename)

        # Optionally set an index (e.g. "Point" or "TestID")
        if index_col is not None:
            df = df.set_index(index_col)

        # Optionally drop incomplete rows
        if dropna:
            df = df.dropna()

        # Force numeric where possible
        df = df.apply(pd.to_numeric, errors="coerce")
        if dropna:
            df = df.dropna()

        self.df = df
        self.columns: list[str] = list(df.columns)
        self.index = df.index
        # Dense 2D NumPy array of shape (n_points, n_vars)
        self.values: np.ndarray = df.to_numpy()

    def __len__(self) -> int:
        """Number of data points (rows)."""
        return self.values.shape[0]

    def __getitem__(self, key: str) -> np.ndarray:
        """
        Column access by key: data['Tt41'].

        Returns the column as a 1D NumPy array.
        """
        return self.df[key].to_numpy()

    def __getattr__(self, name: str) -> np.ndarray:
        """
        Column access by attribute: data.Tt41.

        This is just a convenience on top of __getitem__.
        """
        if name in self.df.columns:
            return self.df[name].to_numpy()
        raise AttributeError(f"'DataSet' object has no attribute '{name}'")

    def row_dict(self, i: int) -> dict:
        """
        Return row i as a dict {column_name: value}.

        Useful if you want to feed a single test point to a model.
        """
        return dict(zip(self.columns, self.values[i, :]))

def write(model, outfile: Optional[str] = None):
    """
    Writes the model data to an Excel file or prints selected data to the console.

    Parameters:
    -----------
    model : object
        The model object containing the data to be written or printed.
    outfile : str, optional
        The output filename (without extension) to write data into an Excel file. 
        If not provided, the data will be printed to the console.
    """
    # Flow station variables to collect data from the model
    ID_FS, station, W, Wc, Tt, Pt, ht, gamma, cp, FAR, Ps, M, V, A, R = ([] for i in range(15))
    # Compressor specific variables
    ID_C, PWR_C, PR, eff_C, Nc_C, beta, ScaleWc_C, ScalePR, ScaleAddEff_C = ([] for i in range(9))
    # Turbine specific variables
    ID_T, ER, PWR_T, eff_T, Nc_T, ScaleWc_T, ScaleER, ScaleAddEff_T  = ([] for i in range(8))
    # Shaft specific variables
    ID_S, Nmech, mechEff, torque, PWR_S = ([] for i in range(5))
    # Burner specific variables
    ID_B, FAR_B, combEff, Wfuel= ([] for i in range(4))
    
    # Extract data from the model components
    for comp in model.comps:
        if hasattr(comp, 'inlet'):  # Flow station data
            ID_FS.append(comp.ID)
            station.append(comp.station)
            W.append(comp.inlet.W)
            Wc.append(comp.inlet.Wc)
            Tt.append(comp.inlet.Tt)
            Pt.append(comp.inlet.Pt)
            ht.append(comp.inlet.ht)
            gamma.append(comp.inlet.gamma)
            cp.append(comp.inlet.cp)
            FAR.append(comp.inlet.FAR)
            Ps.append(comp.inlet.Ps)
            M.append(comp.inlet.M)
            V.append(comp.inlet.V)
            A.append(comp.inlet.A)
            R.append(comp.inlet.R)
        
        # Compressor specific data
        if isinstance(comp, compressor):
            ID_C.append(comp.ID)
            PWR_C.append(comp.power)
            PR.append(comp.PR)
            eff_C.append(comp.eff)
            Nc_C.append(comp.Nc)
            beta.append(comp.beta)
            ScaleWc_C.append(comp.DPscaleWc)
            ScalePR.append(comp.DPscalePR)
            ScaleAddEff_C.append(comp.DPaddEff)

        
        # Turbine specific data
        if isinstance(comp, turbine):
            ID_T.append(comp.ID)
            PWR_T.append(comp.power)
            ER.append(comp.ER)
            eff_T.append(comp.eff)
            Nc_T.append(comp.Nc)
            ScaleWc_T.append(comp.DPscaleWc)
            ScaleER.append(comp.DPscaleER)
            ScaleAddEff_T.append(comp.DPaddEff)

        # Shaft specific data
        if isinstance(comp, shaft):
            ID_S.append(comp.ID)
            Nmech.append(comp.Nmech)
            mechEff.append(comp.mechEff)
            torque.append(comp.torque)
            PWR_S.append(comp.power)

        # Burner specific data
        if isinstance(comp, burner):
            ID_B.append(comp.ID)
            FAR_B.append(comp.FAR)
            combEff.append(comp.eff)
            Wfuel.append(comp.Wf)

    # Create pandas DataFrames for various component data
    flowstations = pd.DataFrame({
        'ID': ID_FS, 'ST': station,'W [kg/s]': W, 'Wc [kg/s]': Wc, 'Tt [K]': Tt, 'Pt [kPa]': [x / 1E3 for x in Pt],
        'ht [kJ/kg]': [x / 1E3 for x in ht], 'gamma': gamma, 'cp [J/(kgK)]': cp, 'R [J/(kgK)]': R,
        'FAR': FAR, 'M': M, 'Ps [kPa]': [x / 1E3 if x is not None else None for x in Ps], 'A [m2]': A, 'V [m/s]': V
    })
    
    compdata = pd.DataFrame({
        'ID': ID_C, 'PWR [kW]': [x / 1E3 for x in PWR_C], 'PR': PR, 'EFF': eff_C, 'Nc [RPM]': Nc_C, 'Beta': beta, 'Map Flow Scalar': ScaleWc_C, 'Map PR Scalar': ScalePR, 'Map Add Eff Scalar': ScaleAddEff_C
    })

    turbdata = pd.DataFrame({
        'ID': ID_T, 'PWR [kW]': [x / 1E3 for x in PWR_T], 'ER': ER, 'EFF': eff_T, 'Nc [RPM]': Nc_T, 'Map Flow Scalar': ScaleWc_T, 'Map ER Scalar': ScaleER, 'Map Add Eff Scalar': ScaleAddEff_T
    })

    shaftdata = pd.DataFrame({
        'ID': ID_S, 'Nmech [RPM]': Nmech, 'Mech.Eff': mechEff, 'Torque [N-m]': torque, 'PWR [kW]': [x / 1E3 for x in PWR_S]
    })
    burnerdata = pd.DataFrame({
        'ID': ID_B, 'FAR': FAR_B, 'Comb.Eff': combEff, 'W.fuel [kg/s]': Wfuel
    })

    paramdata = pd.DataFrame(list(model.params.__dict__.items())[1:], columns=['Name', 'Value'])

    # Write to Excel or print data to console
    if outfile is not None:
        # Function to handle Excel writer setup
        def setWriter() -> pd.ExcelWriter:
            try:
                writer = pd.ExcelWriter(outfile + '.xlsx', engine='openpyxl', mode='a', if_sheet_exists='overlay')
            except PermissionError:
                print(f'WRITE ERROR: Cannot write to {outfile}.xlsx. Please close the file and hit ENTER.')
                input()
                writer = setWriter()
            except FileNotFoundError:
                writer = pd.ExcelWriter(outfile + '.xlsx', engine='openpyxl', mode='w')
            return writer

        writer = setWriter()

        sheet = datetime.now().isoformat().replace(':', '')  # Use current timestamp as sheet name
        # Write DataFrames to the Excel sheet
        flowstations.to_excel(writer, sheet_name=sheet, startrow=1, startcol=0)
        compdata.to_excel(writer, sheet_name=sheet, startrow=1, startcol=17)
        turbdata.to_excel(writer, sheet_name=sheet, startrow=3 + len(ID_C), startcol=17)
        shaftdata.to_excel(writer, sheet_name=sheet, startrow=5 + len(ID_C) + len(ID_T), startcol=17)
        burnerdata.to_excel(writer, sheet_name=sheet, startrow=7 + len(ID_C) + len(ID_T)+ len(ID_S), startcol=17)
        paramdata.to_excel(writer, sheet_name=sheet, startrow=9 + len(ID_C) + len(ID_T) + len(ID_S) + len(ID_B), startcol=17)

        # Set worksheet metadata and close writer
        worksheet = writer.sheets[sheet]
        worksheet.cell(row=1, column=1).value = model.status
        worksheet.cell(row=2, column=18).value = 'COMP.'
        worksheet.cell(row=4 + len(ID_C), column=18).value = 'TURB.'
        worksheet.cell(row=6 + len(ID_C) + len(ID_T), column=18).value = 'SHFT.'
        worksheet.cell(row=8 + len(ID_C) + len(ID_T)+ len(ID_S), column=18).value = 'BRNR.'
        worksheet.cell(row=10 + len(ID_C) + len(ID_T)+ len(ID_S)+ len(ID_B), column=18).value = 'PARAMS.'
        writer.close()
    else:
        # Print flow station data to the console if no file is specified
        print(flowstations[['ID', 'ST', 'W [kg/s]', 'Wc [kg/s]', 'Tt [K]', 'Pt [kPa]']])


class CaseWriter:
    """
    Collect model states across multiple runs and write them in a case-row format.

    Usage pattern
    -------------
    crw = CaseWriter()

    for i, some_setting in enumerate(cases):
        # set up the case
        # ...
        solver()          # or sizeEngine(), tetmatch(), etc.
        crw.update(model, case_label=f"case_{i}")

    crw.write("my_cases")  # produces my_cases.xlsx with one sheet of rows

    The columns are derived from the model structure on the first call to update. Each
    subsequent call assumes the same model structure (same components and IDs).
    """

    def __init__(self, model):
        self._rows = []          # list of dicts, each dict is a case
        self._columns = None     # column ordering fixed from the first row
        self.model=model        # cy23 model instance

    def update(self, case_label: Optional[str] = None):
        """
        Capture the current state of `model` as one row.
        Parameters
        ----------
        case_label : str, optional
            Label for this case. If not given, an integer index is used.
        """
        row = self._extract_row(case_label)

        if self._columns is None:
            # Fix the column order based on the first row
            self._columns = list(row.keys())

        self._rows.append(row)

    def _extract_row(self, case_label: Optional[str]):
        """
        Flatten the model state into a single dictionary suitable for
        case-row output. Column names are formatted for human readability.
        """

        row = {}

        # Case identifier and status
        row["case"] = case_label if case_label is not None else len(self._rows) + 1
        row["status"] = getattr(self.model, "status", None)

        # ------------------------------------------------------------
        # Flow-station quantities (for all components with an inlet)
        # ------------------------------------------------------------
        for comp in self.model.comps:
            if hasattr(comp, "inlet"):
                prefix = f"{comp.ID}.inlet"
                row[f"{prefix} W (kg/s)"] = comp.inlet.W
                row[f"{prefix} Wc (kg/s)"] = comp.inlet.Wc
                row[f"{prefix} Tt (K)"] = comp.inlet.Tt
                row[f"{prefix} Pt (kPa)"] = comp.inlet.Pt / 1e3 if comp.inlet.Pt is not None else None
                row[f"{prefix} ht (kJ/kg)"] = comp.inlet.ht / 1e3 if comp.inlet.ht is not None else None
                row[f"{prefix} gamma"] = comp.inlet.gamma
                row[f"{prefix} cp (J/kg-K)"] = comp.inlet.cp
                row[f"{prefix} R (J/kg-K)"] = comp.inlet.R
                row[f"{prefix} FAR"] = comp.inlet.FAR
                row[f"{prefix} Ps (kPa)"] = comp.inlet.Ps / 1e3 if comp.inlet.Ps is not None else None
                row[f"{prefix} M"] = comp.inlet.M
                row[f"{prefix} V (m/s)"] = comp.inlet.V
                row[f"{prefix} A (m^2)"] = comp.inlet.A

        # ------------------------------------------------------------
        # Compressor quantities
        # ------------------------------------------------------------
        for comp in self.model.comps:
            if isinstance(comp, compressor):
                prefix = f"{comp.ID} (compressor)"

                row[f"{prefix} PWR (kW)"] = comp.power / 1e3 if comp.power is not None else None
                row[f"{prefix} Wc"] = comp.inlet.Wc
                row[f"{prefix} PR"] = comp.PR
                row[f"{prefix} EFF"] = comp.eff
                row[f"{prefix} Nc"] = comp.Nc
                row[f"{prefix} Beta"] = comp.beta

        # ------------------------------------------------------------
        # Turbine quantities
        # ------------------------------------------------------------
        for comp in self.model.comps:
            if isinstance(comp, turbine):
                prefix = f"{comp.ID} (turbine)"

                row[f"{prefix} PWR (kW)"] = comp.power / 1e3 if comp.power is not None else None
                row[f"{prefix} Wc"] = comp.inlet.Wc
                row[f"{prefix} ER"] = comp.ER
                row[f"{prefix} EFF"] = comp.eff
                row[f"{prefix} Nc"] = comp.Nc

        # ------------------------------------------------------------
        # Shaft quantities
        # ------------------------------------------------------------
        for comp in self.model.comps:
            if isinstance(comp, shaft):
                prefix = f"{comp.ID} (shaft)"

                row[f"{prefix} Nmech"] = comp.Nmech
                row[f"{prefix} Mech Eff"] = comp.mechEff
                row[f"{prefix} Torque (N-m)"] = comp.torque
                row[f"{prefix} PWR (kW)"] = comp.power / 1e3 if comp.power is not None else None

        # ------------------------------------------------------------
        # Burner quantities
        # ------------------------------------------------------------
        for comp in self.model.comps:
            if isinstance(comp, burner):
                prefix = f"{comp.ID} (burner)"
                row[f"{prefix} FAR"] = comp.FAR
                row[f"{prefix} Comb Eff"] = comp.eff
                row[f"{prefix} Wfuel (kg/s)"] = comp.Wf

        # ------------------------------------------------------------
        # Nozzle quantities
        # ------------------------------------------------------------
        for comp in self.model.comps:
            if isinstance(comp, nozzle):
                prefix = f"{comp.ID} (nozzle)"
                row[f"{prefix} Fg (kN)"] = comp.Fg/1000
                row[f"{prefix} Cd"] = comp.Cd
                row[f"{prefix} Cv"] = comp.Cv
                row[f"{prefix} Ath (m2)"] = comp.Ath
                row[f"{prefix} Aex (m2)"] = comp.Aex
                row[f"{prefix} Pback (kPa)"] = comp.Pback/1000

        # ------------------------------------------------------------
        # Load quantities
        # ------------------------------------------------------------
        for comp in self.model.comps:
            if isinstance(comp, load):
                prefix = f"{comp.ID} (load)"
                row[f"{prefix} Shaft Power (W)"] = comp.power
                row[f"{prefix} Gen Eff"] = comp.eff
                row[f"{prefix} Elec Power (kWe)"] = comp.kWe


        # ------------------------------------------------------------
        # Global parameters (eng.params)
        # ------------------------------------------------------------
        param_items = list(self.model.params.__dict__.items())[1:]
        for name, value in param_items:
            row[f"param {name}"] = value

        return row

    def to_dataframe(self) -> pd.DataFrame:
        """
        Return the collected cases as a pandas DataFrame.

        Behavior:
        - Inlet A, V, Ps, M columns are dropped if any row has None/NaN.
        - Inlet columns are grouped by variable type (all W, then Wc, Tt, etc.),
        and within each variable group, components follow model instance order.
        - Blank spacing columns separate inlet variable groups, inlet vs. components,
        and components vs. params.
        - For non-converged cases (status != "Converged"), only 'case' and
        'status' are kept; all other columns are left blank.
        """
        if not self._rows:
            return pd.DataFrame()

        df = pd.DataFrame(self._rows)

        # ------------------------------------------------------------
        # 1. Drop inlet A/V/Ps/M columns if any row is missing
        # ------------------------------------------------------------
        conditional_suffixes = [
            " A (m^2)",
            " V (m/s)",
            " Ps (kPa)",
            " M",
        ]

        cols_to_drop = []
        for col in df.columns:
            if any(col.endswith(suffix) for suffix in conditional_suffixes):
                if df[col].isna().any():
                    cols_to_drop.append(col)

        if cols_to_drop:
            df = df.drop(columns=cols_to_drop)

        # ------------------------------------------------------------
        # 2. Group columns: inlet, components, params, with spacing
        # ------------------------------------------------------------
        all_cols = list(df.columns)

        case_cols = [c for c in all_cols if c in ("case", "status")]
        inlet_cols = [c for c in all_cols if ".inlet " in c]
        comp_cols = [
            c for c in all_cols
            if (" (compressor)" in c
                or " (turbine)" in c
                or " (shaft)" in c
                or " (burner)" in c
                or " (nozzle)" in c
                or " (load)" in c)
        ]
        param_cols = [c for c in all_cols if c.startswith("param ")]
        other_cols = [
            c for c in all_cols
            if c not in case_cols + inlet_cols + comp_cols + param_cols
        ]

        # Build inlet columns grouped by variable type,
        # while preserving model instance order via inlet_cols order.
        grouped_inlet_cols: list[str] = []

        if inlet_cols:
            # Define logical variable groups (flows, temps/enthalpy, pressures, gas props, geometry/velocity)
            var_groups = [
                ["W (kg/s)", "Wc (kg/s)"],                 # flows
                ["Tt (K)", "ht (kJ/kg)"],                  # temperatures / enthalpy
                ["Pt (kPa)", "Ps (kPa)"],                  # pressures
                ["gamma", "cp (J/kg-K)", "R (J/kg-K)", "FAR"],  # gas props / mixture
                ["M", "V (m/s)", "A (m^2)"],               # Mach, velocity, area
            ]

            # We'll use "   " (three spaces) as the generic inlet-group separator
            inlet_group_sep = "   "

            # Base ordering: scan inlet_cols (which already reflects model.comps order)
            # for each variable in each group.
            for g_idx, group in enumerate(var_groups):
                group_added_any = False
                for var in group:
                    pattern = f" {var}"
                    for col in inlet_cols:
                        if col.endswith(pattern) and col not in grouped_inlet_cols:
                            grouped_inlet_cols.append(col)
                            group_added_any = True

                # Add a separator after each non-empty group, except we might
                # trim the last one later if it ends the inlet section.
                if group_added_any:
                    grouped_inlet_cols.append(inlet_group_sep)

            # Remove trailing separator if it ended up last
            if grouped_inlet_cols and grouped_inlet_cols[-1] == inlet_group_sep:
                grouped_inlet_cols.pop()

            # Any inlet columns not captured by the groups (unexpected names) go at the end
            for col in inlet_cols:
                if col not in grouped_inlet_cols:
                    grouped_inlet_cols.append(col)

            # Ensure the separator column exists in df as blanks
            if "   " not in df.columns:
                df["   "] = ""

        # Spacing columns between inlet and components, components and params
        section_seps = []
        if grouped_inlet_cols and comp_cols:
            section_seps.append(" ")      # between inlet and components
        if (grouped_inlet_cols or comp_cols) and param_cols:
            section_seps.append("  ")     # between components and params

        for sep in section_seps:
            if sep not in df.columns:
                df[sep] = ""

        # Build final column order
        new_cols: list[str] = []

        # 1) Case identifiers and any "other" columns
        new_cols.extend(case_cols)
        new_cols.extend(other_cols)

        # 2) Inlet data (with internal group separators), then inlet–component separator
        if grouped_inlet_cols:
            new_cols.extend(grouped_inlet_cols)
        if grouped_inlet_cols and comp_cols:
            new_cols.append(" ")

        # 3) Component data, then component–param separator
        if comp_cols:
            new_cols.extend(comp_cols)
        if (grouped_inlet_cols or comp_cols) and param_cols:
            new_cols.append("  ")

        # 4) Parameter data
        if param_cols:
            new_cols.extend(param_cols)

        # Reindex DataFrame with the new column ordering
        df = df[new_cols]

        # ------------------------------------------------------------
        # 3. Blank out data for non-converged cases
        # ------------------------------------------------------------
        if "status" in df.columns:
            status_str = df["status"].astype(str)
            mask_nonconv = status_str != "Converged"

            # Keep 'case' and 'status'; blank everything else
            keep_cols = {"case", "status"}
            # Remove duplicates while preserving column order
            cols_to_clear = list(dict.fromkeys(c for c in df.columns if c not in keep_cols))

            df.loc[mask_nonconv, cols_to_clear] = pd.NA


        return df

    def write(self, outfile: Optional[str] = None, sheet_name: Optional[str] = None):
        """
        Write all collected cases to Excel or print to console.

        Parameters
        ----------
        outfile : str, optional
            Base filename (without extension) for the Excel file.
            If not provided, the DataFrame is printed to the console.
        sheet_name : str, optional
            Name of the Excel sheet. If not given, a timestamp is used.
        """
        df = self.to_dataframe()

        if df.empty:
            print("CaseWriter: no cases collected.")
            return

        if outfile is None:
            print(df)
            return

        if sheet_name is None:
            sheet_name = datetime.now().isoformat().replace(":", "")

        def setWriter() -> pd.ExcelWriter:
            try:
                writer = pd.ExcelWriter(
                    outfile + ".xlsx",
                    engine="openpyxl",
                    mode="a",
                    if_sheet_exists="overlay",
                )
            except PermissionError:
                print(f"WRITE ERROR: Cannot write to {outfile}.xlsx. Please close the file and hit ENTER.")
                input()
                writer = setWriter()
            except FileNotFoundError:
                writer = pd.ExcelWriter(
                    outfile + ".xlsx",
                    engine="openpyxl",
                    mode="w",
                )
            return writer

        writer = setWriter()
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        worksheet = writer.sheets[sheet_name]
        worksheet.cell(row=1, column=1).value = "case-row output"
        writer.close()

    def clear(self):
        """
        Delete data

        This is useful if wanting to dump data into different sheets during 
        single run. Data can be cleared after write command and updated again.
        """

        self._rows = []

# ---------------------------------------------------------------------------
# ParamWriter: case row export for user-selected quantities
# ---------------------------------------------------------------------------
class ParamWriter:
    """
    Lightweight case-row writer that records only selected quantities from
    a model over multiple runs.

    Each quantity is specified as either:

        - A simple name: "FN_kN"
          -> first tries model.params.FN_kN, then model.FN_kN

        - A dotted path from the model root: "params.TSFC_US",
          "HPC.PR", "noz.inlet.M", etc.
          -> interpreted literally as getattr(getattr(...), ...)


    Typical usage
    -------------
    pcw = ParamWriter(model,[
        "FN_kN",                 # model.params.FN_kN
        "TSFC_US",               # model.params.TSFC_US
        "Eff_oall",              # model.params.Eff_oall
        "HPC.PR",                # model.HPC.PR
    ])

    for i, some_setting in enumerate(cases):
        # set case conditions and run model
        # ...
        sizeEngine()
        tetmatch()
        performance()  # update eng.params.*

        pcw.update(eng, case_label=f"case_{i}")

    pcw.write("exampleParams")  # writes exampleParams.xlsx with one row per case

    Behavior
    --------
    - Columns are: "case", "status", then one column per spec in the order
      given at construction (using the column name or alias provided).
    - For non-converged cases (model.status != "Converged"), only "case"
      and "status" are kept; all parameter columns are left blank so the
      row is visible but has no numeric data.
    """

    def __init__(self, model, specs: list):

        self.model=model
        self._rows: list[dict] = []
        self._fields=specs

    def _get_by_path(self, path: str):
        """
        Resolve a dotted path.

        Examples
        --------
        path = "model.params.FN_kN" -> model.params.FN_kN
        path = "model.HPC.PR"       -> model.HPC.PR
        path = "model.noz.inlet.M"  -> model.noz.inlet.M
        """
        obj = self.model
        for name in path.split("."):
            if not hasattr(obj, name):
                return None
            obj = getattr(obj, name)
        return obj


    def update(self, case_label: Optional[str] = None):
        """
        Capture the current state of `model` as a single row.

        Parameters
        ----------
        model : object
            Engine model object (cy23 graph instance).
        case_label : str, optional
            Label for this case. If not given, an integer index is used.
        """
        row: dict = {}

        # Case identifier and status
        row["case"] = case_label if case_label is not None else len(self._rows) + 1

        # Selected quantities
        for path in self._fields:
            row[path] = self._get_by_path(path)

        self._rows.append(row)

    def write(self, outfile: Optional[str] = None, sheet_name: Optional[str] = None):
        """
        Write all collected cases to Excel or print to console.

        Parameters
        ----------
        outfile : str, optional
            Base filename (without extension) for the Excel file.
            If not provided, the DataFrame is printed to the console.
        sheet_name : str, optional
            Name of the Excel sheet. If not given, a timestamp is used.
        """
        df = pd.DataFrame(self._rows)

        if df.empty:
            print("ParamCaseWriter: no cases collected.")
            return

        if outfile is None:
            print(df)
            return

        if sheet_name is None:
            sheet_name = datetime.now().isoformat().replace(":", "")

        def setWriter() -> pd.ExcelWriter:
            try:
                writer = pd.ExcelWriter(
                    outfile + ".xlsx",
                    engine="openpyxl",
                    mode="a",
                    if_sheet_exists="overlay",
                )
            except PermissionError:
                print(f"WRITE ERROR: Cannot write to {outfile}.xlsx. Please close the file and hit ENTER.")
                input()
                writer = setWriter()
            except FileNotFoundError:
                writer = pd.ExcelWriter(
                    outfile + ".xlsx",
                    engine="openpyxl",
                    mode="w",
                )
            return writer

        writer = setWriter()
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        worksheet = writer.sheets[sheet_name]
        worksheet.cell(row=1, column=1).value = "Case"
        writer.close()
    
    def clear(self):
        """
        Delete data

        This is useful if wanting to dump data into different sheets during 
        single run. Data can be cleared after write command and updated again.
        """
        
        self._rows = []

class BatchRun:
    def __init__(self, solver, inputs, outputs):
        self.model=solver.model
        self.solver=solver # Solver to be called to run with these inputs to get outputs.
        self.inputs = inputs # List of input strings.
        self.outputs= inputs + outputs # List of output strings. Inputs are appended.
        self._prev_values = {name: None for name in self.inputs}
        self.max_blank_rows=2

    def _get_by_path(self, path: str):
        """
        Resolve a dotted path. 
        """
        obj = self.model
        for name in path.split("."):
            if hasattr(obj, name):
                obj = getattr(obj, name)
            else:
                obj="UNDEFINED"
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

    def _row_is_empty(self, ws, row_index):
        input_start_col = 2
        return all(
            ws.cell(row=row_index, column=col).value in (None, "")
            for col in range(input_start_col, len(self.inputs) + input_start_col)
        )

    def __call__(self, workbook_name, sheet_name=None):
        self._prev_values = {name: self._get_by_path(name) for name in self.inputs}
        workbook_path=Path(workbook_name)
        wb = load_workbook(workbook_path)
        # Read from and write back to the same worksheet.
        if sheet_name is not None:
            ws = wb[sheet_name]
        else:
            ws=wb.active

        input_count = len(self.inputs)
        input_start_col = 2      # Column B
        spacer_col = input_start_col + input_count
        timestamp_col=spacer_col+1
        status_col= timestamp_col + 1
        output_start_col = status_col + 1  # one blank spacer column
        
        # Write headers for output columns
        ws.cell(row=1, column=spacer_col, value="")
        ws.cell(row=1, column=timestamp_col, value="Timestamp")
        ws.cell(row=1, column=status_col, value="Status")
        for j, out_name in enumerate(self.outputs, start=output_start_col):
            ws.cell(row=1, column=j, value=out_name)

        # Clear Previous results
        for row in range(2, ws.max_row + 1):
            ws.cell(row=row, column=timestamp_col).value = None
            ws.cell(row=row, column=status_col).value = None
            for col in range(output_start_col, output_start_col + len(self.outputs)):
                ws.cell(row=row, column=col).value = None

        # Process each data row
        blank_row_count = 0
        for row_index in range(2, ws.max_row + 1):
            if self._row_is_empty(ws, row_index):
                blank_row_count += 1
                if blank_row_count >= self.max_blank_rows:
                    break
                continue
            else:
                blank_row_count = 0

            for col_index, var_name in enumerate(self.inputs, start=input_start_col):
                value = ws.cell(row=row_index, column=col_index).value

                # Carry forward blank cells from previous row
                if value in (None, ""):
                    value = self._prev_values[var_name]
                else:
                    self._prev_values[var_name] = value

                self._set_by_path(var_name, value)
            
            self.solver()

            # Write outputs into the same row, after one spacer column
            # Include status and timestamp columns
            ws.cell(row=row_index, column=timestamp_col, value=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            ws.cell(row=row_index, column=status_col, value=getattr(self.model, "status", None))
            for j, out_name in enumerate(self.outputs, start=output_start_col):
                ws.cell(row=row_index, column=j, value=self._get_by_path(out_name))

        wb.save(workbook_path)
        wb.close()

def ewma():
    "smoothing function"


class Performance:
    """
    Class to calculate performance parameters for a cy23 engine model.
    """

    def __init__(self, model, arch="jet", params_attr="params"):
        self.model = model
        self.params_attr = params_attr
        self.arch=arch

    def __call__(self):
        if self.arch=="jet":
            self.jetPerf()
        else:
            raise NotImplementedError(f"Performance calculation for architecture '{self.arch}' is not implemented.")

    def jetPerf(self):
        eng = self.model
        amb = None
        splitter_component = None
        thrust_nozzles = []
        fuel_components = []

        for comp in eng.comps:
            if isinstance(comp, ambient):
                amb = comp
            if isinstance(comp, splitter):
                splitter_component = comp
            if isinstance(comp, nozzle):
                thrust_nozzles.append(comp)
            if isinstance(comp, burner):
                fuel_components.append(comp)
        
        params = getattr(eng, self.params_attr)

        N_TO_LBF = 0.224809
        KG_PER_S_TO_LBM_PER_HR = 7936.64
        KG_TO_LBM=2.20462

        if amb is not None:

            # Net thrust from all specified nozzles
            FN_nozzles = 0.0
            EGT = 0.0
            pwr_nozzles = 0.0
            V9max = 0.0
            for noz in thrust_nozzles:
                EGT = max(EGT, noz.inlet.Tt)
                V9max = max(V9max, noz.outlet.V)
                FN_nozzles += noz.Fg
                pwr_nozzles += noz.power

            FN = FN_nozzles - amb.W * amb.V
            lbf = FN * N_TO_LBF
            kN = FN/1000.0
            params.FN_kN = kN
            params.FN_lbf = lbf
            params.EGT_degC = EGT-273.15
            params.V9max = V9max
            params.noz_kW = pwr_nozzles/1000.0


            # Total fuel flow & fuel power
            fuel_power = 0.0
            Wf_total = 0.0
            for comp in fuel_components:
                Wf_total += comp.Wf
                fuel_power += comp.Wf * comp.LHV

            # Avoid division by zero in degenerate cases
            if lbf != 0.0:
                params.TSFC_SI = Wf_total * 1000/ kN
                params.TSFC_US = Wf_total * KG_PER_S_TO_LBM_PER_HR / lbf
            else:
                params.TSFC_SI = float("inf")
                params.TSFC_US = float("inf")

            # Specific thrust
            params.SFn_SI=FN/amb.W
            params.SFn_US=lbf/(amb.W*KG_TO_LBM)

            # Overall efficiency
            TP = FN * amb.V

            if fuel_power != 0.0:
                params.Eff_oall = TP / fuel_power
            else:
                params.Eff_oall = 0.0

            # BPR if splitter exist
            if splitter_component is not None:
                params.BPR=splitter_component.BPR

class UndefinedTableArgument(Exception):
    """Custom exception for undefined table argument errors."""
    pass
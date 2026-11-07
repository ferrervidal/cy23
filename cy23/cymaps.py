"""
CY23.CYMAPS

This module provides functions and classes for loading, processing, and
interpolating compressor and turbine performance maps used by cycle models.

Key Components & Features:
---------------
1. **preProcess**: loads csv compressor map data and resamples at standardized number of 
beta lines (20 by default).
2. **interpMap**: Functions to interpolate processed maps at additional
   corrected speeds and convert them into rectangular spline data. This may be called by 
   the user to interpolate maps at additional speeds from the provided data.
3. **processCmap & processTmap**: Creates interpolators to interpolate maps as desired.
4. **Map Objects**: Callable compressor and turbine map classes that provide
   flow, pressure-ratio, efficiency, and swirl values to the cycle model.
5. **Map Serialization**: Functions to save processed maps as pickle files and
   load them for later use.

Usage:
------
This module is used by the cycle-modeling package to convert CSV performance
maps into callable compressor and turbine map objects. 
"""


from __future__ import annotations

from pathlib import Path
import pickle

import numpy as np
import pandas as pd
import scipy as sp
from scipy.interpolate import PchipInterpolator


_MAP_COLUMNS = 5


class BadMap(Exception):
    """Raised when map data cannot be processed safely."""


def _read_numeric_map(mapfile):
    """Read a CSV map and remove rows containing missing values."""
    try:
        data = pd.read_csv(mapfile, delimiter=",", index_col=None).to_numpy(
            dtype=float
        )
    except (OSError, ValueError, TypeError) as exc:
        raise BadMap(f"Could not read numeric map data from {mapfile!r}.") from exc

    data = data[~np.isnan(data).any(axis=1)]
    if data.size == 0:
        raise BadMap(f"Map {mapfile!r} contains no complete numeric rows.")
    return data


def _require_strictly_increasing(x, description):
    """Validate an interpolation coordinate and return it as a 1-D array."""
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size < 2 or not np.all(np.diff(x) > 0.0):
        raise BadMap(f"{description} must contain at least two strictly increasing values.")
    return x


def _write_pickle(mapfile, value):
    """Write the processed map beside the source CSV using a safe suffix change."""
    output = Path(mapfile).with_suffix(".pkl")
    with output.open("wb") as mapfile_handle:
        pickle.dump(value, mapfile_handle)


def flowFix(Wc):
    """Enforce non-decreasing corrected flow along a speedline."""
    values = np.asarray(Wc)
    if values.ndim != 1:
        values = values.reshape(-1)
    if values.size == 0:
        return values.copy()
    return np.maximum.accumulate(values)


def preProcess(mymap, betadisc=20.0):
    """Read a compressor map and resample each speedline by beta.

    The input CSV is expected to contain corrected speed, corrected flow,
    pressure ratio, and efficiency in its first four columns.  Speedlines are
    expected to be ordered from low to high flow.
    """
    data = _read_numeric_map(mymap)
    if data.shape[1] < 4:
        raise BadMap("A compressor map must have at least four columns.")

    NrT, Wc, PR, Eta = data[:, :4].T
    speeds = np.unique(NrT)

    try:
        beta_count = int(round(float(betadisc)))
    except (TypeError, ValueError) as exc:
        raise BadMap("betadisc must be a positive number.") from exc
    if beta_count < 2:
        raise BadMap("betadisc must produce at least two beta points.")

    new_betas = np.arange(1.0, beta_count + 1.0)
    beta_column = new_betas[:, None]
    characteristics = []

    for speed in speeds:
        mask = NrT == speed
        this_pr = np.asarray(PR[mask], dtype=float)
        this_eta = np.asarray(Eta[mask], dtype=float)
        this_wc = np.asarray(flowFix(Wc[mask]), dtype=float)

        # Extend each speedline to PR=1 so all lines share the same low-PR end.
        if np.min(this_pr) > 1.0:
            this_pr = np.append(this_pr, 1.0)
            this_eta = np.append(this_eta, 0.0)
            this_wc = np.append(this_wc, np.max(this_wc))

        min_flow = np.min(this_wc)
        start_pr = this_pr[np.isclose(this_wc, min_flow)]
        if start_pr.size != 1:
            raise BadMap(
                "Map has multiple values with the same minimum mass flow on one speedline."
            )

        # Parameterize the speedline by normalized distance in the flow/PR plane.
        pr_range = np.max(np.abs(this_pr - start_pr[0]))
        wc_range = np.max(np.abs(this_wc - min_flow))
        if pr_range == 0.0 and wc_range == 0.0:
            raise BadMap(f"Speedline {speed!r} contains no usable variation.")

        dist_pr = (this_pr - start_pr[0]) / (pr_range or 1.0)
        dist_wc = (this_wc - min_flow) / (wc_range or 1.0)
        distance = np.hypot(dist_wc, dist_pr)
        max_distance = np.max(distance)
        if max_distance == 0.0:
            raise BadMap(f"Speedline {speed!r} contains no usable variation.")

        parameter = distance / max_distance * new_betas.size
        order = np.argsort(parameter, kind="stable")
        parameter = parameter[order]
        this_wc = this_wc[order]
        this_pr = this_pr[order]
        this_eta = this_eta[order]

        # PCHIP requires unique x coordinates. Repeated points do not add
        # information and otherwise cause a cryptic SciPy exception.
        unique_parameter, unique_indices = np.unique(parameter, return_index=True)
        if unique_parameter.size < 2:
            raise BadMap(f"Speedline {speed!r} has insufficient unique points.")
        this_wc = this_wc[unique_indices]
        this_pr = this_pr[unique_indices]
        this_eta = this_eta[unique_indices]

        wc_fit = PchipInterpolator(unique_parameter, this_wc, extrapolate=True)
        pr_fit = PchipInterpolator(unique_parameter, this_pr, extrapolate=True)
        eta_fit = PchipInterpolator(unique_parameter, this_eta, extrapolate=True)

        speed_column = np.full((new_betas.size, 1), speed, dtype=float)
        characteristics.append(
            np.column_stack(
                (speed_column[:, 0], wc_fit(new_betas), pr_fit(new_betas),
                 eta_fit(new_betas), beta_column[:, 0])
            )
        )

    return np.vstack(characteristics)


def interpMap(mymap, intspeeds):
    """Interpolate a processed compressor map at additional speed values."""
    map_array = np.asarray(mymap, dtype=float)
    if map_array.ndim != 2 or map_array.shape[1] != _MAP_COLUMNS:
        raise BadMap("A processed compressor map must be a two-dimensional five-column array.")

    requested_speeds = np.asarray(intspeeds, dtype=float).reshape(-1)
    if requested_speeds.size == 0:
        return map_array.copy()

    speeds = _require_strictly_increasing(np.unique(map_array[:, 0]), "Map speeds")
    betas = np.unique(map_array[:, -1])
    interpolated = []

    for beta in betas:
        mask = map_array[:, -1] == beta
        this_map = map_array[mask]
        if this_map.shape[0] != speeds.size:
            raise BadMap("Every beta must have one point at every map speed.")

        # Sort by speed before fitting, independent of the input row order.
        this_map = this_map[np.argsort(this_map[:, 0], kind="stable")]
        wc_fit = PchipInterpolator(speeds, this_map[:, 1], extrapolate=True)
        pr_fit = PchipInterpolator(speeds, this_map[:, 2], extrapolate=True)
        eta_fit = PchipInterpolator(speeds, this_map[:, 3], extrapolate=True)
        interpolated.append(
            np.column_stack(
                (requested_speeds, wc_fit(requested_speeds),
                 pr_fit(requested_speeds), eta_fit(requested_speeds),
                 np.full(requested_speeds.size, beta))
            )
        )

    # Avoid duplicate rows when a requested speed already exists in the map.
    additions = np.vstack(interpolated)
    combined = np.vstack((map_array, additions))
    _, unique_indices = np.unique(combined[:, [0, 4]], axis=0, return_index=True)
    combined = combined[np.sort(unique_indices)]
    return combined[np.argsort(combined[:, 0], kind="stable")]


def outformat(newChic):
    """Convert a processed map to rectangular arrays for SciPy splines."""
    map_array = np.asarray(newChic, dtype=float)
    if map_array.ndim != 2 or map_array.shape[1] != _MAP_COLUMNS:
        raise BadMap("A processed map must be a two-dimensional five-column array.")

    speeds = np.unique(map_array[:, 0])
    betas = np.unique(map_array[:, -1])
    expected_size = speeds.size * betas.size
    if map_array.shape[0] != expected_size:
        raise BadMap("Processed map does not form a complete rectangular speed/beta grid.")

    # Explicit indexing makes the result independent of the input row order.
    speed_index = {value: index for index, value in enumerate(speeds)}
    beta_index = {value: index for index, value in enumerate(betas)}
    rows = np.array([speed_index[value] for value in map_array[:, 0]])
    cols = np.array([beta_index[value] for value in map_array[:, -1]])
    if np.unique(np.column_stack((rows, cols)), axis=0).shape[0] != expected_size:
        raise BadMap("Processed map contains duplicate speed/beta pairs.")

    flow = np.empty((speeds.size, betas.size))
    pressure_ratio = np.empty_like(flow)
    efficiency = np.empty_like(flow)
    flow[rows, cols] = map_array[:, 1]
    pressure_ratio[rows, cols] = map_array[:, 2]
    efficiency[rows, cols] = map_array[:, 3]
    return speeds, betas, flow, pressure_ratio, efficiency


def processCmap(mapfile, interpspeeds=None, betadisc=20.0):
    """Process a compressor map and save the resulting callable map object."""
    proc = preProcess(mapfile, betadisc)
    if interpspeeds is not None:
        proc = interpMap(proc, interpspeeds)

    speeds, betas, flow, pressure_ratio, efficiency = outformat(proc)
    flow_chic = sp.interpolate.RectBivariateSpline(speeds, betas, flow, kx=1, ky=1)
    pr_chic = sp.interpolate.RectBivariateSpline(speeds, betas, pressure_ratio, kx=1, ky=1)
    eff_chic = sp.interpolate.RectBivariateSpline(speeds, betas, efficiency, kx=1, ky=1)
    _write_pickle(mapfile, compressorMap(flow_chic, pr_chic, eff_chic))


def processTmap(mapfile):
    """Process a turbine map in column format and save the map object."""
    data = _read_numeric_map(mapfile)
    if data.shape[1] < 4:
        raise BadMap("A turbine map must have at least four columns.")

    NrT, ER, Wc, Eta = data[:, :4].T
    Swirl = data[:, 4] if data.shape[1] > 4 else np.zeros_like(Eta)
    speeds = np.unique(NrT)

    er_interp = np.linspace(0.0, np.max(ER), 30)
    characteristics = []
    for speed in speeds:
        mask = NrT == speed
        order = np.argsort(ER[mask], kind="stable")
        this_er = _require_strictly_increasing(
            ER[mask][order], f"Expansion ratios at speed {speed!r}"
        )
        this_wc = Wc[mask][order]
        this_eta = Eta[mask][order]
        this_swirl = Swirl[mask][order]

        wc_fit = PchipInterpolator(this_er, this_wc, extrapolate=True)
        eta_fit = PchipInterpolator(this_er, this_eta, extrapolate=True)
        swirl_fit = PchipInterpolator(this_er, this_swirl, extrapolate=True)
        characteristics.append(
            np.column_stack(
                (np.full(er_interp.size, speed), er_interp,
                 wc_fit(er_interp), eta_fit(er_interp), swirl_fit(er_interp))
            )
        )

    processed = np.vstack(characteristics)
    speeds, expansion_ratios, flow, efficiency, swirl = _turbine_outformat(processed)
    flow_chic = sp.interpolate.RectBivariateSpline(speeds, expansion_ratios, flow, kx=1, ky=1)
    eff_chic = sp.interpolate.RectBivariateSpline(speeds, expansion_ratios, efficiency, kx=1, ky=1)
    swirl_chic = sp.interpolate.RectBivariateSpline(speeds, expansion_ratios, swirl, kx=1, ky=1)
    _write_pickle(mapfile, turbineMap(flow_chic, eff_chic, swirl_chic))


def _turbine_outformat(processed):
    """Convert the internal rectangular turbine representation to arrays."""
    speeds = np.unique(processed[:, 0])
    expansion_ratios = np.unique(processed[:, 1])
    expected_size = speeds.size * expansion_ratios.size
    if processed.shape[0] != expected_size:
        raise BadMap("Processed turbine map does not form a rectangular grid.")
    rows = np.repeat(np.arange(speeds.size), expansion_ratios.size)
    cols = np.tile(np.arange(expansion_ratios.size), speeds.size)
    return (
        speeds,
        expansion_ratios,
        processed[:, 2].reshape(speeds.size, expansion_ratios.size),
        processed[:, 3].reshape(speeds.size, expansion_ratios.size),
        processed[:, 4].reshape(speeds.size, expansion_ratios.size),
    )


class compressorMap:
    def __init__(self, flow, PR, eff):
        self.flow = flow
        self.PR = PR
        self.eff = eff

    def __call__(self, Nmap, beta):
        beta = float(beta)
        evaluation_beta = max(1.0, beta)
        Wc = float(self.flow(Nmap, evaluation_beta))
        PR = float(self.PR(Nmap, evaluation_beta))
        eff = float(self.eff(Nmap, evaluation_beta))

        # Preserve the intended linear extension below beta=1.
        if beta < 1.0:
            zeta = abs(beta) if abs(beta) < 1.0 else 1.0 / abs(beta)
            Wc *= zeta
            PR = 1.0 + (PR - 1.0) * zeta

        return Wc, PR, eff


class turbineMap:
    def __init__(self, flow, eff, swirl):
        self.flow = flow
        self.eff = eff
        self.swirl = swirl

    def __call__(self, Nmap, ER):
        return (
            float(self.flow(Nmap, ER)),
            float(self.eff(Nmap, ER)),
            float(self.swirl(Nmap, ER)),
        )


def loadmap(pkl):
    """Load a pickled map object, accepting a path with or without ``.pkl``."""
    path = Path(pkl)
    if path.suffix.lower() != ".pkl":
        path = path.with_name(path.name + ".pkl")
    with path.open("rb") as mapfile:
        return pickle.load(mapfile)

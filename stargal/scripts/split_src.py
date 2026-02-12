"""
Here lies the different classes and functions that are used to run SPLIT.
"""

#--- Import packages ---
import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from . import statsplot as lsp

#--- sed grid class ---
# The class used to organize a sed grid thanks to their physical parameters and their Ids.
class SED_GRID:
    def __init__(self, list_path):
        """Load SED list."""
        sed_list = []
        with open(list_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                filepath = line.split()[0]
                sp_type = filepath.split('/')[-1].replace('.sed', "")
                sed_list.append(sp_type)
        self.sed_list = sed_list
        self.sed_grid = None

    def build(self):
        """Build the sed grid [Id, Teff, logg, FeH] from filenames."""
        pattern = r"Teff([-\d\.]+)_logg([-\d\.]+)_FeH([-\d\.]+)"
        sed_grid = []
        for Id, s in enumerate(self.sed_list):
            match = re.search(pattern, s)
            if match:
                Teff, logg, FeH = match.groups()
                Teff, logg, FeH = float(Teff), float(logg), float(FeH)
                sed_grid.append([Id, Teff, logg, FeH])
        self.sed_grid = np.array(sed_grid)
        return self.sed_grid

    def plot(self, cmap='gnuplot'):
        """Plot 3D SED grid with color depending on FeH."""
        if self.sed_grid is None:
            raise ValueError("Call make_sed_grid before plot.")

        Teff, logg, FeH = self.sed_grid[:, 1], self.sed_grid[:, 2], self.sed_grid[:, 3]

        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(projection='3d')

        # Couleur selon la métallicité
        sc = ax.scatter(
            Teff, logg, FeH,
            c=FeH, cmap=cmap,
            alpha=0.8, s=10, edgecolor="none"
        )

        # Colorbar
        cbar = plt.colorbar(sc, ax=ax, pad=0.1, shrink=0.8)
        cbar.set_label('[Fe/H]', rotation=270, labelpad=15)

        ax.set_xlabel('Teff [K]')
        ax.set_ylabel('log(g)')
        ax.set_zlabel('[Fe/H]')
        ax.set_title("SED Grid Colored by Metallicity")

        # Optionnel : meilleure orientation par défaut
        ax.view_init(elev=20, azim=-45)

        plt.tight_layout()
        plt.show()


#--- star pdf class ---
# The class used to project the fitting parameter (Chi2) of a source to each sed of the 
# precedently computed sed grid thanks to the sed Ids. Then the chi2_grid can be marginalised
# and turn to probability values. 


class STAR_PDF:
    def __init__(self, sed_grid):
        if sed_grid is None:
            raise ValueError("SED_GRID: call sed_grid_obj.build() before using STAR_PDF.")

        self.sed_grid = sed_grid   # base grid
        self.len_sed_grid = len(self.sed_grid)
        self.chi2_grid = None

    def build_chi2_grid(self, chi2_array):
        """Attach chi2 array to SED grid."""
        if len(chi2_array) != self.len_sed_grid:
            raise ValueError("chi2 array length does not match sed_grid length.")

        self.chi2_grid = np.column_stack((self.sed_grid, chi2_array))
        return self.chi2_grid

    def project_min_chi2(self, xaxis="Teff", yaxis=None, to_prob=False, fixed=None):
        """
        Project Chi² parallel to one or two axes.
        Uses self.chi2_grid automatically.
        """
        if self.chi2_grid is None:
            raise ValueError("Call build_chi2_grid before projecting.")

        grid = self.chi2_grid
        columns = ["Id", "Teff", "logg", "FeH", "Chi2"]
        col_idx = {k: i for i, k in enumerate(columns)}

        # 1) Apply fixed filters
        mask = np.ones(len(grid), dtype=bool)
        if fixed:
            for key, val in fixed.items():
                mask &= np.isclose(grid[:, col_idx[key]], val)
            grid = grid[mask]
            if len(grid) == 0:
                raise ValueError("No grid point matches the given fixed conditions.")

        xcol = col_idx[xaxis]

        if yaxis:   # ----------- 2D PROJECTION -------------
            ycol = col_idx[yaxis]

            X_vals = np.unique(grid[:, xcol])
            Y_vals = np.unique(grid[:, ycol])
            Z = np.full((len(Y_vals), len(X_vals)), np.nan)

            for i, xv in enumerate(X_vals):
                for j, yv in enumerate(Y_vals):
                    mask_xy = np.isclose(grid[:, xcol], xv) & np.isclose(grid[:, ycol], yv)
                    if np.any(mask_xy):
                        Z[j, i] = np.nanmin(grid[mask_xy, -1])

            if to_prob:
                Zmin = np.nanmin(Z)
                Z = np.exp(-0.5 * (Z - Zmin))

            return X_vals, Y_vals, Z

        else:       # ----------- 1D PROJECTION -------------
            X_vals = np.unique(grid[:, xcol])
            Z = np.full(len(X_vals), np.nan)

            for i, xv in enumerate(X_vals):
                mask_x = np.isclose(grid[:, xcol], xv)
                if np.any(mask_x):
                    Z[i] = np.nanmin(grid[mask_x, -1])

            if to_prob:
                Zmin = np.nanmin(Z)
                Z = np.exp(-0.5 * (Z - Zmin))

            return X_vals, Z

    def plot_min_chi2(self, xaxis="Teff", yaxis=None, to_prob=False,
                      full=False, cmap='gnuplot', fixed=None):

        if self.chi2_grid is None:
            raise ValueError("Call build_chi2_grid before plotting.")

        grid = self.chi2_grid  # just to avoid rewriting this everywhere

        # ---------- FULL CORNER-PLOT ----------
        if full:
            axes = ["Teff", "logg", "FeH"]
            n = len(axes)

            # Uniform axes
            limits = {}
            for ax in axes:
                limits[ax] = np.unique(grid[:, {"Teff": 1, "logg": 2, "FeH": 3}[ax]])
            
            fig = plt.figure(figsize=(8, 8))
            gs = plt.GridSpec(n, n, wspace=0.0, hspace=0.0)

            for i in range(n):  # line
                for j in range(n):  # col
                    if j > i:
                        continue
                    ax = fig.add_subplot(gs[i, j])

                    if i == j:
                        # --- PDF 1D ---
                        X, Z = self.project_min_chi2(xaxis=axes[j], to_prob=to_prob, fixed=fixed)
                        ax.plot(X, Z, color='black', lw=1)
                        # ax.fill_between(X, Z, color='C0', alpha=0.4)
                        ax.set_xlim(limits[axes[j]].min(), limits[axes[j]].max())
                        ax.set_ylim(Z.min() * 0.95, Z.max() * 1.05)
                        if i < n - 1:
                            ax.set_xticklabels([])
                        if i != 0:
                            ax.yaxis.tick_right()
                        else:
                            ax.set_xlabel(axes[j])

                    else:
                        # --- PDF 2D ---
                        X, Y, Z = self.project_min_chi2(xaxis=axes[j], yaxis=axes[i], fixed=fixed)
                        if to_prob:
                            Z = np.exp(-0.5 * (Z - np.nanmin(Z)))
                        im = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap)
                        ax.set_xlim(limits[axes[j]].min(), limits[axes[j]].max())
                        ax.set_ylim(limits[axes[i]].min(), limits[axes[i]].max())
                        if i < n - 1:
                            ax.set_xticklabels([])
                        else:
                            ax.set_xlabel(axes[j])
                        if j > 0:
                            ax.set_yticklabels([])
                        else:
                            ax.set_ylabel(axes[i])

            # Colorbar
            cbar_ax = fig.add_axes([1, 0.15, 0.02, 0.7])
            cbar = fig.colorbar(im, cax=cbar_ax)
            cbar.set_label("Probability" if to_prob else "min(χ²)")
            plt.subplots_adjust(left=0.08, right=0.9, top=0.9, bottom=0.08)
            plt.show()
            return


        # --- Single plot mode
        if yaxis:
            X, Y, Z = self.project_min_chi2(xaxis, yaxis, to_prob=to_prob, fixed=fixed)
            plt.figure(figsize=(6, 4))
            plt.pcolormesh(X, Y, Z, shading='auto')
            plt.colorbar(label="Probability" if to_prob else "min(Chi²)")
            plt.xlabel(xaxis)
            plt.ylabel(yaxis)
            plt.title(f"Minimized Chi² map")
            plt.show()
        else:
            X, Z = self.project_min_chi2(xaxis, to_prob=to_prob, fixed=fixed)
            plt.figure(figsize=(6, 3))
            plt.plot(X, Z, marker='o', lw=1)
            plt.xlabel(xaxis)
            plt.ylabel("Probability" if to_prob else "min(Chi²)")
            plt.title("Minimized Chi² PDF")
            plt.grid(True)
            plt.show()

    def plot_full_chi2(self, cmap='gnuplot'):
        """Plot 3D SED grid with color depending on FeH."""
        grid = self.chi2_grid  # just to avoid rewriting this everywhere
        # print(grid)
        Teff, logg, FeH = grid[:, 1], grid[:, 2], grid[:, 3]
        Chi2 = grid[:, 4]
        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(projection='3d')
        # Coulor as a function of Chi2
        sc = ax.scatter(
            Teff, logg, FeH,
            c=Chi2, cmap=cmap,
            alpha=0.8, s=10, edgecolor="none"
        )

        # Colorbar
        cbar = plt.colorbar(sc, ax=ax, pad=0.1, shrink=0.8)
        cbar.set_label('Chi2', rotation=270, labelpad=15)

        ax.set_xlabel('Teff [K]')
        ax.set_ylabel('log(g)')
        ax.set_zlabel('[Fe/H]')
        ax.set_title("SED Grid Colored by Metallicity")
 
        # Optionnel : meilleure orientation par défaut
        ax.view_init(elev=20, azim=-45)
 
        plt.tight_layout()
        plt.show()

    def best_output(self):
        grid = self.chi2_grid
        best_sed_id = np.argmin(grid[:,4])
        chi2_best = min(grid[:,4])
        Teff_best = grid[best_sed_id, 1]
        logg_best = grid[best_sed_id, 2]
        FeH_best = grid[best_sed_id, 3]
        
        chisec_grid = np.delete(grid, best_sed_id, axis=0)
        sec_best_sed_id = np.argmin(chisec_grid[:,4])
        sec_chi2_best = min(chisec_grid[:,4])
        sec_Teff_best = chisec_grid[sec_best_sed_id, 1]
        sec_logg_best = chisec_grid[sec_best_sed_id, 2]
        sec_FeH_best = chisec_grid[sec_best_sed_id, 3]
        sec_best_sed_id = sec_best_sed_id + 1 if sec_best_sed_id >= best_sed_id else sec_best_sed_id
        return [best_sed_id, Teff_best, logg_best, FeH_best, chi2_best, 
                sec_best_sed_id, sec_Teff_best, sec_logg_best, sec_FeH_best, sec_chi2_best]


#--- star pdf analyzer class---
# The class use to compute the different statistical metrics of a star pdf
# in 1, 2 or 3D. The main tool to determine the source type from star pdfs.

from scipy.signal import find_peaks
from scipy.interpolate import interp1d, griddata, RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from numba import njit

class STAR_PDF_ANALYZER:
    """
    Analyzer for chi2_grid (Nx5 array: [Id, Teff, logg, FeH, Chi2]).
    Option A: interpolation is done on chi2 (not on probabilities).
    - to_prob_3D: if True, compute global P3 (full-grid PDF) once at init
    - to_prob: if True, projectors return normalized probabilities per projection
    """
    def __init__(self, chi2_grid, to_prob=True, to_prob_3D=False,
                oversample=1.0, nmax=800):

        chi2_grid = np.asarray(chi2_grid)
        if chi2_grid.ndim != 2 or chi2_grid.shape[1] < 5:
            raise ValueError("chi2_grid must be shape (N,>=5)")

        self.grid = chi2_grid
        self.to_prob = bool(to_prob)
        self.to_prob_3D = bool(to_prob_3D)

        self.col = {"Id": 0, "Teff": 1, "logg": 2, "FeH": 3, "Chi2": 4}

        # precomputation
        self.axis_vals = {}
        self.axis_idx = {}

        for ax in ["Teff", "logg", "FeH"]:
            vals = np.unique(self.grid[:, self.col[ax]])
            vals.sort()
            self.axis_vals[ax] = vals

            # map grid values -> integer indices
            idx = np.searchsorted(vals, self.grid[:, self.col[ax]])
            self.axis_idx[ax] = idx.astype(np.int32)

        self._proj1d_cache = {}
        self._proj2d_cache = {}

        # global chi2
        chi2 = self.grid[:, self.col["Chi2"]]
        self.chi2_min_global = float(np.nanmin(chi2))

        self.P3 = None
        if self.to_prob_3D:
            self._make_prob3d()

        self._regular_cache = {}
        self._oversample = float(oversample)
        self._nmax = int(nmax)

    # --------------------------
    # Internal helpers: PDFs
    # --------------------------
    def _build_fixed_mask(self, fixed):
        if not fixed:
            return np.ones(len(self.grid), dtype=bool)

        mask = np.ones(len(self.grid), dtype=bool)
        for ax, val in fixed.items():
            i = np.where(np.isclose(self.axis_vals[ax], val))[0]
            if len(i) == 0:
                return np.zeros(len(self.grid), dtype=bool)
            mask &= (self.axis_idx[ax] == i[0])
        return mask

    def _make_prob3d(self):
        """Compute full-grid probability P3 normalized."""
        chi2 = self.grid[:, self.col["Chi2"]].astype(float)
        dchi2 = chi2 - np.nanmin(chi2)
        with np.errstate(over="ignore", invalid="ignore"):
            P = np.exp(-0.5 * dchi2)
        P[np.isnan(P)] = 0.0
        s = P.sum()
        if s <= 0:
            P = np.zeros_like(P)
        else:
            P = P / s
        self.P3 = P
        return P
    # --------------------------
    # Regular grid utilities
    # --------------------------
    def _unique_sorted_axis(self, axis):
        vals = np.unique(self.grid[:, self.col[axis]])
        vals = np.sort(vals)
        return vals

    def _is_regular_axis(self, vals, rtol=1e-8):
        """Return True if axis values are regularly spaced (within rtol)."""
        if vals.size < 3:
            return True
        d = np.diff(vals)
        # consider only positive diffs
        dpos = d[d > 0]
        if dpos.size == 0:
            return True
        return np.allclose(dpos, dpos[0], rtol=rtol, atol=0)

    def _optimal_regular_grid(self, vals):
        """Build regular grid using dx = min positive difference (pas minimum observed)."""
        vals = np.sort(np.unique(vals))
        if vals.size < 2:
            return vals.copy()
        d = np.diff(vals)
        dpos = d[d > 0]
        if dpos.size == 0:
            return vals.copy()
        dx_min = float(np.min(dpos))
        xmin, xmax = float(vals[0]), float(vals[-1])
        # number of steps based on dx_min and oversampling factor
        n = int(np.ceil(self._oversample * (xmax - xmin) / dx_min)) + 1
        n = max(2, min(n, self._nmax))
        grid = np.linspace(xmin, xmax, n)
        return grid

    def _get_regular_axis(self, axis):
        """Return cached regular grid for axis or build it. Also return flag regular."""
        if axis in self._regular_cache:
            return self._regular_cache[axis]["vals"], self._regular_cache[axis]["regular"]

        vals = self._unique_sorted_axis(axis)
        regular = self._is_regular_axis(vals)
        if regular:
            grid = vals.copy()
        else:
            grid = self._optimal_regular_grid(vals)
        self._regular_cache[axis] = {"vals": grid, "regular": regular}
        return grid, regular

    # --------------------------
    # Regularize (interpolate) 1D and 2D on chi2 then convert if needed
    # --------------------------
    def _regularize_1d_chi2(self, X, P, axis):
        good = ~np.isnan(P)
        if not np.any(good):
            return np.array([]), np.array([])

        Xg = X[good]
        Pg = P[good]

        grid, is_reg = self._get_regular_axis(axis)

        f = interp1d(Xg, Pg, bounds_error=False, fill_value=np.nan)
        return grid, f(grid)

    def _regularize_2d_chi2(self, X, Y, P, xaxis, yaxis):
        Xr, _ = self._get_regular_axis(xaxis)
        Yr, _ = self._get_regular_axis(yaxis)

        # interpolation séparable (beaucoup plus rapide que griddata)
        Px = np.array([
            np.interp(Xr, X, row, left=np.nan, right=np.nan)
            for row in P
        ])

        Pr = np.array([
            np.interp(Yr, Y, Px[:, i], left=np.nan, right=np.nan)
            for i in range(len(Xr))
        ]).T

        return Xr, Yr, Pr

    # --------------------------
    # Projection on original grid (min chi2 over other dims)
    # --------------------------
    def _project_min_chi2_1d_raw(self, axis, fixed=None):
        idx = self.axis_idx[axis]
        chi2 = self.grid[:, self.col["Chi2"]]

        mask = self._build_fixed_mask(fixed)
        idx = idx[mask]
        chi2 = chi2[mask]

        n = len(self.axis_vals[axis])
        out = np.full(n, np.inf)
        np.minimum.at(out, idx, chi2)

        out[out == np.inf] = np.nan
        return self.axis_vals[axis], out


    def _project_min_chi2_2d_raw(self, xaxis, yaxis, fixed=None):
        xi = self.axis_idx[xaxis]
        yi = self.axis_idx[yaxis]
        chi2 = self.grid[:, self.col["Chi2"]]

        mask = self._build_fixed_mask(fixed)
        xi = xi[mask]
        yi = yi[mask]
        chi2 = chi2[mask]

        nx = len(self.axis_vals[xaxis])
        ny = len(self.axis_vals[yaxis])

        P = np.full((ny, nx), np.inf)
        np.minimum.at(P, (yi, xi), chi2)

        P[P == np.inf] = np.nan
        return self.axis_vals[xaxis], self.axis_vals[yaxis], P

    # --------------------------
    # Public projectors that return either chi2 or probability on regularized grid
    # --------------------------
    def project_1d(self, axis, to_prob=None, fixed=None):
        """
        Return Xr, P where:
         - Xr: regularized axis grid
         - P: if to_prob True -> normalized probability (sum=1); else -> chi2 (min chi2 per X)
        """
        key = (axis, bool(to_prob), None if fixed is None else tuple(sorted(fixed.items())))
        if key in self._proj1d_cache:
            return self._proj1d_cache[key]
        if to_prob is None:
            to_prob = self.to_prob

        Xraw, Pchi2_raw = self._project_min_chi2_1d_raw(axis, fixed=fixed)
        if Xraw.size == 0:
            return np.array([]), np.array([])

        Xr, Pchi2 = self._regularize_1d_chi2(Xraw, Pchi2_raw, axis)

        if Xr.size == 0:
            return np.array([]), np.array([])

        if to_prob:
            # convert chi2 -> probability with local minimum
            if np.all(np.isnan(Pchi2)):
                return Xr, np.zeros_like(Pchi2)
            dchi = Pchi2 - np.nanmin(Pchi2)
            with np.errstate(over="ignore", invalid="ignore"):
                P = np.exp(-0.5 * dchi)
            P[np.isnan(P)] = 0.0
            s = np.nansum(P)
            if s > 0:
                P = P / s
            else:
                P = np.zeros_like(P)
            self._proj1d_cache[key] = (Xr, P)
            return Xr, P
        else:
            self._proj1d_cache[key] = (Xr, P)
            return Xr, Pchi2

    def project_2d(self, xaxis, yaxis, to_prob=None, fixed=None):
        """
        Return Xr, Yr, P where:
         - Xr, Yr: regularized axis grids
         - P: if to_prob True -> normalized probability (sum ~1); else -> chi2 map (min chi2 per cell)
        """
        key = (xaxis, yaxis, bool(to_prob), None if fixed is None else tuple(sorted(fixed.items())))
        if key in self._proj2d_cache:
            return self._proj2d_cache[key]
        if to_prob is None:
            to_prob = self.to_prob

        Xraw, Yraw, Pchi2_raw = self._project_min_chi2_2d_raw(xaxis, yaxis, fixed=fixed)
        if Xraw.size == 0 or Yraw.size == 0:
            return np.array([]), np.array([]), np.array([[]])

        Xr, Yr, Pchi2 = self._regularize_2d_chi2(Xraw, Yraw, Pchi2_raw, xaxis, yaxis)

        if Xr.size == 0 or Yr.size == 0:
            return np.array([]), np.array([]), np.array([[]])

        if to_prob:
            if np.all(np.isnan(Pchi2)):
                return Xr, Yr, np.zeros_like(Pchi2)
            dchi = Pchi2 - np.nanmin(Pchi2)
            with np.errstate(over="ignore", invalid="ignore"):
                P = np.exp(-0.5 * dchi)
            P[np.isnan(P)] = 0.0
            s = np.nansum(P)
            if s > 0:
                P = P / s
            else:
                P = np.zeros_like(P)
            self._proj2d_cache[key] = (Xr, Yr, P)
            return Xr, Yr, P
        else:
            self._proj2d_cache[key] = (Xr, Yr, P)
            return Xr, Yr, Pchi2

    # --------------------------
    # Metrics (1D)
    # --------------------------
    def npeaks_1d(self, axis, threshold=0.2, fixed=None):
        _, P = self.project_1d(axis, to_prob=self.to_prob, fixed=fixed)
        if P.size == 0:
            return -99

        arr = P if self.to_prob else -P
        m = np.nanmax(arr)
        if not np.isfinite(m) or m <= 0:
            return 0

        peaks, _ = find_peaks(arr, height=threshold * m)
        return int(peaks.size)

    def peakratio_1d(self, axis, fixed=None):
        _, P = self.project_1d(axis, to_prob=self.to_prob, fixed=fixed)
        if P.size == 0:
            return -99

        mx = np.nanmax(P)
        return float(np.nanmean(P) / mx) if mx > 0 else -99

    def std_1d(self, axis, fixed=None):
        X, P = self.project_1d(axis, to_prob=self.to_prob, fixed=fixed)
        if P.size == 0:
            return -99

        if self.to_prob:
            W = P
        else:
            dchi = P - np.nanmin(P)
            W = np.exp(-0.5 * dchi)
            W[~np.isfinite(W)] = 0.0

        s = W.sum()
        if s <= 0:
            return -99
        W /= s
        mu = np.sum(X * W)
        return float(np.sqrt(np.sum(W * (X - mu) ** 2)))

    # --------------------------
    # Metrics (2D)
    # --------------------------
    def cov_2d(self, xaxis, yaxis, fixed=None):
        X, Y, P = self.project_2d(xaxis, yaxis, to_prob=True, fixed=fixed)
        if P.size == 0:
            return (-99, -99, -99, -99)

        # adapt proba to the regularized grid
        dx = np.gradient(X)
        dy = np.gradient(Y)
        DX, DY = np.meshgrid(dx, dy)
        W = P * DX * DY

        tot = W.sum()
        if tot <= 0:
            return (-99, -99, -99, -99)
        # make prob grid
        W /= tot
        Xg, Yg = np.meshgrid(X, Y)

        # Compute covariances
        mx = np.sum(W * Xg)
        my = np.sum(W * Yg)

        cov_xx = np.sum(W * (Xg - mx) ** 2)
        cov_yy = np.sum(W * (Yg - my) ** 2)
        cov_xy = np.sum(W * (Xg - mx) * (Yg - my))

        cov_mat = np.array([[cov_xx, cov_xy],
                            [cov_xy, cov_yy]])

        eig = np.linalg.eigvalsh(cov_mat)
        eig = np.maximum(eig, 1e-300)
        axis_ratio = float(np.sqrt(eig[1] / eig[0]))

        return float(cov_xx), float(cov_yy), float(cov_xy), axis_ratio

    def area_highprob_2d(self, xaxis, yaxis, threshold=0.5, fixed=None):
        X, Y, P = self.project_2d(xaxis, yaxis, to_prob=True, fixed=fixed)
        if P.size == 0:
            return -99.0

        Pmax = np.nanmax(P)
        if Pmax <= 0:
            return 0.0
        # adapt proba to the regularized grid
        dx = np.gradient(X)
        dy = np.gradient(Y)
        DX, DY = np.meshgrid(dx, dy)
        A = DX * DY
        
        valid = np.isfinite(P)
        high = P >= threshold * Pmax

        Atot = np.sum(A[valid])
        if Atot <= 0:
            return -99.0

        return float(np.sum(A[high & valid]) / Atot)


    # --------------------------
    # Global metrics
    # --------------------------
    def entropy(self):
        if self.P3 is None:
            # if P3 not precomputed, compute from chi2 on the fly
            self._make_prob3d()
        P = np.array(self.P3, dtype=float)
        mask = P > 0
        if np.sum(mask) == 0:
            return -99.0
        H = -np.sum(P[mask] * np.log(P[mask]))
        return float(H)

    def concentration(self, delta_chi2=1.0):
        chi2 = self.grid[:, self.col["Chi2"]].astype(float)
        if self.P3 is None:
            self._make_prob3d()
        mask = chi2 <= (self.chi2_min_global + float(delta_chi2))
        if np.sum(mask) == 0:
            return 0.0
        return float(np.sum(self.P3[mask]))

    def delta_chi2(self):
        chi2 = np.array(self.grid[:, self.col["Chi2"]], dtype=float)
        if chi2.size < 2:
            return -99.0
        sorted_chi = np.sort(chi2)
        return float(sorted_chi[1] - sorted_chi[0])

    # --------------------------
    # Main entry: compute all features
    # --------------------------
    def compute_all(self, fixed=None):
        """
        Compute all requested features and return as dict.
        fixed: optional dict to restrict projection (e.g. {"FeH":0.0})
        """
        feats = {}

        for ax in ["Teff", "logg", "FeH"]:
            feats[f"npeaks_{ax}"] = self.npeaks_1d(ax, threshold=0.2, fixed=fixed)
            feats[f"peakratio_{ax}"] = self.peakratio_1d(ax, fixed=fixed)
            feats[f"std_{ax}"] = self.std_1d(ax, fixed=fixed)

        pairs = [("Teff", "logg"), ("Teff", "FeH"), ("logg", "FeH")]
        for x, y in pairs:
            cov_xx, cov_yy, cov_xy, axis_ratio = self.cov_2d(x, y, fixed=fixed)
            feats[f"cov_{x}_{y}_xx"] = cov_xx
            feats[f"cov_{x}_{y}_yy"] = cov_yy
            feats[f"cov_{x}_{y}_xy"] = cov_xy
            feats[f"axisratio_{x}_{y}"] = axis_ratio
            feats[f"area_{x}_{y}_p50"] = self.area_highprob_2d(x, y, threshold=0.5, fixed=fixed)
            feats[f"area_{x}_{y}_p90"] = self.area_highprob_2d(x, y, threshold=0.9, fixed=fixed)

        feats["entropy"] = self.entropy()
        feats["concentration_dchi2_1"] = self.concentration(delta_chi2=1.0)
        feats["delta_chi2"] = self.delta_chi2()
        return feats


#--- sample analyzer class
# Run star_pdf_analyzer metrics on a sample

from joblib import Parallel, delayed

class SAMPLE_ANALYZER:
    def __init__(self, chi2_path, sed_grid_path, nrows=None, n_jobs=1):

        # --- Load raw file as strings ---
        raw = pd.read_csv(chi2_path, sep=r"\s+", header=None,
                        dtype="str", comment="#", nrows=nrows)
        raw = raw.to_numpy()
        self.ids = raw[:, 0].astype(np.int64)
        self.chi2_values = raw[:, 1:].astype(float)

        # SED grid (constant for all sources)
        self.sed_grid_obj = SED_GRID(sed_grid_path)
        self.sed_grid = self.sed_grid_obj.build()

        self.n_jobs = int(n_jobs)

    def best_fit_distribution(self, n_rows=None):

        star_pdf = STAR_PDF(self.sed_grid)
        data = []

        N = n_rows if n_rows is not None else len(self.ids)

        for i in range(N):
            _ = star_pdf.build_chi2_grid(self.chi2_values[i])
            df = star_pdf.best_output()
            data.append([self.ids[i]] + df)

        columns = [
            'source_id',
            'best_sed_id', 'Teff_best', 'logg_best', 'FeH_best', 'chi2_best',
            'sec_best_sed_id', 'sec_Teff_best', 'sec_logg_best', 'sec_FeH_best', 'sec_chi2_best'
        ]

        self.sample_df = pd.DataFrame(data, columns=columns)
        return self.sample_df

    def stat_feature_distribution(self, n_rows=None):

        N = n_rows if n_rows is not None else len(self.ids)
        star_pdf = STAR_PDF(self.sed_grid)

        # ---- worker function 
        def onesource(i):
            chi2_grid = star_pdf.build_chi2_grid(self.chi2_values[i])
            analyzer = STAR_PDF_ANALYZER(chi2_grid)
            feats = analyzer.compute_all()
            feats["source_id"] = self.ids[i]
            return feats

        # ---- run
        if self.n_jobs == 1:
            rows = [onesource(i) for i in range(N)]
        else:
            rows = Parallel(n_jobs=self.n_jobs, backend="loky")(
                delayed(onesource)(i) for i in range(N)
            )

        # ---- build dataframe with fixed column order
        columns = ["source_id"] + [k for k in rows[0].keys() if k != "source_id"]
        self.feature_df = pd.DataFrame(rows, columns=columns)

        return self.feature_df

    def plot_best(self, axis='chi2_best'):
        plt.figure()
        plt.hist(self.sample_df[axis].values, bins=50)
        plt.xlabel(axis)
        plt.ylabel("N")
        plt.show()

#--- SPLIT ---

from scipy.integrate import quad
class SPLIT_TRAINING:
    """
    Prepare the statistical feature distributions for the classification.
    """
    def list_to_pdf(self, data, interp=False, bins=100, recenter_bins=True, xrange=None, smooth_eps=1e-9, normalize=True):
        """
        Compute the pdf from a data array.

        Parameters
        ----------
        data: list or np.darray
            Array to compute the pdf from.
        interp: Bool or str.
            Interpolation method to build the pdf.
        bins: int or list or np.darray.
            Bins used to compute the histogram used in <scipy.interpolate.{interp}> methods.
        recenter_bins: bool.
            If True, recenter the bins over the bin edges. Useful when using a list in bins
        xrange: tuple.
            Range to compute the pdf within.
        smooth_eps: float.
            Zero probability value.
        
        Returns
        -------
        pdf: pdf function.
        """
        if interp in ('scott', 'silverman') or isinstance(interp,(float,int)) and interp!=False:
            kde = gaussian_kde(data, bw_method=interp)
            def pdf(x):
                return kde(x)
            return pdf

        if recenter_bins==True and not isinstance(bins, int):
            new_bins = shift_bins(bins)
            counts, _ = np.histogram(data, bins=new_bins, range=xrange, density=True)
            x = bins
        else:
            counts, bin_edges = np.histogram(data, bins=bins, range=xrange, density=True)
            x=0.5 * (bin_edges[:-1] + bin_edges[1:])
        
        counts = counts + smooth_eps
        # --- Interpolation ---
        if interp == 'pchip':
            pdf_raw = sip.PchipInterpolator(x, counts, extrapolate=False)
        else:
            pdf_raw = sip.interp1d(x, counts, kind=interp, bounds_error=False, fill_value=0.0)

        # --- Normalization ---
        if normalize==True:
            xmin, xmax = min(x), max(x)
            Z, _ = quad(pdf_raw, xmin, xmax, limit=1000)

            if Z <= 0 or not np.isfinite(Z):
                raise ValueError("PDF normalization failed (integral <= 0).")

            def pdf(x):
                return pdf_raw(x) / Z
        else:
            pdf = pdf_raw

        return pdf

    def pdfs_from_arrays(self, arrays, **ltpkwargs):
        """
        Return the pdf from two arrays and compute the weight if specified.

        Parameters
        ----------
        arrays: tuple of (list or np.darray)
            Arrays to compute the pdfs from.
        ltpkwargs: args of <list_to_pdf> method.
        Returns
        -------
        pdfs: list of pdf methods.
        """
        pdfs = []
        for dist in arrays:
            pdf = self.list_to_pdf(dist, **ltpkwargs)
            pdfs.append(pdf)

        return pdfs

    def compute_weight(self, source_feats, source_labels, method='f_classif'):
        """
        Compute weight by comparing feature distributions between several classes.
        The weight determine if the feature is relevant to be used to compare the classes.

        Parameters
        ----------
        source_feats: list of (list or np.darray)
            Array of dimension N_sources x N_features.
            Each element corresponds to a source, and each element of a row is a statistical feature.
        source_labels: list of (str, int...).
            The list of label corresponding to each source.
        method: str.
            The method to use to compute weight (fisher or f_classic)
        Returns
        -------
        List of weights of length N_features.
        """
        source_feats = np.asarray(source_feats)
        source_labels = np.asarray(source_labels)

        if method == 'fisher':
            classes = np.unique(source_labels)
            n_features = source_feats.shape[1]
            weights = np.zeros(n_features)
            global_mean = np.mean(source_feats, axis=0)
            for j in range(n_features):
                num = 0.0  # inter-class variance
                den = 0.0  # intra-class variance
                for c in classes:
                    mask = source_labels == c
                    x_c = source_feats[mask, j]
                    if len(x_c) < 2:
                        continue
                    mean_c = np.mean(x_c)
                    var_c = np.var(x_c)
                    num += len(x_c) * (mean_c - global_mean[j])**2
                    den += len(x_c) * var_c
                weights[j] = num / den if den > 0 else 0.0

        elif method == 'f_classif':
            f_statistic, _ = f_classif(source_feats, source_labels)
            weights = f_statistic

        else:
            raise NotImplementedError(
                f"Weight method '{method}' not recognized."
            )

        # --- normalisation ---
        weights = np.nan_to_num(weights, nan=0.0, posinf=0.0)
        if np.sum(weights) > 0:
            weights /= np.sum(weights)
        return weights

    def pdfs_from_sample(self, training_df, label_column, scols=None, exclude_scols=False, 
                        compute_weight=False, weight_method='fisher', selectKbest=None, 
                        custo_bins_interp=None, **ltpkwargs):
        """
        Compute the pdf for each class, for each features.

        Parameters
        ----------
        dfs: pandas.dataframe.
            df to compute the pdfs from. It must contain at leats a label column,
            and the features columns to compute the pdfs from. Each row is a source.
        label_column: str.
            Column holding the label attributed to each source. 
        scols: str or list of.
            Selected columns to compute the pdfs from is [{col_name}],
            or to ignore if ~[{col_name}].
        exclude_scols: bool.
            If True, scols are now the columns to exclude from the computations.
        compute_weight: bool.
            If true, return the weight of each statistical feature.
        custo_bins_interp: tuple (str, list) or list of.
            If specified, for each col in custo_bins_interp[0], list_to_pdf uses the bins in custo_bins_interp[1].
        ltpkwargs: args of <list_to_pdf> method.
        Returns
        -------
        pdfs_dict: Dict
            pdf functions, with weights if True.
        """

        # --- Select feature columns ---
        if scols is not None:
            if exclude_scols==True:
                cols = [c for c in training_df.columns
                            if c not in scols
                            and c != label_column]
            elif exclude_scols==False:
                cols = [c for c in training_df.columns
                            if c in scols
                            and c != label_column]
                
        else:
            cols = training_df.drop(columns=[label_column]).columns.tolist()
        # --- extract features / labels ---
        features_values = training_df[cols].values.astype(float) # Nsources(rows) x Nfeatures(cols)
        source_labels = training_df[label_column].values # Nsources
        unique_labels = np.unique(source_labels)
        # --- compute weights ---
        if compute_weight==True:
            weights = self.compute_weight(features_values, source_labels, method=weight_method)
            # select the k best weight features from weights
            if isinstance(selectKbest, int):
                idx = np.argsort(weights)[::-1][:selectKbest]
                cols = [cols[i] for i in idx]
                weights = weights[idx]
                features_values = features_values[:, idx]

                weights /= np.sum(weights)
        else:
            weights = np.ones(len(cols)) / len(cols)

        # --- Build pdf dictionary ---
        # len(source_labels)*N_features = N_pdfs
        custo_map = None
        if custo_bins_interp is not None:
            custo_map = {}
            for item in custo_bins_interp:
                if len(item) != 3:
                    raise ValueError(
                        "Each element of custo_bins_interp must be (col, bins, interp)"
                    )
                col, bins, interp = item
                custo_map[col] = {"bins": bins, "interp": interp}

        pdf_dict = {}
        for i, col in enumerate(cols):
            pdfs = {}

            for label in unique_labels:
                mask = source_labels == label
                data = features_values[mask, i]

                local_kwargs = ltpkwargs
                # --- Apply custom bins / interp if specified ---
                if custo_map is not None and col in custo_map:
                    if custo_map[col]["bins"] is not None:
                        local_kwargs["bins"] = custo_map[col]["bins"]
                    if custo_map[col]["interp"] is not None:
                        local_kwargs["interp"] = custo_map[col]["interp"]
                pdfs[label] = self.list_to_pdf(data, **local_kwargs)
            pdf_dict[col] = {"pdfs": pdfs, "weight": weights[i]}

        return pdf_dict
        

class SPLIT:
    """
    Class where appends the main part of split. 

    Init with the pdf dict created with <SPLIT_TRAINING>.
    """
    def __init__(self, pdf_dict):
        self.pdf_dict = pdf_dict
        self.features = list(pdf_dict.keys())
        self.classes = list(
            next(iter(pdf_dict.values()))['pdfs'].keys()
        )

    def feature_loglikelihood(self, x, feature, eps=1e-12):
        """
        Compute log p(x | class) for a given feature.
        """
    
        pdfs = self.pdf_dict[feature]['pdfs']

        logL = []

        for c in self.classes:
            p = pdfs[c](x)
            p = np.maximum(p, eps)
            logL.append(np.log(p))

        return np.array(logL)

    def class_prob(self, onesource, priors=None, use_weights=False, eps=1e-12):
        """
        Compute the probability of a source to belong to each class using <self.post_prob>.

        Parameters:
        -----------
        onesource: the source to evaluate the class on. Can be a list of features corresponding to the pdfs

        priors: None or list of float of length = len(pdfs).
            If specified, take priors in account.
        use_weights: bool
            If True, uses the feature weights specified in pdfs.
        
        Return:
        -------
        probs: list of floats
            List of probabilies to belong to the classes.
            
        """
        n_classes = len(self.classes)

        if priors is None:
            priors = np.ones(n_classes) / n_classes
        else:
            priors = np.asarray(priors)

        log_post = np.log(priors)

        for feature in self.features:

            x = onesource[feature]

            logL_feat = self.feature_loglikelihood(x, feature, eps)

            if use_weights:
                w = self.pdf_dict[feature]['weight']
                log_post += w * logL_feat
            else:
                log_post += logL_feat

        # softmax stable
        max_log = np.max(log_post)
        exp_log = np.exp(log_post - max_log)
        probs = exp_log / np.sum(exp_log)

        return probs

    def classify_sample(self, df, priors=None, use_weights=False):

        results = []

        for _, row in df.iterrows():
            probs = self.class_prob(
                row,
                priors=priors,
                use_weights=use_weights
            )
            results.append(probs)

        posterior = np.vstack(results)
        
        return pd.DataFrame(
            posterior,
            columns=[f"P_{c}" for c in self.classes],
            index=df.index
        )

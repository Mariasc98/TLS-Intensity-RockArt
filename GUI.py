"""
TLS 2D Raster Transform GUI
----------------------------

Tkinter GUI to run and compare PCA / FastICA / NMF / MNF on 4 band raster dataset:
R, G, B, Intensity (single-band GeoTIFFs)

This script contains:
  - GUI logic
  - I/O helpers for viewing/saving/export
  - Visualization

The transformation algorithms are imported from:
  - PCA.py
  - fastICA.py
  - NMF.py
  - MNF.py

Note that while the number of components can be selected within the GUI, the parameters of 
PCA, ICA, MNF and NMF used for the analysis within the MA Thesis are set as the defaults.
Therefore, to adapt the parameters to specific case studies these need to be changed within 
this GUI script.

Script developed by Maria Sotomayor Chicote
"""

import os
import csv
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import rasterio

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Import existing algorithms
from PCA import pca_from_geotiff
from fastICA import fastica_from_geotiff
from NMF import nmf_from_geotiff
from MNF import mnf_from_geotiff

# Utils

def stretch01(arr, p_low=2, p_high=98):
    m = np.isfinite(arr)
    if not np.any(m):
        return arr
    lo, hi = np.nanpercentile(arr[m], [p_low, p_high])
    return np.clip((arr - lo) / (hi - lo + 1e-12), 0, 1)


def read_band(path):
    with rasterio.open(path) as src:
        band = src.read(1).astype(np.float32)
        profile = src.profile
        nodata = src.nodata
    return band, profile, nodata


def write_multiband_tif(path, arrays, profile):
    arrays = np.asarray(arrays)
    if arrays.ndim != 3:
        raise ValueError("arrays must have shape (n_bands, H, W)")

    out_profile = profile.copy()
    out_profile.update(
        dtype="float32",
        count=int(arrays.shape[0]),
        nodata=np.nan,
        compress="lzw"
    )

    with rasterio.open(path, "w", **out_profile) as dst:
        for k in range(arrays.shape[0]):
            dst.write(arrays[k].astype(np.float32), k + 1)

# GUI

class TransformApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("TLS 2D Raster Transformations: PCA / ICA / NMF / MNF")
        self.geometry("1250x880")

        # File inputs
        self.paths = {
            "R": tk.StringVar(value=""),
            "G": tk.StringVar(value=""),
            "B": tk.StringVar(value=""),
            "Intensity": tk.StringVar(value=""),
        }

        # Transform toggles
        self.apply_pca = tk.BooleanVar(value=True)
        self.apply_ica = tk.BooleanVar(value=False)
        self.apply_nmf = tk.BooleanVar(value=False)
        self.apply_mnf = tk.BooleanVar(value=False)

        # Parameters
        self.n_components = tk.IntVar(value=4)

        # Stored results
        self.base_profile = None
        self.originals = None
        self.last_pca = None
        self.last_ica = None
        self.last_nmf = None
        self.last_mnf = None

        self._build_ui()

    def _build_ui(self):
        left = ttk.Frame(self)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)

        ttk.Label(left, text="Input rasters (single-band GeoTIFFs)", font=("Segoe UI", 10, "bold")).pack(anchor="w")

        for ch in ["R", "G", "B", "Intensity"]:
            row = ttk.Frame(left)
            row.pack(fill=tk.X, pady=3)
            ttk.Label(row, text=ch, width=10).pack(side=tk.LEFT)
            ttk.Entry(row, textvariable=self.paths[ch], width=52).pack(side=tk.LEFT, padx=4)
            ttk.Button(row, text="Browse", command=lambda c=ch: self.browse(c)).pack(side=tk.LEFT)

        opt = ttk.LabelFrame(left, text="Options")
        opt.pack(fill=tk.X, pady=10)

        ttk.Label(opt, text="Number of components:").grid(row=0, column=0, sticky="w", padx=6, pady=6)
        ttk.Spinbox(opt, from_=1, to=4, textvariable=self.n_components, width=6).grid(row=0, column=1, padx=6, pady=6)

        ttk.Checkbutton(opt, text="Apply PCA", variable=self.apply_pca).grid(row=1, column=0, sticky="w", padx=6)
        ttk.Checkbutton(opt, text="Apply ICA (FastICA)", variable=self.apply_ica).grid(row=2, column=0, sticky="w", padx=6)
        ttk.Checkbutton(opt, text="Apply NMF", variable=self.apply_nmf).grid(row=3, column=0, sticky="w", padx=6)
        ttk.Checkbutton(opt, text="Apply MNF", variable=self.apply_mnf).grid(row=4, column=0, sticky="w", padx=6)

        btns = ttk.Frame(left)
        btns.pack(fill=tk.X, pady=8)

        ttk.Button(btns, text="Run Analysis", command=self.on_run).pack(fill=tk.X, pady=4)
        ttk.Button(btns, text="Save Components to Multiband GeoTIFF", command=self.on_save_components).pack(fill=tk.X, pady=4)
        ttk.Button(btns, text="Export PCA eigenvalues (CSV)", command=self.on_export_pca).pack(fill=tk.X, pady=4)
        ttk.Button(btns, text="Export MNF eigenvalues (CSV)", command=self.on_export_mnf).pack(fill=tk.X, pady=4)

        # Plot area
        self.fig = plt.Figure(figsize=(9, 10), constrained_layout=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def browse(self, channel):
        path = filedialog.askopenfilename(
            title=f"Select {channel} raster",
            filetypes=[("GeoTIFF", "*.tif;*.tiff"), ("All files", "*.*")]
        )
        if path:
            self.paths[channel].set(path)

    def _check_inputs(self):
        r = self.paths["R"].get().strip()
        g = self.paths["G"].get().strip()
        b = self.paths["B"].get().strip()
        i = self.paths["Intensity"].get().strip()
        if not (r and g and b and i):
            raise ValueError("Please select all four rasters: R, G, B, Intensity.")
        return r, g, b, i

    def on_run(self):
        try:
            r, g, b, i = self._check_inputs()
            n = int(self.n_components.get())
            if n < 1 or n > 4:
                raise ValueError("Number of components must be between 1 and 4 (because you have 4 bands).")

            # Reset last results
            self.last_pca = None
            self.last_ica = None
            self.last_nmf = None
            self.last_mnf = None

            # Read originals for display + store base profile
            R, profile, _ = read_band(r)
            G, _, _ = read_band(g)
            B, _, _ = read_band(b)
            I, _, _ = read_band(i)

            self.base_profile = profile
            self.originals = [("R", R), ("G", G), ("B", B), ("Intensity", I)]

            # Run selected transforms (from imports). You can change the parameters here.

            if self.apply_pca.get():
                pca_img, pca_model, eigenvalues, ratio = pca_from_geotiff(
                    r_path=r, g_path=g, b_path=b, i_path=i,
                    out_path=None,
                    n_components=n,
                    whiten=False,
                    svd_solver="auto",
                    random_state=None
                )
                self.last_pca = {"img": pca_img, "model": pca_model, "eigenvalues": eigenvalues, "ratio": ratio}

            if self.apply_ica.get():
                ica_img, mixing, unmixing = fastica_from_geotiff(
                    r_path=r, g_path=g, b_path=b, i_path=i,
                    algorithm="parallel",
                    fun="cube",
                    tol=1e-5,
                    whiten="unit-variance",
                    max_iter=2000,
                    n_components=n,
                    random_state=0,
                    out_path=None
                )
                self.last_ica = {"img": ica_img, "mixing": mixing, "unmixing": unmixing}

            if self.apply_nmf.get():
                nmf_img, nmf_model, Wmat, Hmat = nmf_from_geotiff(
                    r_path=r, g_path=g, b_path=b, i_path=i,
                    n_components=n,
                    init="nndsvd",
                    max_iter=2000,
                    random_state=0,
                    alpha_W=0.0,
                    alpha_H=0.0,
                    l1_ratio=0.0,
                    out_path=None
                )
                self.last_nmf = {"img": nmf_img, "model": nmf_model, "W": Wmat, "H": Hmat}

            if self.apply_mnf.get():
                mnf_img, mnf_eigs, W_noise, pca_model = mnf_from_geotiff(
                    r_path=r, g_path=g, b_path=b, i_path=i,
                    n_components=n,
                    out_path=None
                )
                self.last_mnf = {"img": mnf_img, "eigenvalues": mnf_eigs, "W_noise": W_noise, "model": pca_model}

            self.plot_results()

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def plot_results(self):
        if self.originals is None:
            return

        self.fig.clf()

        blocks = []
        blocks.append(("Original", [arr for _, arr in self.originals]))

        if self.last_pca:
            blocks.append(("PCA", [self.last_pca["img"][k] for k in range(self.last_pca["img"].shape[0])]))
        if self.last_ica:
            blocks.append(("ICA", [self.last_ica["img"][k] for k in range(self.last_ica["img"].shape[0])]))
        if self.last_nmf:
            blocks.append(("NMF", [self.last_nmf["img"][k] for k in range(self.last_nmf["img"].shape[0])]))
        if self.last_mnf:
            blocks.append(("MNF", [self.last_mnf["img"][k] for k in range(self.last_mnf["img"].shape[0])]))

        n_rows = len(blocks)
        n_cols = max(len(blocks[0][1]), max(len(b[1]) for b in blocks))

        gs = self.fig.add_gridspec(n_rows, n_cols, hspace=0.25, wspace=0.15)

        for r_idx, (label, images) in enumerate(blocks):
            for c_idx in range(n_cols):
                ax = self.fig.add_subplot(gs[r_idx, c_idx])
                ax.axis("off")
                if c_idx < len(images):
                    ax.imshow(stretch01(images[c_idx]), cmap="gray")
                    if label == "Original":
                        ch_name = self.originals[c_idx][0] if c_idx < len(self.originals) else f"{c_idx + 1}"
                        ax.set_title(f"{label}: {ch_name}")
                    else:
                        ax.set_title(f"{label} {c_idx + 1}")

        self.canvas.draw()

    def on_save_components(self):
        if self.base_profile is None:
            messagebox.showinfo("No data", "Run analysis first.")
            return

        folder = filedialog.askdirectory(title="Select folder to save component GeoTIFF(s)")
        if not folder:
            return

        saved = []

        try:
            if self.last_pca:
                outp = os.path.join(folder, "PCA_components.tif")
                write_multiband_tif(outp, self.last_pca["img"], self.base_profile)
                saved.append(outp)

            if self.last_ica:
                outp = os.path.join(folder, "ICA_components.tif")
                write_multiband_tif(outp, self.last_ica["img"], self.base_profile)
                saved.append(outp)

            if self.last_nmf:
                outp = os.path.join(folder, "NMF_components.tif")
                write_multiband_tif(outp, self.last_nmf["img"], self.base_profile)
                saved.append(outp)

            if self.last_mnf:
                outp = os.path.join(folder, "MNF_components.tif")
                write_multiband_tif(outp, self.last_mnf["img"], self.base_profile)
                saved.append(outp)

        except Exception as e:
            messagebox.showerror("Write error", f"Failed to save GeoTIFF(s):\n{e}")
            return

        if saved:
            messagebox.showinfo("Saved", "Saved files:\n" + "\n".join(saved))
        else:
            messagebox.showinfo("Nothing to save", "No transform results available to save.")

    def on_export_pca(self):
        if not self.last_pca:
            messagebox.showinfo("No PCA results", "Run analysis with PCA enabled first.")
            return

        folder = filedialog.askdirectory(title="Select folder to save PCA eigenvalues CSV")
        if not folder:
            return

        out_csv = os.path.join(folder, "pca_eigenvalues.csv")
        eigenvalues = self.last_pca["eigenvalues"]
        ratio = self.last_pca["ratio"]

        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["component_index", "eigenvalue", "explained_variance_ratio"])
            for idx, (ev, er) in enumerate(zip(eigenvalues, ratio), start=1):
                w.writerow([idx, float(ev), float(er)])

        messagebox.showinfo("Saved", f"PCA eigenvalues saved to:\n{out_csv}")

    def on_export_mnf(self):
        if not self.last_mnf:
            messagebox.showinfo("No MNF results", "Run analysis with MNF enabled first.")
            return

        folder = filedialog.askdirectory(title="Select folder to save MNF eigenvalues CSV")
        if not folder:
            return

        out_csv = os.path.join(folder, "mnf_eigenvalues.csv")
        eigenvalues = self.last_mnf["eigenvalues"]

        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["component_index", "mnf_eigenvalue"])
            for idx, ev in enumerate(eigenvalues, start=1):
                w.writerow([idx, float(ev)])

        messagebox.showinfo("Saved", f"MNF eigenvalues saved to:\n{out_csv}")


if __name__ == "__main__":
    app = TransformApp()
    app.mainloop()




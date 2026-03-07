# TLS Intensity Analysis for Rock Art Enhancement

**Version 1.0**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview 

This repository contains the Python-based workflows developed for the MA Thesis in [Digital and Computational Archaeology](https://dca.uni-koeln.de/) at the University of Cologne: **"From 3D Point Clouds to Raster Images: A Python-based Approach for Rock Art Enhancement through Intensity Data of Terrestrial Laser Scanning (TLS)"**.

The workflows enable the analysis of TLS intensity data for rock art documentation, enhancement, and interpretation through two complementary approaches: 3D spectral index analysis and 2D raster-based statistical transformations.

Within the framework of the MA Thesis, these workflows were tested on TLS datasets from Peñablanca (Philippines) and Wadi Sura
(Egypt), demonstrating that intensity values can reveal faint or invisible motifs, separate superimposed pigments, and differentiate
original rock art from later graffiti. The results obtained within this work are not universally generalizable, as intensity behaviour is 
context-dependent and relies on a series of external factors such as acquisition parameters, wavelenght selection, surface properties 
or environmental conditions. Therefore, more testing of these workflows in further case studies is needed. 
The aim of this repository is to make available the scripts developed within the MA Thesis to enable reproducibility and further testing of these
methods, contributing future research in TLS intensity analysis.  

---

## Context

Terrestrial Laser Scanning has become a fundamental tool for the non-invasive documentation of rock art, primarily due to its capacity to capture high-resolution
three-dimensional geometry. However, beyond geometric information, TLS also records intensity values that reflect the interaction between laser wavelenght and
surface material properties. Recent research demonstrated the ability of near-infrared (NIR) TLS intensity values to reveal hidden carbon-based pigments
(Jalandoni, Winans et al. 2021). However, despite their analytical potential, TLS intensity values remain underexplored in rock art research and are rarely integrated into standardized, reproducible analytical workflows. 

Within this framework, the MA Thesis mentioned above aimed to address this research gap by examining TLS intensity values as an analytical and archaeological source of information, rather than merely a technical product of three-dimensional recording. 

---

## Features

### 3D Spectral Index Analysis (`3D_indices.py`)

Compute and visualize spectral indices directly on 3D point clouds:
- **IRR** (Intensity-to-Red Ratio)
- **IRD** (Intensity-Red Difference)
- **RNI** (Red-to-Intensity Ratio)
- **NRAI** (Normalized Rock Art Index)

### 2D Raster-based Analysis

Four statistical transformation methods are available for analyzing four-band raster datasets (R,G,B,Intensity):

- **Principal Component Analysis** (PCA)- (`PCA.py`)
- **Independent Component Analysis** (ICA)- (`fastICA.py`)
- **Minimum Noise Fraction** (MNF)- (`MNF.py`)
- **Non-Negative Matrix Factorization** (NMF)- (`NMF.py`)

### Graphical User Interface (GUI) (`GUI.py`)

A unified interface that integrates all four transformation methods, allowing the user to:

- Load single-band GeoTIFFs for R, G, B, and Intensity channels.
- Select multiple transformation methods simultaneously.
- Adjust the number of components for each method.
- Compare outputs from all methods.
- Export results as multi-band GeoTIFFs.
- Save eigenvalues and explained variance as CSV files.

---

## Repository Structure

```
TLS-Intensity-RockArt/
│
├── README.md                # This file
├── LICENSE                  # MIT License
├── environment.yml          # Conda environment specification
├── 3D_indices.py            # 3D spectral index calculations (IRR, IRD, RNI, NRAI)
├── GUI.py                   # Graphical interface for 2D raster transformations
├── PCA.py                   # Principal Component Analysis
├── fastICA.py               # Independent Component Analysis (FastICA algorithm)
├── MNF.py                   # Minimum Noise Fraction
└── NMF.py                   # Non-Negative Matrix Factorization
```
---

## Quick Start

### Prerequisites

- Anaconda3

### Installation

1. **Clone the repository**
``` 
git clone https://github.com/Mariasc98/TLS-Intensity-RockArt
```

2. **Create the conda environment**
```
cd TLS-Intensity-RockArt
```

```
conda env create -f environment.yml
```

3. **Activate the environment**
```
conda activate intensity_rockart
```

---

## Input Data Formats

### 2D Workflow Input (rasters)

Provide aligned single-band GeoTIFFs:

- R.tif (Red)
- G.tif (Green)
- B.tif (Blue)
- Intensity.tif (TLS intensity)

Important:

- All rasters must share the same extent, resolution, and pixel alignment.
- NoData values are supported (the workflow uses a valid-pixel mask).
- Best practice is to generate them using the same rasterization parameters (cell size, interpolation,etc). For the
rasterization of the point cloud data software such as [CloudCompare](https://www.cloudcompare.org/) can be used.

### 3D Workflow Input (point cloud)

.e57 point cloud with XYZ, RGB and intensity values.

---

## How to Run

### 1. Run the 2D scripts without the GUI

You can run the transformation algorithms individually via the standalone scripts. First, edit the input paths and parameters
at the end of each script. Then run the script.

Alternatively you can run directly the selected script from the activated conda environment:
```
python PCA.py
python fastICA.py
python MNF.py
python NMF.py
```

### 2. Run the GUI (2D raster tool)

Within the GUI you can select the transformation algorithms to run, number of components as well as optional outputs. However, the
specific parameters of each transform should be changed manually within the script. If you do not change them, the parameters used 
within the analysis carried out for the MA Thesis are the default. 
You can run the script directly, or from the activated environment:
```
python GUI.py
```
Workflow:

1. **Load Bands** (4 GeoTIFFs).

2. **Select Methods** (PCA/ICA/MNF/NMF). You can select one or multiple methods at once for side-by-side comparison.

3. **Set number of components**.

4. **Run** to generate outputs and preview component layers.

5. **Export** (optional):
   - GeoTIFF component stacks per method.
   - CSV for PCA/MNF (eigenvalues/ explained variance).

### 3. Run the 3D indices workflow

Inside the script, set:

- the path to your .e57 file.
- the index to compute/visualize.

---

## Outputs

### 2D methods

- On-screen previews of components (grayscale with consistent stretching for comparability).
- Optional exported GeoTIFFs:
  - a multi-band GeoTIFF, one band per component per method.
- Optional CSV exports:
  - eigenvalues/ explained variance ratios for PCA and MNF.

### 3D indices

- point cloud visualization with computed scalar field (index values).

---

## License

This project is licensed under the MIT License - see the `LICENSE`file for details.

---

## Citation

If you use this code, please cite as follows:
```
@software{sotomayor2026tlsintensity,
  author = {Sotomayor Chicote, Maria},
  title = {TLS Intensity Analysis for Rock Art Enhancement},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/Mariasc98/TLS-Intensity-RockArt}
}
```

## References

Jalandoni, A, W. Winans, and M. Willis. 2021. "Intensity values of Terrestrial Laser Scans 
Reveal Hidden Black Rock Art Pigment". *Remote Sensing* 13 (7): 1357. DOI: https//:doi.org/10.3390/rs13071357

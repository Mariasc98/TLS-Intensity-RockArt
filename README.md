# TLS Intensity Analysis for Rock Art Enhancement

**Version 1.0**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview 

This repository contains the Python-based workflows developed for the MA Thesis in Digital and Computational Archaeology
at the University of Cologne: **"From 3D Point Clouds to Raster Images: A Python-based Approach for Rock Art Enhancement
through Intensity Data of Terrestrial Laser Scanning (TLS)"**

The workflows enable the analysis of TLS intensity data for rock art documentation, enhancement, and interpretation. Two
complementary analytical approaches are implemented:

1. **3D Point Cloud Analysis** -
Spectral indices applied directly to 3D point clouds (`3D_indices.py`)
2. **2D Raster-based Analysis** -
Statistical transformations (PCA, ICA, MNF, NMF) on four-band raster datasets (R,G,B,Intensity), with two ways to run them:
-**Individual scripts** for each transformation method (`PCA.py´), (`fastICA.py´), (`MNF.py´) and (`NMF.py´)
-**Graphical User Interface (GUI)** (`GUI.py´) for easy application and comparison of all methods

Within the framework of the MA Thesis, these workflows were tested on TLS datasets from Peñablanca (Philippines) and Wadi Sura
(Egypt), demonstrating that intensity values can reveal faint or invisible motifs, separate superimposed pigments, and differentiate
original rock art from later graffiti. The results obtained within this work are not universally generalizable, as intensity behaviour is 
context-dependent and relies on a series of external factors such as acquisition parameters, wavelenght selection, surface properties 
or environmental conditions. Therefore, more testing of these workflows in further case studies is needed. 
The aim of this repository is to make available these scripts to enable reproducibility and further testing of these
methods, contributing future research in TLS intensity analysis.  also to other case studies...

## Features

### 1. 3D Spectral Index Analysis (`3D_indices.py`)

Compute and visualize spectral indices directly on 3D point clouds:
- **IRR** (Intensity-to-Red Ratio)
- **IRD** (Intensity-Red Difference)
- **RNI** (Red-to-Intensity Ratio)
- **NRAI** (Normalized Rock Art Index)

### 2. 2D Raster-based Analysis

Four statistical transformation methods are available for analyzing four-band rasters (R,G,B,Intensity):

| Method | Description | Script |







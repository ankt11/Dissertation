# MSc Dissertation  
## Reconstructing Ocean Heat Content Anomalies with Random Forest 

---

## Table of Contents

- [Project Overview](#-project-overview)
- [Research Motivation](#-research-motivation)
- [Research Question](#-research-question)
- [Technologies Used](#-technologies-used)
- [Methodological Approach](#-methodological-approach)
- [How the Project Developed](#-how-the-project-developed)
- [Challenges & Limitations](#-challenges--limitations)
- [Problem Addressed](#-problem-addressed)
- [Intended Use](#-intended-use)
- [Key Outputs](#-key-outputs)
- [Credits](#-credits)

---

## Project Overview

This dissertation investigates whether machine learning methods can reconstruct **detrended and de-seasoned ocean heat content (OHC) anomalies** under sparse observational conditions.

The work follows the Experiment A protocol of the MapEval4OceanHeat (ME4OH) initiative.

The project was conducted in collaboration with the UK Met Office as part of the MSc in Applied Data Science and Statistics at the University of Exeter.

All modelling, spatial processing, and evaluation were implemented entirely in **R**.


---

## Research Motivation

Monitoring ocean heat content is critical because:

- The ocean absorbs over 90% of excess heat from anthropogenic climate change.
- Ocean heat is a key indicator of global warming.
- Historical ocean observations are spatially sparse, especially pre-Argo.

Traditional mapping approaches rely heavily on statistical interpolation and physical assumptions. This project explores whether **data-driven models** can offer complementary reconstruction capability — particularly in data-sparse regions.

---

## Research Question

> Can machine learning methods reliably reconstruct detrended ocean heat content anomalies under realistic sparse-observation scenarios?

---

## Technologies Used

### Programming
- **R** – End-to-end modelling, spatial analysis, validation, and visualisation.

### Key Packages
- `ranger` – Efficient Random Forest implementation  
- `terra` – Raster and gridded data processing  
- `sf` – Spatial vector data handling  
- `ncdf4` / `terra` – NetCDF data manipulation  
- `ggplot2` – Visualisation  

### Why These Tools?
- **Random Forest** handles nonlinear relationships without strong parametric assumptions.
- **R ecosystem** supports reproducible spatial statistics workflows.
- **NetCDF compatibility** is essential for working with climate model outputs.

---

## Methodological Approach

### Data Preparation

- Used ME4OH “model truth” ocean heat content fields.
- Applied realistic sampling masks to simulate historical observation sparsity.
- Target variable: detrended, de-seasoned OHC anomalies.

### Two-Step Random Forest Framework

**Step 1 — Seasonal Structure (T₀ + Climatology)**  
Modelled large-scale mean and seasonal structure.

**Step 2 — Residual Mapping**  
Modelled spatial residual anomalies after removing seasonal effects.

This decomposition was chosen to:

- Reduce model complexity  
- Separate structured seasonal behaviour from anomaly dynamics  
- Diagnose where spatial artefacts emerged  

### Validation Strategy

- Compared predictions against known model truth.
- Evaluated residual maps and spatial error patterns.
- Analysed performance in sparse vs dense regions.

---

## How the Project Developed

This project was developed as part of a supervised collaboration with the UK Met Office and University of Exeter, contributing to the broader ME4OH intercomparison effort.

The goal was to:

- Diagnose how ML behaves in physically constrained geospatial systems.
- Evaluate robustness under controlled data sparsity.
- Understand limitations of purely data-driven reconstructions.

---

## Challenges & Limitations

### Data Sparsity
Sparse regions led to:

- Spatial artefacts  
- Over-smoothed reconstructions  
- Reduced reliability in poorly observed basins  

### Non-Physical Learning
Machine learning models:

- Do not inherently respect conservation laws  
- May introduce physically implausible structures  

### Validation Complexity
Because the target variable was already anomalies, care was required to avoid mis-specification or double-detrending.

---

## Problem Addressed

> How can we reliably reconstruct historical ocean heat content when direct observations are limited?

Reliable reconstructions are critical for:

- Climate monitoring  
- Attribution studies  
- Policy-relevant climate indicators  

---

## Intended Use

This repository is intended to:

- Demonstrate a reproducible ML pipeline for spatial climate reconstruction.
- Provide insight into strengths and weaknesses of Random Forest in geospatial anomaly mapping.
- Serve as a technical portfolio piece for roles in:
  - Climate analytics  
  - Environmental modelling  
  - Geospatial machine learning  

---

## Key Outputs

- Global maps of reconstructed OHC anomalies  
- Climatology and baseline structure diagnostics  
- Residual spatial error visualisations  
- NetCDF prediction outputs for scientific comparison  

---

## Credits

- UK Met Office (Supervisory Collaboration)  
- University of Exeter (MSc Programme)  
- Mapping Evaluation For Ocean Heat (ME4OH) evaluation framework  

---

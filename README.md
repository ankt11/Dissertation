# Reconstructing Ocean Heat Content Anomalies with Random Forests

MSc Applied Data Science and Statistics, University of Exeter (2025)  
Supervised collaboration with the UK Met Office

[Full project walkthrough](eclectic-phoenix-32ab7d.netlify.app)

---

## Overview

This project applies Random Forest regression to reconstruct detrended, de-seasonalised ocean heat content (OHC) anomalies from sparse in-situ observations which is a core challenge in historical climate monitoring. It was developed as part of the ME4OH (Mapping Evaluation for Ocean Heat) intercomparison initiative using synthetic profile data from the OFAM3 ocean model.

Two modelling strategies were compared: a **Baseline RF** that predicts anomalies directly, and a **Two-Step RF** that separates seasonal and residual components before reconstruction.

---

## Key Results

The Baseline RF outperformed the Two-Step variant across all metrics on the held-out test set (20% split):

| Model | RMSE | MAE | Bias | R² |
|---|---|---|---|---|
| Baseline RF | 0.463 | 0.236 | 0.001 | 0.656 |
| Two-Step RF | 0.511 | 0.271 | 0.000 | 0.583 |

Findings: Both models captured basin-scale anomaly structure, but reconstruction skill was constrained by observational sampling density. The Two-Step decomposition introduced spatial artefacts due to variance misallocation which shows when decomposition helps versus hinders in geospatial ML.

---

## Approach

Features were engineered from profile metadata: latitude, longitude, decimal year, and annual sinusoidal harmonics to represent spatio-temporal structure. The target variable, Layer 1 OHC anomalies (0–286.6 m depth), was used directly as detrended and de-seasonalised under Experiment A.

The two-step framework separated seasonal structure (T₀ + climatology) from residual anomaly mapping, with the intent of reducing model complexity at each stage. Validation used held-out test profiles plus spatial and temporal diagnostics across Pre-Argo (1993–2004) and Argo (2005–2014) eras.

Full-field reconstructions were produced over the North Atlantic (0°–60°N, 280°–360°E) on a 0.5° grid and compared against withheld ME4OH model truth fields.

---

## Repository

```
├── train_rf_models.R       # Trains Baseline and Two-Step RF models
├── rf_fullfield.R          # Generates full-field reconstructions over the North Atlantic
└── README.md
```

*Additional scripts and outputs will be added as the repository develops.*

---

## Tools

R · `ranger` · `terra` · `sf` · `ncdf4` · `ggplot2`

---

## Limitations

Reconstruction quality degrades in data-sparse regions regardless of modelling approach. Spatial error is driven more by observation density than model choice. Random Forest does not enforce physical conservation laws, which can introduce implausible structures in poorly sampled basins.

---

## Credits

Supervisors: Donata Giglio, James Salter, Matt Palmer  
UK Met Office · University of Exeter · ME4OH initiative

---


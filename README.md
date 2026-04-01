# Project Title: Optimal Field Sampling Framework for Disease Surveillance of Clonally Propagated Plants
<u> </u>
## 📖 Associated Manuscripts

This repository accompanies:

* 📄 Manuscript under review — Journal of Theoretical Biology (JTB)
* 📄 Preprint:
“Spatial dynamics and decision-theoretic surveillance in vegetatively propagated crops”
https://www.biorxiv.org/content/10.64898/2025.12.09.693344v1

## 🎯 Overview

This repository implements a mechanistic, spatially explicit Delay Differential Equation (DDE) model for pathogen dynamics in vegetatively propagated crops.

The framework integrates:

* 🌱 Within-farm epidemiology (S–L–I compartments)
* 🌍 Spatial coupling across farms (distance-based kernel)
* ⏳ Fixed latency with survival (true DDE, not ODE approximation)
* 🧠 Behavioral feedback in planting decisions
* 🔁 Seasonal resetting dynamics
* ⚖️ Exact population conservation
* 📊 Decision-theoretic optimization of surveillance effort

## 🧠 Key Scientific Contributions

### 1. Delay-driven epidemic dynamics

* Explicit latency ($\tau$) modeled via true delay differential equations
* Survival-adjusted progression using exponential decay
* Demonstrates non-trivial impact of latency on epidemic growth and detectability

### 2. Spatial transmission structure

* Farms interact via a distance-decay kernel:
  
  $$𝑊_{ij} = 𝑒^{- \gamma D_{ij}}$$
​
* Captures localized clustering and long-range spread

### 3. Behavioral recruitment

* Nonlinear farmer response:
  
$$f(\psi) = \frac{\psi}{ (1+ \kappa \psi)}$$

* Models adaptive planting under perceived infection risk

### 4. Analytical reproduction number

The model derives a latency-dependent reproduction number:

 $$R_0(\tau) = \frac{\beta}{\gamma} (1+e^{-\gamma \tau}) + e^{-\gamma \tau}$$

### 5. Decision-theoretic surveillance

A full Bayesian decision framework is implemented to compute:

* Optimal sampling effort $𝑛$
* Expected Net Gain of Sampling (ENGS)
* Trade-offs between:
* Epidemiological risk
* Economic cost
* Detection uncertainty

## 📂 Repository Structure

```bash
.
├── main_model.R            # Core DDE simulation
├── /data/                  # Input or generated datasets
├── /outputs/               # Simulation outputs
├── /figures/               # Figures for manuscript
├── /spatial/               # GIS-related scripts
├── /sensitivity/           # Global sensitivity analysis
├── /optimization/          # Decision-theoretic sampling
└── README.md
```

### ▶️ Running the Code

#### 1️⃣ Install dependencies

```R
install.packages(c(
  "deSolve",
  "tidyverse",
  "reshape2",
  "lhs",
  "sensitivity",
  "sf",
  "gstat",
  "geodata",
  "viridis",
  "patchwork",
  "ggspatial",
  "extraDistr"
))
```
#### 2️⃣ Run simulation

```R
source("main_model.R")
```

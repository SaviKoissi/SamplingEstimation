** Project Title: Optimal Field Sampling Framework for Disease Surveillance of Clonally Propagated Plants**

📖 Associated Manuscripts

This repository accompanies:

* 📄 Manuscript under review — Journal of Theoretical Biology (JTB)
* 📄 Preprint:
“Spatial dynamics and decision-theoretic surveillance in vegetatively propagated crops”
https://www.biorxiv.org/content/10.64898/2025.12.09.693344v1

🎯 Overview

This repository implements a mechanistic, spatially explicit Delay Differential Equation (DDE) model for pathogen dynamics in vegetatively propagated crops.

The framework integrates:

* 🌱 Within-farm epidemiology (S–L–I compartments)
* 🌍 Spatial coupling across farms (distance-based kernel)
* ⏳ Fixed latency with survival (true DDE, not ODE approximation)
* 🧠 Behavioral feedback in planting decisions
* 🔁 Seasonal resetting dynamics
* ⚖️ Exact population conservation
* 📊 Decision-theoretic optimization of surveillance effort

🧠 Key Scientific Contributions
1. Delay-driven epidemic dynamics

* Explicit latency ($\tau$) modeled via true delay differential equations
* Survival-adjusted progression using exponential decay
* Demonstrates non-trivial impact of latency on epidemic growth and detectability

2. Spatial transmission structure
* Farms interact via a distance-decay kernel:
  
$$𝑊_{ij} = 𝑒^{- \gamma D_{ij}}$$
​
* Captures localized clustering and long-range spread

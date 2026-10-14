# R₀ Finder

**Automated derivation of the basic reproduction number (R₀) from ODE compartmental models and Petri nets.**

R₀ Finder is a desktop application that

- computes R₀ symbolically from the topology of a **Stochastic Petri Net (SPN)** or a **Variable Arc Weight Petri Net (VAPN)** using the next-generation-matrix method;
- maps a system of **ordinary differential equations** into an SPN or a VAPN while preserving the biological structure of the model (compartments, flows, catalysts, conservation), so that the R₀ obtained from the net equals the R₀ of the ODE system;
- computes R₀ directly from the ODE system, resolves the disease-free equilibrium automatically (including rank-deficient systems), and reports normalized sensitivity (elasticity) indices of R₀ with respect to every parameter.

The methodology is described in

> Reckell, T., Sterner, B., & Jevtić, P. (2025). The basic reproduction number for Petri net models: A next-generation matrix approach. *Applied Sciences*, 15(23), 12827. https://doi.org/10.3390/app152312827

and in a second manuscript describing the software (in review, *PeerJ Computer Science*). The software's R₀ construction generalizes the published one; see the paper for details.

## Requirements

- Python 3.9 or later
- `sympy`, `matplotlib`, `pillow`
- `tkinter` (ships with most Python installations; on some Linux distributions install `python3-tk`)

```bash
pip install sympy matplotlib pillow
```

## Running

```bash
python R0_finder.py
```

From the home screen choose **ODE**, **SPN** or **VAPN**.

- **ODE**: enter one equation per state variable (or load a saved `.json` model), then **Calculate R₀**, or **Map to SPN** / **Map to VAPN** to generate a Petri net in a new tab.
- **SPN / VAPN**: draw places, transitions and arcs, edit propensity functions or arc weights, then **Calculate R₀**.
- In the calculation window: tick the infected compartments, optionally mark free variables or add constraints, **Auto-Detect DFE**, then **Calculate R₀**. Answer any prompt about an ambiguous transition (new infection vs. migration) or an ambiguous proportional flow (merge vs. separate). **Sensitivity for R₀** opens the elasticity analysis.

Example models (SIR, SIRS, SEIR, SEEIR, SVEIR, COVID-19, vector-borne, two-patch and others) are included as `.json` files.

## Files

| File | Contents |
|---|---|
| `R0_logic.py` | Core engine: Petri net classes, ODE parser and mapper, DFE solver, R₀ and sensitivity computation |
| `R0_finder.py` | Graphical interface |
| `models/` | Example ODE, SPN and VAPN models |

## Citation

If you use R₀ Finder, please cite the *Applied Sciences* paper above and the software release (Zenodo DOI to be added on release).

## License

Copyright (c) 2026 [Trevor Reckell].

R₀ Finder is released under the **PolyForm Noncommercial License 1.0.0** (see `LICENSE`). You may use, modify and redistribute it free of charge for **non-commercial purposes**, including academic research, teaching, and use by universities, non-profit organizations and government agencies. **Any commercial use requires a separate license.** For commercial licensing, contact Skysong Innovations (Arizona State University): https://skysonginnovations.com.

**Patent pending.** The methods implemented in this software are the subject of a pending patent application.

## Funding

Development was supported in part by NIH grants 5R01GM131405-02 and DMS-1615879.

## Authors

Trevor Reckell, Beckett Sterner, Petar Jevtić — Arizona State University.

# BCQM-Programs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17191306.svg)](https://doi.org/10.5281/zenodo.17191306)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17242311.svg)](https://doi.org/10.5281/zenodo.17242311) (Analytical Proof)

Analytical and computational companions for the BCQM project.

This repository now includes **`bcqm_tools/`**, a *drop‑in, reproducible, and supported* toolset that
supersedes the original ad‑hoc scripts in `double_slit/` and `toy_experiment/`. Those legacy folders
remain for transparency, but new users should start with **`bcqm_tools`**.

---

## 📦 What is `bcqm_tools/`?
A self-contained folder you can run **without installing anything** (or optionally install locally).
It provides:

- A **CLI** with three programs:
  - `ramsey` – computes D(t)=exp(-∫_0^t Γ) , the bound F_opt ≤ 1/2 (1+D), and the **W** threshold.
  - `gkls` – minimal, phase‑covariant GKLS qubit simulator (dephasing ± relaxation); plots Bloch components.
  - `galton` – seeded Galton channels toy model (repeatable histogram).
- **YAML configs** (zero typing of parameters).
- **Deterministic outputs** (seeded RNG) + **metadata JSON** per run.
- **Tests** (`pytest`) for quick smoke checks.

> For full details, see `bcqm_tools/README.md`.

---

## 🚀 Quickstart (repo‑root commands)

```bash
# 1) Create & activate a local virtual environment
python -m venv .venv
# macOS/Linux
source .venv/bin/activate
# Windows (PowerShell)
# .venv\Scripts\Activate.ps1

# 2) Install minimal dependencies for the tools
pip install -r bcqm_tools/requirements.txt

# 3) Run examples (no install required)
python -m bcqm_tools.cli --help

# Ramsey (constant gamma): writes CSV/PNG/JSON under bcqm_tools/examples/
python -m bcqm_tools.cli ramsey --gamma 0.1 --fstar 0.9 --tmax 5e-5 --dt 2e-7 --out bcqm_tools/examples/ramsey

# GKLS qubit (phase‑covariant)
python -m bcqm_tools.cli gkls --gamma-phi 0.1 --gamma-relax 0.0 --tmax 5e-5 --dt 5e-7 --out bcqm_tools/examples/gkls

# Galton channels (seeded)
python -m bcqm_tools.cli galton --rows 10 --n 1000 --seed 42 --out bcqm_tools/examples/galton

# YAML-driven (zero typing)
python -m bcqm_tools.run_from_yaml --config bcqm_tools/configs/ramsey_constant.yml
python -m bcqm_tools.run_from_yaml --config bcqm_tools/configs/gkls_basic.yml
python -m bcqm_tools.run_from_yaml --config bcqm_tools/configs/galton.yml

# (Optional) Smoke tests
pytest -q bcqm_tools/tests
```

---

## 🗂 Repository structure

```
BCQM-Programs/
├─ bcqm_tools/                 # NEW: reproducible, supported tools (CLI + YAML + tests + docs)
│  ├─ README.md
│  ├─ requirements.txt
│  ├─ configs/                 # ready-to-run YAMLs
│  ├─ examples/                # outputs written here (CSV/PNG/JSON)
│  ├─ tests/                   # pytest smoke tests
│  ├─ __init__.py  cli.py  run_from_yaml.py  ramsey.py  gkls.py  galton.py  utils.py
│  └─ pyproject.toml  LICENSE  CITATION.cff  .gitignore
│
├─ double_slit/                # legacy scripts (kept for transparency)
├─ toy_experiment/             # legacy scripts (kept for transparency)
│
├─ CITATION.cff
├─ README.md                   # (this file)
├─ requirements.txt            # (optional; repo-wide requirements if used)
├─ .gitignore
├─ .gitattributes
└─ .DS_Store
```

**Why this structure?** So a new reader can clone the repo, run the tools immediately from the root,
and reproduce the key figures/numerics without touching legacy code.

---

## 📑 Citing

Please see the repository’s `CITATION.cff`. If you archive the repository or releases to Zenodo, add the DOI here.

---

## 🔁 Legacy vs. `bcqm_tools`

- The **legacy** folders (`double_slit/`, `toy_experiment/`) were minimal, exploratory scripts.  
- **`bcqm_tools/`** is the production path: CLI, configs, seeded runs, metadata, and tests.
  If you are reviewing or reproducing the BCQM results, start there.

---

## (Original README content)

# BCQM Programs

The code in this repository is the actual code used in investigating the two time axes.
it has not been changed.

Simulation code supporting the **Boundary-Condition Quantum Mechanics (BCQM)** framework.  
This repository contains two core models:

- **Toy Experiment** (`toy_experiment/`):  
  A minimal BCQM simulator demonstrating collapse dynamics, horizon effects, and ensemble reproducibility.  
- **Double-Slit Simulation** (`double_slit/`):  
  Numerical model of interference, decoherence, and collapse horizons in the double-slit setup.

---

## Repository Structure

```
bcqm-programs/
├─ toy_experiment/
│  ├─ bcqm_toy_experiment.py
│  ├─ *.csv                  # example output data
│  └─ README.md              # details for this module
├─ double_slit/
│  ├─ bcqm_double_slit.py
│  ├─ bcqm_double_slit_sim.py
│  └─ README.md              # details for this module
├─ requirements.txt
├─ .gitignore
└─ README.md                 # (this file)
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<your-username>/bcqm-programs.git
cd bcqm-programs
```

Create a virtual environment and install dependencies:

```bash
python3 -m venv venv
source venv/bin/activate   # or `venv\Scripts\activate` on Windows
pip install -r requirements.txt
```

---

## Usage

### Toy Experiment
```bash
cd toy_experiment
python bcqm_toy_experiment.py
```

### Double-Slit Simulation
```bash
cd double_slit
python bcqm_double_slit_sim.py
```

Both programs generate data/figures into their respective directories.

---

## Requirements

- Python ≥3.8  
- Dependencies: `numpy`, `scipy`, `matplotlib`, `pandas`

---

## License

MIT License (or another license if you prefer).  
You are free to use, modify, and distribute this code with attribution.

---

## Citation

If you use this code in academic work, please cite:

> Ferguson, P. M. (2025). *Boundary-Condition Quantum Mechanics (BCQM): Simulation Programs*. GitHub.  
> https://github.com/<your-username>/bcqm-programs
---

# BCQM Programs

This repository contains the two small simulation programs used alongside the BCQM preprint.

- **double_slit_sim/** – double-slit ensemble simulation producing interference/incoherent outputs.
- **bcqm_toy_experiment/** – toy collapse dynamics producing damping/coherence figures.

Each subfolder includes:
- A short README explaining how the figures were generated,
- The exact `outputs/*.png` used in the paper,
- A `PROVENANCE.txt` noting environment and date.

> These scripts retain “fast/stub” branches used during exploratory runs. They’re documented in each README.

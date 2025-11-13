# bcqm_tools — drop-in, runnable tools (no install required)

Copy this **bcqm_tools/** folder into your repository root (e.g. `BCQM-Programs/`) and run the commands below
from the repo root. These tools reproduce the Ramsey/W bound, a phase-covariant GKLS qubit, and a seeded
Galton toy model. Outputs (CSV/PNG/JSON) are written under `bcqm_tools/examples/` with run metadata.

## 0) One-time setup
```bash
python -m venv .venv
# macOS/Linux
source .venv/bin/activate
# Windows (PowerShell)
# .venv\Scripts\Activate.ps1

pip install -r bcqm_tools/requirements.txt
```

## 1) Run (no install)

```bash
# Show help
python -m bcqm_tools.cli --help

# Ramsey (constant gamma) — writes CSV/PNG/JSON under bcqm_tools/examples/
python -m bcqm_tools.cli ramsey --gamma 0.1 --fstar 0.9 --tmax 5e-5 --dt 2e-7 --out bcqm_tools/examples/ramsey

# GKLS qubit (phase-covariant; pure dephasing example)
python -m bcqm_tools.cli gkls --gamma-phi 0.1 --gamma-relax 0.0 --tmax 5e-5 --dt 5e-7 --out bcqm_tools/examples/gkls

# Galton channels (seeded)
python -m bcqm_tools.cli galton --rows 10 --n 1000 --seed 42 --out bcqm_tools/examples/galton
```

## 2) YAML-driven (no typing of params)
```bash
python -m bcqm_tools.run_from_yaml --config bcqm_tools/configs/ramsey_constant.yml
python -m bcqm_tools.run_from_yaml --config bcqm_tools/configs/gkls_basic.yml
python -m bcqm_tools.run_from_yaml --config bcqm_tools/configs/galton.yml
```

## 3) Quick smoke tests
```bash
pytest -q bcqm_tools/tests
```

## 4) Optional local install (to get a `bcqm` command)
```bash
cd bcqm_tools
pip install -e .
bcqm --help
```

## Reproducibility
- All stochastic steps are **seeded**.
- Every run writes a `*.json` **metadata** file (parameters + Python/platform versions).
- Tests cover constant-γ Ramsey and pure-dephasing decay.


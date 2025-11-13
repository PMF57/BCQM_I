# BCQM Toy Sandbox — Reproducible Methods

This toy model provides a minimal, reproducible demonstration of **Boundary-Condition Quantum Mechanics (BCQM)** principles.  
It shows how one framework can simultaneously capture:

1. **Resonant eigenstate selection** (via integrated weights).  
2. **Apparent randomness of detection positions** (via fine "pin" structure inside a smooth Gaussian envelope).

---

## Model Setup

- **Space:** 1D grid, x ∈ [−L/2, L/2], with L = 60.0 and N = 256 points.  
- **Envelope:** Gaussian \( G(x) = \exp(-x^2 / (2\sigma^2)) \), with σ = 10.0.  
- **Outcomes:** M = 3, equally prepared.  
- **Sub-waves per outcome:** 3 plane waves each.  
  - Wave numbers: \( k_{i,n} = k0_i + n \, dk \), with \( dk = 0.3 \), \( k0s = [0.6, 0.9, 1.2] \).  
  - Random phases drawn per run.  
- **Propagation:** Free one-step via \( e^{-i k^2 T / (2m)} \) with T = 0.8, \( m = 1.0 \).  

---

## Coherent vs Dephased Densities

- **Coherent density:** \( P_\text{total}(x) = \sum_i |\psi_i(x)|^2 \) including sub-wave interference.  
- **Dephased baseline:** incoherent sum of squared basis modes (pins removed).  

---

## Measurement Per Run

1. **Resonant outcome selection:**  
   - Compute weights \( w_i = \int |\psi_i(x)|^2 dx \).  
   - Choose outcome *i* with probability \( w_i \).  
2. **Detection position:**  
   - Draw hit \( x_\text{hit} \) from conditional density \( |\psi_i(x)|^2 \).  

---

## Monte Carlo Simulation

- **Runs:** R = 2000.  
- **Random seeds:**  
  - 314159 for single-run structures.  
  - 271828 for Monte Carlo loop.  
- **Outputs:**  
  - Outcome frequencies.  
  - Detection positions.  
  - Densities (coherent vs dephased).  

---

## Files Produced

- `bcqm_toy_densities.csv` — coherent and dephased densities.  
- `bcqm_outcome_counts.csv` — histogram of selected outcomes.  
- `bcqm_positions.csv` — sampled detection positions.  

### Plots
- `plot_A_conditional_and_total.png` — conditional densities & total.  
- `plot_B_coherent_vs_dephased.png` — coherent vs dephased comparison.  
- `plot_C_outcome_histogram.png` — outcome histogram across runs.  
- `plot_D_positions_vs_profiles.png` — detection positions vs profiles.  

---

## How to Run

### Requirements
- Python ≥ 3.8  
- Packages: `numpy`, `matplotlib`, `pandas`  

Install dependencies (if needed):  
```bash
pip install numpy matplotlib pandas
```

### Running
1. Navigate to the `BCQM_toy_experiment` folder.  
2. Run the main script:  
   ```bash
   python bcqm_toy.py
   ```
3. Outputs (CSV + plots) will be saved in the same directory.  

---

## Interpretation

- The **Gaussian envelope** + **sub-waves** = amorphous quantum state.  
- **Pins** arise from sub-wave interference.  
- **Resonant selection** → which eigenstate outcome is realised.  
- **Conditional sampling** → where the particle appears.  
- **Dephasing** → removes pins but preserves envelope, mirroring collapse removing \( t^- \) coherence while preserving \( t^+ \) populations.  

---

## Notes

- Fully self-contained, uses only standard quantum ingredients.  
- Modifying seeds or wave numbers alters pin structure but not the overall logic.  
- Serves as a lightweight entry point for exploring BCQM ideas.  

---

## Citation

If you use this code or data, please cite:

- Ferguson, P.M. (2025). *Boundary-Condition Quantum Mechanics (BCQM) Toy Sandbox*. Zenodo. https://doi.org/10.5281/zenodo.17191307  
- GitHub Repository: https://github.com/PMF57/BCQM-Programs

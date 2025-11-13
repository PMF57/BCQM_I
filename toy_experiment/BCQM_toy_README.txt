
BCQM Toy Sandbox — Reproducible Methods (Sub-waves, Pins, Resonant Selection, Random Appearance)
===============================================================================================

Purpose
-------
Demonstrate in a minimal, reproducible way how a single model can account for:
(1) State (eigenvalue) selection by resonant weighting; and
(2) Apparent randomness of detection position via fine "pin" structure inside a smooth envelope.

Amorphousness
-------------
The particle is not localized initially. It is represented by a Gaussian envelope G(x) with width sigma = 10.0.
Each internal outcome's wavefunction is a coherent superposition of a few plane-wave sub-components under this envelope.
This captures the "amorphousness": a delocalized state with multiple channels present simultaneously.

Model Summary
-------------
Space: 1D grid on x in [-L/2, L/2], L = 60.0, N = 256.
Envelope: G(x) = exp(-x^2 / (2 sigma^2)), sigma = 10.0.

Internal outcomes: M = 3, equal preparation weights.

Sub-waves per outcome: 3 plane-wave components with wave numbers
k_i,n = k0_i + n * dk for n in {-2..1},
with dk = 0.3 and central k0s = [0.6, 0.9, 1.2].
Random phases are drawn independently per run for each component.

Propagation: one-step free propagation via k-space phase exp(-i k^2 T / (2 m_eff)), with
T = 0.8, m_eff = 1.0.

Coherent vs Dephased:
- Coherent total density P_total(x): sum over outcomes of |psi_i(x)|^2 with sub-wave interference present.
- Dephased baseline P_deph(x): sum of |basis_mode(x)|^2 over all sub-waves (removes sub-wave cross terms), normalized.

Measurement (per run)
---------------------
1) Resonant outcome selection:
   Compute weights w_i = ∫ |psi_i(x)|^2 dx for i=0..M-1 (normalized to 1). Choose outcome i with probability w_i.
2) Position appearance:
   Given outcome i, draw detection position x_hit from the conditional density |psi_i(x)|^2.

Monte Carlo
-----------
R = 2000 runs. Independent random phases each run (seeded with 271828 for reproducibility), outcome frequencies and positions recorded.

Files Produced
--------------
- /mnt/data/bcqm_toy_densities.csv
  Columns: x, P_total_coherent, P_dephased_baseline.

- /mnt/data/bcqm_outcome_counts.csv
  Columns: outcome_index, count. (Histogram of chosen outcomes across R runs.)

- /mnt/data/bcqm_positions.csv
  Column: x_hit. (Detected positions across R runs.)

Plots
-----
- /mnt/data/plot_A_conditional_and_total.png: Single-run conditional densities per outcome, plus coherent total.
- /mnt/data/plot_B_coherent_vs_dephased.png: Coherent total vs dephased baseline (pins removed).
- /mnt/data/plot_C_outcome_histogram.png: Outcome histogram across R runs.
- /mnt/data/plot_D_positions_vs_profiles.png: Detected positions vs coherent & dephased profiles.

Reproducibility Notes
---------------------
- Random number generators: main seed = 314159 for single-run structures; 271828 for Monte Carlo loop.
- All parameters are specified above; changing seeds, modes_per_state, dk, or k0s alters pin structure but not the logic.
- This toy is intentionally minimal and uses only standard linear quantum ingredients; no extra variables beyond the wavefunction's sub-components are introduced.

Interpretation
--------------
- The Gaussian envelope + superposed sub-waves is the "amorphous" quantum state.
- Sub-wave interference creates fine "pins" inside the envelope.
- Resonant selection (integrated conditional weights) models which eigenstate the apparatus reports.
- Conditional sampling from |psi_i|^2 models where the particle "appears."
- Dephasing removes the pins but not the envelope, mirroring collapse erasing coherence (t^-) but leaving populations (t^+).

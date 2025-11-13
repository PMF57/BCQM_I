#!/usr/bin/env python3
"""
BCQM Toy Sandbox — Sub-waves ("pins"), Resonant State Selection, Random Appearance
==================================================================================

Purpose
-------
Demonstrate, in a minimal and reproducible way, how a single quantum-wave toy model can account for:
1) State (eigenvalue) selection via resonant weighting over internal outcomes, and
2) Apparent randomness of detection position via fine "pin" structure (coherent sub-wave interference)
   inside a smooth envelope.
"""

import os
import numpy as np
import pandas as pd

# Optional plotting (set to True to produce PNGs)
ENABLE_PLOTS = False
if ENABLE_PLOTS:
    import matplotlib.pyplot as plt

# -------------------------------
# Parameters
# -------------------------------

# Space and envelope (amorphousness)
L = 60.0
N = 256
sigma = 10.0

# Internal outcomes (eigenmodes)
M = 3
prep_weights = np.array([1.0, 1.0, 1.0], dtype=float)

# Sub-waves per outcome
modes_per_state = 3
dk = 0.30
k0s = np.array([0.6, 0.9, 1.2])

# Propagation parameters
T = 0.8
m_eff = 1.0

# Monte Carlo runs
R = 2000

# Seeds for reproducibility
SEED_SINGLE_RUN = 314159
SEED_MONTE_CARLO = 271828

# Output directory
OUT_DIR = os.path.join(os.path.dirname(__file__), "outputs")


def gaussian_envelope(x, sigma):
    return np.exp(-0.5 * (x/sigma)**2)


def propagate_free(field_x, x, T, m_eff):
    """Free propagation for time T in 1D using FFT."""
    N = x.size
    dx = x[1] - x[0]
    dk_fft = 2*np.pi/(N*dx)
    k_fft = dk_fft * (np.arange(N) - N//2)
    Psi_k = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(field_x)))
    phase = np.exp(-1j * (k_fft**2) * T / (2*m_eff))
    Psi_k *= phase
    out = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Psi_k)))
    return out


def normalized_total(psi_components, x):
    total_intensity = np.trapz(np.sum(np.abs(psi_components)**2, axis=0), x)
    return psi_components / np.sqrt(total_intensity)


def sample_from_pdf(x, pdf, rng):
    w = pdf / np.trapz(pdf, x)
    w_disc = w / w.sum()
    idx = rng.choice(len(x), p=w_disc)
    return x[idx]


def build_basis_modes(x, sigma, T, m_eff, k0s, dk, modes_per_state):
    N = x.size
    M = k0s.size
    k_indices = np.arange(-(modes_per_state//2), modes_per_state//2 + 1)
    G = gaussian_envelope(x, sigma)
    basis_modes = np.zeros((M, modes_per_state, N), dtype=complex)
    for i in range(M):
        ks = k0s[i] + dk * k_indices
        for n, k_val in enumerate(ks):
            mode_field = G * np.exp(1j * k_val * x)
            basis_modes[i, n] = propagate_free(mode_field, x, T, m_eff)
    return basis_modes


def dephased_baseline(basis_modes, x):
    P = np.sum(np.abs(basis_modes)**2, axis=(0, 1))
    return P / np.trapz(P, x)


def build_coherent_components_single_run(basis_modes, prep_weights, rng, x):
    M, modes_per_state, N = basis_modes.shape
    psi_components = np.zeros((M, N), dtype=complex)
    for i in range(M):
        phases = rng.uniform(0, 2*np.pi, size=modes_per_state)
        coeffs = np.exp(1j * phases)
        field = np.tensordot(coeffs, basis_modes[i], axes=(0, 0))
        psi_components[i] = prep_weights[i] * field
    psi_components = normalized_total(psi_components, x)
    return psi_components


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    x = np.linspace(-L/2, L/2, N)
    basis = build_basis_modes(x, sigma, T, m_eff, k0s, dk, modes_per_state)
    P_deph = dephased_baseline(basis, x)

    # Single-run snapshot
    rng_single = np.random.default_rng(SEED_SINGLE_RUN)
    psi_run = build_coherent_components_single_run(basis, prep_weights, rng_single, x)
    P_i = np.abs(psi_run)**2
    P_total = np.sum(P_i, axis=0)
    P_total /= np.trapz(P_total, x)

    densities_path = os.path.join(OUT_DIR, "bcqm_toy_densities.csv")
    pd.DataFrame({"x": x, "P_total_coherent": P_total, "P_dephased_baseline": P_deph}).to_csv(densities_path, index=False)

    # Monte Carlo
    outcome_counts = np.zeros(M, dtype=int)
    positions = []
    rng_mc = np.random.default_rng(SEED_MONTE_CARLO)

    for _ in range(R):
        psi_r = build_coherent_components_single_run(basis, prep_weights, rng_mc, x)
        P_i_r = np.abs(psi_r)**2
        w_i = np.trapz(P_i_r, x, axis=1)
        w_i /= w_i.sum()
        i_out = rng_mc.choice(M, p=w_i)
        outcome_counts[i_out] += 1
        x_hit = sample_from_pdf(x, P_i_r[i_out], rng_mc)
        positions.append(x_hit)

    outcomes_path = os.path.join(OUT_DIR, "bcqm_outcome_counts.csv")
    positions_path = os.path.join(OUT_DIR, "bcqm_positions.csv")

    pd.DataFrame({"outcome_index": np.arange(M), "count": outcome_counts}).to_csv(outcomes_path, index=False)
    pd.DataFrame({"x_hit": positions}).to_csv(positions_path, index=False)

    print("Done. Outputs in:", OUT_DIR)


if __name__ == "__main__":
    main()

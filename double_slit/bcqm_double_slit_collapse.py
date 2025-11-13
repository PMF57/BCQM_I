#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BCQM Double-Slit with Selective Damping (Collapse = Absorption of a Dominant Wave)
-------------------------------------------------------------------------------
- Propagates a 1D double-slit q-wave using the angular spectrum method.
- Decomposes the field into symmetric (S) and antisymmetric (AS) parts.
- At z ≈ z_damp, applies *local* damping only to the *dominant* component (S or AS).
- Continues propagation; saves diagnostics and the final screen profile.

Outputs (saved under ./out/):
- snapshots_profiles.png
- dominance_map.png
- damping_diagnostics.png
- final_screen_profile.png

Dependencies: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft, fftfreq
import os

# --------------------------- Propagation Core ---------------------------------
def angular_spectrum_step(psi, x, wavelength, dz):
    """
    One paraxial propagation step dz via angular spectrum method (1D).
    """
    k = 2*np.pi / wavelength
    dx = x[1] - x[0]
    N = x.size
    fx = fftfreq(N, d=dx)
    kx = 2*np.pi*fx
    kz = np.sqrt(np.maximum(0.0, k**2 - kx**2))  # evanescent filtered
    Psi = fft(psi)
    Psi *= np.exp(1j * kz * dz)
    return ifft(Psi)

def symmetric_antisymmetric(psi):
    """
    Decompose psi(x) into symmetric (S) and antisymmetric (AS) parts assuming
    x is a symmetric, uniformly spaced grid (x[i] = -x[-i-1]).
    """
    psi_rev = psi[::-1]
    psi_S  = 0.5 * (psi + psi_rev)
    psi_AS = 0.5 * (psi - psi_rev)
    return psi_S, psi_AS

# ----------------------- BCQM Collapse (Selective Damping) --------------------
def apply_selective_damping(psi, alpha):
    """
    Apply spatially selective damping to the *dominant* component at each x.
    - Split into symmetric/antisymmetric parts: psi_S, psi_AS
    - Build boolean masks where each is locally dominant
    - Dampen the dominant part by factor (1 - alpha) with 0 < alpha < 1

    Returns:
        psi_damped : recombined field after damping
        mask_S     : boolean array where S dominates
        mask_AS    : boolean array where AS dominates
        psi_S_pre, psi_AS_pre : components before damping (for diagnostics)
    """
    psi_S, psi_AS = symmetric_antisymmetric(psi)
    IS  = np.abs(psi_S)**2
    IAS = np.abs(psi_AS)**2
    mask_AS = IAS > IS
    mask_S  = ~mask_AS

    psi_S_d  = psi_S * (1 - alpha * mask_S)
    psi_AS_d = psi_AS * (1 - alpha * mask_AS)
    psi_d    = psi_S_d + psi_AS_d
    return psi_d, mask_S, mask_AS, psi_S, psi_AS

# --------------------------- Initial Conditions -------------------------------
def initialize_double_slit(x, slit_sep=0.40e-3, slit_w=0.08e-3):
    """
    Two Gaussian slits centered at ±slit_sep/2 with width slit_w.
    Normalized to unit total intensity.
    """
    x0 = slit_sep/2
    slit1 = np.exp(-0.5*((x - x0)/slit_w)**2)
    slit2 = np.exp(-0.5*((x + x0)/slit_w)**2)
    psi0 = slit1 + slit2
    # Normalize
    norm = np.trapz(np.abs(psi0)**2, x)
    if norm > 0:
        psi0 /= np.sqrt(norm)
    return psi0

# ------------------------------- Simulation -----------------------------------
def run_simulation(
    wavelength=532e-9,   # 532 nm (green)
    x_span=1.5e-3,       # half-width of computational window (meters)
    Nx=2048,             # spatial samples
    z_total=1.2,         # total propagation distance (meters)
    Nz=600,              # number of z-steps
    slit_sep=0.40e-3,    # slit separation (meters)
    slit_w=0.08e-3,      # slit Gaussian width (meters)
    z_damp=0.60,         # where the damping (collapse) is applied (meters)
    alpha=0.35           # damping strength (0<alpha<1)
):
    """
    Simulate q-wave through a double slit and apply BCQM-style collapse by
    selectively damping the locally dominant component at z ≈ z_damp.
    """
    os.makedirs("out", exist_ok=True)

    # Spatial and longitudinal grids
    x = np.linspace(-x_span, x_span, Nx)
    z_grid = np.linspace(0.0, z_total, Nz)
    dz = z_grid[1] - z_grid[0]

    # Initial field just after the mask
    psi = initialize_double_slit(x, slit_sep, slit_w)

    # Storage for diagnostics
    snapshots = []        # tuples (z, psi, psi_S, psi_AS)
    dominance_maps = []   # tuples (z, dom_mask) with dom_mask in {0(S),1(AS)}
    damp_applied = False
    diag = None

    for zi, z in enumerate(z_grid):
        # Save a few longitudinal snapshots
        if zi % max(1, Nz//6) == 0:
            psi_S, psi_AS = symmetric_antisymmetric(psi)
            snapshots.append((z, psi.copy(), psi_S.copy(), psi_AS.copy()))
            dom_mask = (np.abs(psi_AS)**2 > np.abs(psi_S)**2).astype(float)
            dominance_maps.append((z, dom_mask))

        # Apply damping once, at z ≈ z_damp
        if (not damp_applied) and (z >= z_damp):
            psi, mask_S, mask_AS, psi_S_pre, psi_AS_pre = apply_selective_damping(psi, alpha)
            damp_applied = True
            diag = dict(z=z, mask_S=mask_S, mask_AS=mask_AS,
                        psi_S_pre=psi_S_pre, psi_AS_pre=psi_AS_pre,
                        psi_post=psi.copy())

        # Propagate forward
        if zi < Nz - 1:
            psi = angular_spectrum_step(psi, x, wavelength, dz)

    # Final screen intensity (normalized)
    I_final = np.abs(psi)**2
    I_final /= np.trapz(I_final, x)

    # ------------------------------- Plots ------------------------------------
    # 1) Snapshots along z
    nrows = len(snapshots)
    fig, axs = plt.subplots(nrows, 1, figsize=(9, 10), sharex=True)
    if nrows == 1:
        axs = [axs]
    for ax, (z, psi_z, psiS_z, psiAS_z) in zip(axs, snapshots):
        ax.plot(x*1e3, np.abs(psi_z)**2, label=f"|ψ|² @ z={z:.2f} m")
        ax.plot(x*1e3, np.abs(psiS_z)**2, linestyle='--', label="|ψ_S|²")
        ax.plot(x*1e3, np.abs(psiAS_z)**2, linestyle=':', label="|ψ_AS|²")
        ax.legend(fontsize=8, ncol=3)
        ax.set_ylabel("Intensity")
    axs[-1].set_xlabel("x (mm)")
    fig.suptitle("Propagation snapshots (before/after damping)")
    fig.tight_layout()
    fig.savefig("out/snapshots_profiles.png", dpi=200)
    plt.close(fig)

    # 2) Dominance map
    Zs = [z for z,_ in dominance_maps]
    D  = np.stack([d for _,d in dominance_maps], axis=0)
    fig, ax = plt.subplots(figsize=(8,4))
    im = ax.imshow(D, aspect="auto",
                   extent=[x[0]*1e3, x[-1]*1e3, Zs[-1], Zs[0]])
    ax.set_title("Dominance map (1=AS dominant, 0=S dominant)")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("z (m)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig("out/dominance_map.png", dpi=200)
    plt.close(fig)

    # 3) Damping diagnostics at the collapse plane
    if damp_applied and diag is not None:
        fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
        axs[0].plot(x*1e3, np.abs(diag["psi_S_pre"])**2, label="|ψ_S|² (pre)")
        axs[0].plot(x*1e3, np.abs(diag["psi_AS_pre"])**2, label="|ψ_AS|² (pre)")
        axs[0].legend()
        axs[0].set_title(f"Components at damping plane z≈{diag['z']:.3f} m")

        axs[1].plot(x*1e3, diag["mask_S"].astype(int), label="S dominant mask")
        axs[1].plot(x*1e3, diag["mask_AS"].astype(int), label="AS dominant mask")
        axs[1].legend()

        axs[2].plot(x*1e3, np.abs(diag["psi_post"])**2, label="|ψ|² after damping")
        axs[2].legend()
        axs[2].set_xlabel("x (mm)")
        fig.tight_layout()
        fig.savefig("out/damping_diagnostics.png", dpi=200)
        plt.close(fig)

    # 4) Final screen profile
    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(x*1e3, I_final)
    ax.set_title(f"Final screen intensity at z={z_total} m")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("Normalized intensity")
    fig.tight_layout()
    fig.savefig("out/final_screen_profile.png", dpi=200)
    plt.close(fig)

    print("Done. Figures saved to ./out")

    # Return useful objects for programmatic use
    return {
        "x": x, "I_final": I_final, "snapshots": snapshots,
        "dominance_maps": dominance_maps, "damp_applied": damp_applied
    }

# ------------------------------- CLI ------------------------------------------
if __name__ == "__main__":
    # Default run; adjust parameters by editing below if needed.
    run_simulation()

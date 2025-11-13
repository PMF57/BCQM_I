#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft, fftfreq
from scipy.signal import find_peaks
from mpl_toolkits.mplot3d import Axes3D  # noqa

def angular_spectrum_propagate(psi0, x, wavelength, z):
    k = 2*np.pi / wavelength
    dx = x[1] - x[0]
    N = x.size
    fx = fftfreq(N, d=dx)
    kx = 2*np.pi*fx
    kz_sq = k**2 - kx**2
    kz = np.sqrt(np.clip(kz_sq, 0.0, None))
    H = np.exp(1j * kz * z)
    return ifft(fft(psi0) * H)

def gaussian_beam_envelope(x, w0):
    return np.exp(-(x**2)/(w0**2))

def rectangular_slit(x, center, width):
    return ((x >= center - width/2) & (x <= center + width/2)).astype(float)

def main():
    # --- Parameters (tuned for ~7 fringes) ---
    wavelength = 633e-9     # meters
    slit_sep   = 120e-6
    slit_width = 20e-6
    beam_waist = 0.6e-3
    x_span     = 2.0e-3
    nx         = 2048
    z_screen   = 1.4        # meters
    nz         = 220

    # Grid
    x = np.linspace(-x_span, x_span, nx)
    z_slices = np.linspace(0.0, z_screen, nz)

    # Initial fields
    beam = gaussian_beam_envelope(x, beam_waist)
    A = rectangular_slit(x, -slit_sep/2, slit_width)
    B = rectangular_slit(x, +slit_sep/2, slit_width)
    psi_A0, psi_B0 = beam*A, beam*B
    psi0 = psi_A0 + psi_B0
    psi_S0  = (psi_A0 + psi_B0)/np.sqrt(2)
    psi_AS0 = (psi_A0 - psi_B0)/np.sqrt(2)

    # Allocate
    I_total = np.zeros((nz, nx))
    I_A = np.zeros_like(I_total)
    I_B = np.zeros_like(I_total)
    I_S = np.zeros_like(I_total)
    I_AS = np.zeros_like(I_total)
    interf = np.zeros_like(I_total)

    # Propagate
    for i, z in enumerate(z_slices):
        psi   = angular_spectrum_propagate(psi0,   x, wavelength, z)
        psi_A = angular_spectrum_propagate(psi_A0, x, wavelength, z)
        psi_B = angular_spectrum_propagate(psi_B0, x, wavelength, z)
        psi_S = angular_spectrum_propagate(psi_S0, x, wavelength, z)
        psi_AS= angular_spectrum_propagate(psi_AS0,x, wavelength, z)

        I_total[i,:] = np.abs(psi)**2
        I_A[i,:]     = np.abs(psi_A)**2
        I_B[i,:]     = np.abs(psi_B)**2
        I_S[i,:]     = np.abs(psi_S)**2
        I_AS[i,:]    = np.abs(psi_AS)**2
        interf[i,:]  = 2*np.real(psi_A * np.conj(psi_B))

    # Output dir
    os.makedirs("out", exist_ok=True)

    # 3D surface
    X, Z = np.meshgrid(x*1e3, z_slices)
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Z, I_total, linewidth=0, antialiased=True)
    ax.set_xlabel("x (mm)"); ax.set_ylabel("z (m)"); ax.set_zlabel("Intensity")
    ax.set_title("Double-slit |ψ(x,z)|²")
    fig.savefig("out/intensity_3d_surface.png", dpi=200)
    plt.close(fig)

    # Probability channels & interference
    vmax = np.percentile(I_total, 99.5)
    fig, axs = plt.subplots(3,1, figsize=(10,9), sharex=True)
    axs[0].imshow(I_A, aspect="auto", extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)
    axs[0].set_title("|ψ_A|²")
    axs[1].imshow(I_B, aspect="auto", extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)
    axs[1].set_title("|ψ_B|²")
    axs[2].imshow(interf, aspect="auto", extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]])
    axs[2].set_title("2 Re(ψ_A ψ_B*)")
    for a in axs: a.set_ylabel("z (m)")
    axs[-1].set_xlabel("x (mm)")
    fig.savefig("out/channels_probability.png", dpi=200)
    plt.close(fig)

    # Symmetric / antisymmetric
    fig, axs = plt.subplots(2,1, figsize=(10,6), sharex=True)
    axs[0].imshow(I_S, aspect="auto", extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)
    axs[0].set_title("|ψ_S|² (symmetric)")
    axs[1].imshow(I_AS, aspect="auto", extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)
    axs[1].set_title("|ψ_AS|² (antisymmetric)")
    fig.savefig("out/sym_antisym_components.png", dpi=200)
    plt.close(fig)

    # Dominance map
    dominance = (I_AS > I_S).astype(float)
    plt.imshow(dominance, aspect="auto", extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=1)
    plt.title("Dominance (0=symmetric, 1=antisymmetric)")
    plt.xlabel("x (mm)"); plt.ylabel("z (m)")
    plt.savefig("out/dominance_map.png", dpi=200)
    plt.close()

    # Screen profile
    I_screen = I_total[-1,:]
    I_norm = I_screen / I_screen.max()
    peaks, _ = find_peaks(I_norm, prominence=0.02)
    c = peaks[np.argmin(np.abs(x[peaks]))] if peaks.size>0 else None
    hl = []
    if c is not None:
        left = peaks[peaks<c][-3:]
        right= peaks[peaks>c][:3]
        hl = np.concatenate([left,[c],right])
    plt.plot(x*1e3, I_screen)
    if len(hl)>0: plt.plot(x[hl]*1e3, I_screen[hl], "o")
    plt.xlabel("x (mm)"); plt.ylabel("Intensity")
    plt.title(f"Screen profile at z={z_screen} m")
    plt.savefig("out/screen_profile.png", dpi=200)
    plt.close()

    print("Simulation done. PNGs saved in ./out")

if __name__ == "__main__":
    main()

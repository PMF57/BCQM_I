#!/usr/bin/env python3\
# bcqm_double_slit_sim.py\
# Double-slit propagation with channel separation (A/B), interference term,\
# symmetric/antisymmetric "eigenmodes", and 3D intensity surface.\
#\
# Dependencies: numpy, scipy, matplotlib\
# Run: python bcqm_double_slit_sim.py [--wavelength 633e-9 --slit-sep 120e-6 --slit-width 20e-6 --z-screen 1.4 --nz 220]\
\
import os\
import argparse\
import numpy as np\
from numpy.fft import fft, ifft, fftfreq, fftshift\
import matplotlib.pyplot as plt\
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401\
from scipy.signal import find_peaks\
\
def angular_spectrum_propagate(psi0, x, wavelength, z):\
    """\
    1D angular spectrum propagation from z=0 to z>0 (paraxial-friendly exact scalar solution).\
    psi0: field at z=0 (array over x)\
    x: coordinate array (meters)\
    wavelength: wavelength (meters)\
    z: propagation distance (meters)\
    """\
    k = 2*np.pi / wavelength\
    dx = x[1] - x[0]\
    N = x.size\
\
    fx = fftfreq(N, d=dx)  # spatial frequency (cycles/m)\
    kx = 2*np.pi*fx        # angular spatial frequency\
    kz_sq = k**2 - kx**2\
    # Evanescent handling: clamp negative kz^2 to zero phase advance (or damp). We choose to zero evanescent phase (no growth).\
    kz = np.sqrt(np.clip(kz_sq, 0.0, None))\
    H = np.exp(1j * kz * z)\
\
    Psi0 = fft(psi0)\
    Psiz = Psi0 * H\
    psz = ifft(Psiz)\
    return psz\
\
def gaussian_beam_envelope(x, w0):\
    return np.exp(- (x**2) / (w0**2))\
\
def rectangular_slit(x, center, width):\
    return ((x >= center - width/2) & (x <= center + width/2)).astype(float)\
\
def gaussian_slit(x, center, width):\
    # width ~ 1/e half-width\
    return np.exp(- ((x - center)**2) / (2*(width/2.355)**2))  # width interpreted as FWHM\
\
def build_double_slit_field(x, wavelength, slit_sep, slit_width, beam_waist, slit_shape="rect"):\
    """\
    Build initial field at z=0: incident Gaussian beam times two slits.\
    slit_sep: center-to-center separation\
    slit_width: physical width (rect) or FWHM (gaussian)\
    """\
    beam = gaussian_beam_envelope(x, beam_waist)\
\
    cA = -slit_sep/2\
    cB = +slit_sep/2\
\
    if slit_shape == "rect":\
        A = rectangular_slit(x, cA, slit_width)\
        B = rectangular_slit(x, cB, slit_width)\
    else:\
        A = gaussian_slit(x, cA, slit_width)\
        B = gaussian_slit(x, cB, slit_width)\
\
    psi_A0 = beam * A\
    psi_B0 = beam * B\
    psi0 = psi_A0 + psi_B0\
    return psi0, psi_A0, psi_B0\
\
def mode_decomposition_symmetric_antisymmetric(psi_A0, psi_B0):\
    # Normalize by sqrt(2) for orthonormal-like splitting\
    psi_S0  = (psi_A0 + psi_B0)/np.sqrt(2)\
    psi_AS0 = (psi_A0 - psi_B0)/np.sqrt(2)\
    return psi_S0, psi_AS0\
\
def ensure_outdir(path="out"):\
    if not os.path.exists(path):\
        os.makedirs(path)\
    return path\
\
def main():\
    parser = argparse.ArgumentParser(description="Double-slit BCQM-style channels & modes simulation.")\
    parser.add_argument("--wavelength", type=float, default=633e-9, help="Wavelength (m). Default 633 nm (HeNe).")\
    parser.add_argument("--slit-sep", type=float, default=120e-6, help="Slit separation center-to-center (m).")\
    parser.add_argument("--slit-width", type=float, default=20e-6, help="Slit width (m) if rect; FWHM if gaussian.")\
    parser.add_argument("--slit-shape", type=str, default="rect", choices=["rect","gauss"], help="Slit transmission profile.")\
    parser.add_argument("--beam-waist", type=float, default=0.6e-3, help="Incident Gaussian 1/e radius w0 (m).")\
    parser.add_argument("--x-span", type=float, default=2.0e-3, help="Transverse simulation half-span (m). Total width=2*x-span.")\
    parser.add_argument("--nx", type=int, default=2048, help="Number of transverse samples.")\
    parser.add_argument("--z-screen", type=float, default=1.4, help="Screen distance (m).")\
    parser.add_argument("--nz", type=int, default=220, help="Number of z-slices from 0 to z-screen.")\
    parser.add_argument("--peak-prominence", type=float, default=0.02, help="Prominence fraction for fringe peak picking.")\
    args = parser.parse_args()\
    
#     class _Args: pass
#     args = _Args()
#     args.wavelength     = 633e-9
#     args.slit_sep       = 120e-6
#     args.slit_width     = 20e-6
#     args.slit_shape     = "rect"     # or "gauss"
#     args.beam_waist     = 0.6e-3
#     args.x_span         = 2.0e-3
#     args.nx             = 2048
#     args.z_screen       = 1.4
#     args.nz             = 220
#     args.peak_prominence= 0.02
    
\
    outdir = ensure_outdir()\
\
    # Grid\
    x = np.linspace(-args.x_span, args.x_span, args.nx)\
    z_slices = np.linspace(0.0, args.z_screen, args.nz)\
\
    # Fields at z=0\
    psi0, psi_A0, psi_B0 = build_double_slit_field(\
        x,\
        args.wavelength,\
        args.slit_sep,\
        args.slit_width,\
        args.beam_waist,\
        slit_shape=args.slit_shape\
    )\
\
    # Sym/antisym decomposition at aperture\
    psi_S0, psi_AS0 = mode_decomposition_symmetric_antisymmetric(psi_A0, psi_B0)\
\
    # Propagate all components across z-slices\
    I_total = np.zeros((args.nz, args.nx))\
    I_A = np.zeros_like(I_total)\
    I_B = np.zeros_like(I_total)\
    I_S = np.zeros_like(I_total)\
    I_AS = np.zeros_like(I_total)\
    interf = np.zeros_like(I_total, dtype=float)  # 2 Re(psi_A * psi_B*)\
\
    for i, z in enumerate(z_slices):\
        psi_z   = angular_spectrum_propagate(psi0,   x, args.wavelength, z)\
        psi_Az  = angular_spectrum_propagate(psi_A0, x, args.wavelength, z)\
        psi_Bz  = angular_spectrum_propagate(psi_B0, x, args.wavelength, z)\
        psi_Sz  = angular_spectrum_propagate(psi_S0, x, args.wavelength, z)\
        psi_ASz = angular_spectrum_propagate(psi_AS0,x, args.wavelength, z)\
\
        I_total[i,:] = np.abs(psi_z)**2\
        I_A[i,:]     = np.abs(psi_Az)**2\
        I_B[i,:]     = np.abs(psi_Bz)**2\
        I_S[i,:]     = np.abs(psi_Sz)**2\
        I_AS[i,:]    = np.abs(psi_ASz)**2\
        interf[i,:]  = 2.0 * np.real(psi_Az * np.conj(psi_Bz))\
\
    # 3D surface of total intensity\
    X, Z = np.meshgrid(x*1e3, z_slices)  # x in mm for readability\
    fig = plt.figure(figsize=(10,6))\
    ax = fig.add_subplot(111, projection='3d')\
    surf = ax.plot_surface(X, Z, I_total, linewidth=0, antialiased=True)\
    ax.set_xlabel("x (mm)")\
    ax.set_ylabel("z (m)")\
    ax.set_zlabel("Intensity (arb.)")\
    ax.set_title("Double-slit intensity |\uc0\u968 (x,z)|\'b2")\
    plt.tight_layout()\
    fig.savefig(os.path.join(outdir, "intensity_3d_surface.png"), dpi=200)\
    plt.close(fig)\
\
    # Channels (probability |A|^2, |B|^2) and interference term\
    vmax = np.percentile(I_total, 99.5)\
    fig, axs = plt.subplots(3,1, figsize=(10,10), sharex=True)\
    im0 = axs[0].imshow(I_A, aspect='auto', extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)\
    axs[0].set_title("|\uc0\u968 _A|\'b2 (probability channel A)")\
    im1 = axs[1].imshow(I_B, aspect='auto', extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)\
    axs[1].set_title("|\uc0\u968 _B|\'b2 (probability channel B)")\
    im2 = axs[2].imshow(interf, aspect='auto', extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]])\
    axs[2].set_title("Interference term 2 Re(\uc0\u968 _A \u968 _B*)")\
    for a in axs:\
        a.set_ylabel("z (m)")\
    axs[-1].set_xlabel("x (mm)")\
    plt.tight_layout()\
    fig.savefig(os.path.join(outdir, "channels_probability.png"), dpi=200)\
    plt.close(fig)\
\
    # Symmetric / Antisymmetric components\
    fig, axs = plt.subplots(2,1, figsize=(10,7), sharex=True)\
    im0 = axs[0].imshow(I_S, aspect='auto', extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)\
    axs[0].set_title("|\uc0\u968 _S|\'b2 (symmetric mode)")\
    im1 = axs[1].imshow(I_AS, aspect='auto', extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=vmax)\
    axs[1].set_title("|\uc0\u968 _AS|\'b2 (antisymmetric mode)")\
    for a in axs:\
        a.set_ylabel("z (m)")\
    axs[-1].set_xlabel("x (mm)")\
    plt.tight_layout()\
    fig.savefig(os.path.join(outdir, "sym_antisym_components.png"), dpi=200)\
    plt.close(fig)\
\
    # Dominance map: 0 if |\uc0\u968 _S|\'b2 >= |\u968 _AS|\'b2, 1 otherwise\
    dominance = (I_AS > I_S).astype(float)\
    fig, ax = plt.subplots(figsize=(10,4))\
    im = ax.imshow(dominance, aspect='auto', extent=[x[0]*1e3, x[-1]*1e3, z_slices[-1], z_slices[0]], vmin=0, vmax=1)\
    ax.set_title("Dominance map (0: symmetric, 1: antisymmetric)")\
    ax.set_ylabel("z (m)")\
    ax.set_xlabel("x (mm)")\
    plt.tight_layout()\
    fig.savefig(os.path.join(outdir, "dominance_map.png"), dpi=200)\
    plt.close(fig)\
\
    # Screen profile with peak finding (mark central + 3 either side if available)\
    I_screen = I_total[-1,:]\
    # Normalize for robust prominence thresholding\
    I_norm = I_screen / (I_screen.max() + 1e-12)\
    peaks, props = find_peaks(I_norm, prominence=args.peak_prominence)\
    # Sort by x and find central (closest to x=0)\
    if peaks.size > 0:\
        central_idx = peaks[np.argmin(np.abs(x[peaks]))]\
        # Get 3 either side if available\
        left = peaks[peaks < central_idx]\
        right = peaks[peaks > central_idx]\
        left_picks = left[-3:] if left.size >= 3 else left\
        right_picks = right[:3] if right.size >= 3 else right\
        highlight = np.concatenate([left_picks, np.array([central_idx]), right_picks])\
        highlight.sort()\
    else:\
        highlight = np.array([], dtype=int)\
\
    fig, ax = plt.subplots(figsize=(10,4))\
    ax.plot(x*1e3, I_screen, lw=1.5)\
    if highlight.size > 0:\
        ax.plot(x[highlight]*1e3, I_screen[highlight], 'o', ms=6)\
    ax.set_xlabel("x (mm)")\
    ax.set_ylabel("Intensity at screen (arb.)")\
    ax.set_title(f"Screen profile at z = \{args.z_screen\} m (peaks marked)")\
    plt.tight_layout()\
    fig.savefig(os.path.join(outdir, "screen_profile.png"), dpi=200)\
    plt.close(fig)\
\
    # Console guidance\
    print(f"Done. Outputs saved to: \{outdir\}")\
    print("Files:")\
    for f in ["intensity_3d_surface.png",\
              "channels_probability.png",\
              "sym_antisym_components.png",\
              "dominance_map.png",\
              "screen_profile.png"]:\
        print("  -", f)\
\
    # Quick tip to adjust fringe count:\
    print("\\nTips:")\
    print("- Increase slit separation (--slit-sep) or screen distance (--z-screen) to get more tightly spaced fringes.")\
    print("- Decrease slit width (--slit-width) to sharpen fringes.")\
    print("- Use --slit-shape gauss if you prefer Gaussian apertures.")\
\
if __name__ == "__main__":\
    main()

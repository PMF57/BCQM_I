from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt

def bloch_ode_step(r: np.ndarray, dt: float, gamma_phi: float, gamma_relax: float) -> np.ndarray:
    rx, ry, rz = r
    Gx = -(gamma_phi + 0.5*gamma_relax) * rx
    Gy = -(gamma_phi + 0.5*gamma_relax) * ry
    Gz = -gamma_relax * (rz + 1.0)
    return r + dt * np.array([Gx, Gy, Gz], dtype=float)

def simulate_gkls(t: np.ndarray, r0: np.ndarray = np.array([1.0, 0.0, 0.0]), gamma_phi: float = 0.0, gamma_relax: float = 0.0) -> np.ndarray:
    r = np.zeros((len(t), 3), dtype=float)
    r[0] = r0
    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        r[i] = bloch_ode_step(r[i-1], dt, gamma_phi, gamma_relax)
    return r

def plot_bloch_components(t: np.ndarray, r: np.ndarray, out_png: str) -> None:
    plt.figure()
    plt.plot(t, r[:,0], label="r_x")
    plt.plot(t, r[:,1], label="r_y")
    plt.plot(t, r[:,2], label="r_z")
    plt.xlabel("t (s)")
    plt.ylabel("Bloch components")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()

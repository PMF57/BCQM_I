from __future__ import annotations
import numpy as np
from typing import Callable, Optional, Tuple
import matplotlib.pyplot as plt

def compute_D(t: np.ndarray, gamma: float | None = None, Gamma_t: Optional[Callable[[float], float]] = None) -> np.ndarray:
    """Coherence factor D(t) = exp(-∫_0^t Γ(s) ds)."""
    t = np.asarray(t, dtype=float)
    if gamma is None and Gamma_t is None:
        raise ValueError("Provide either constant gamma or Gamma_t function")
    if gamma is not None and Gamma_t is not None:
        raise ValueError("Provide only one of gamma or Gamma_t")
    if gamma is not None:
        return np.exp(-gamma * t)
    G = np.zeros_like(t)
    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        G[i] = G[i-1] + 0.5 * dt * (Gamma_t(t[i]) + Gamma_t(t[i-1]))
    return np.exp(-G)

def F_opt_from_D(D: np.ndarray) -> np.ndarray:
    return 0.5 * (1.0 + D)

def W_from_D(t: np.ndarray, D: np.ndarray, F_star: float) -> Tuple[float | None, int]:
    bound = F_opt_from_D(D)
    mask = bound <= F_star
    if not np.any(mask):
        return None, -1
    idx = int(np.argmax(mask))
    return float(t[idx]), idx

def plot_ramsey(t: np.ndarray, D: np.ndarray, F_star: float, out_png: str) -> None:
    Fbound = F_opt_from_D(D)
    plt.figure()
    plt.plot(t, D, label="D(t)")
    plt.plot(t, Fbound, label="0.5*(1+D(t))")
    plt.axhline(F_star, linestyle="--", label="F*")
    plt.xlabel("t (s)")
    plt.ylabel("value")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()

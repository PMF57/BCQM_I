from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

def galton_hits(rows: int, n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.binomial(rows, 0.5, size=n)

def galton_histogram(hits: np.ndarray, rows: int) -> Tuple[np.ndarray, np.ndarray]:
    import numpy as np
    bins = np.arange(-0.5, rows+1.5, 1.0)
    counts, edges = np.histogram(hits, bins=bins, density=False)
    centers = 0.5*(edges[1:] + edges[:-1])
    return centers, counts

def plot_galton(centers, counts, rows: int, out_png: str) -> None:
    plt.figure()
    plt.bar(centers, counts, width=0.9, align="center")
    plt.xlabel("Bin index (0..rows)")
    plt.ylabel("Counts")
    plt.title(f"Galton distribution (rows={rows})")
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()

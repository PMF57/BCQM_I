from __future__ import annotations
import argparse, os
import numpy as np
from .ramsey import compute_D, F_opt_from_D, W_from_D, plot_ramsey
from .gkls import simulate_gkls, plot_bloch_components
from .galton import galton_hits, galton_histogram, plot_galton
from .utils import save_metadata, env_info

def main():
    parser = argparse.ArgumentParser(prog="bcqm", description="BCQM reproducible tools")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_ramsey = sub.add_parser("ramsey", help="Compute D(t), F_opt, and W bound")
    p_ramsey.add_argument("--gamma", type=float, default=None, help="Constant dephasing rate (1/s)")
    p_ramsey.add_argument("--fstar", type=float, required=True, help="Target fidelity threshold F* in (0.5,1)")
    p_ramsey.add_argument("--tmax", type=float, default=50e-6, help="Max time (s)")
    p_ramsey.add_argument("--dt", type=float, default=0.2e-6, help="Time step (s)")
    p_ramsey.add_argument("--out", type=str, default="bcqm_tools/examples/ramsey", help="Output prefix path")

    p_gkls = sub.add_parser("gkls", help="Simulate GKLS qubit (phase-covariant)")
    p_gkls.add_argument("--gamma-phi", type=float, default=0.0, help="Pure dephasing rate (1/s)")
    p_gkls.add_argument("--gamma-relax", type=float, default=0.0, help="Amplitude damping rate (1/s)")
    p_gkls.add_argument("--tmax", type=float, default=50e-6, help="Max time (s)")
    p_gkls.add_argument("--dt", type=float, default=0.5e-6, help="Time step (s)")
    p_gkls.add_argument("--out", type=str, default="bcqm_tools/examples/gkls", help="Output prefix path")

    p_galton = sub.add_parser("galton", help="Run Galton channels toy model")
    p_galton.add_argument("--rows", type=int, default=10, help="Number of rows (steps)")
    p_galton.add_argument("--n", type=int, default=1000, help="Number of particles")
    p_galton.add_argument("--seed", type=int, default=0, help="RNG seed")
    p_galton.add_argument("--out", type=str, default="bcqm_tools/examples/galton", help="Output prefix path")

    args = parser.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    if args.cmd == "ramsey":
        t = np.arange(0.0, args.tmax + 1e-15, args.dt)
        D = compute_D(t, gamma=args.gamma)
        Fbound = F_opt_from_D(D)
        W, idx = W_from_D(t, D, args.fstar)
        np.savetxt(f"{args.out}_t.csv", t, delimiter=",")
        np.savetxt(f"{args.out}_D.csv", D, delimiter=",")
        np.savetxt(f"{args.out}_Fbound.csv", Fbound, delimiter=",")
        plot_ramsey(t, D, args.fstar, f"{args.out}.png")
        meta = {"task":"ramsey","gamma":args.gamma,"fstar":args.fstar,"tmax":args.tmax,"dt":args.dt,"W":None if W is None else float(W),"env":env_info()}
        save_metadata(f"{args.out}.json", meta)
        print(f"W result: {W if W is not None else 'not reached within tmax'}")

    elif args.cmd == "gkls":
        t = np.arange(0.0, args.tmax + 1e-15, args.dt)
        r = simulate_gkls(t, gamma_phi=args.gamma_phi, gamma_relax=args.gamma_relax)
        np.savetxt(f"{args.out}_t.csv", t, delimiter=",")
        np.savetxt(f"{args.out}_r.csv", r, delimiter=",")
        plot_bloch_components(t, r, f"{args.out}.png")
        meta = {"task":"gkls","gamma_phi":args.gamma_phi,"gamma_relax":args.gamma_relax,"tmax":args.tmax,"dt":args.dt,"env":env_info()}
        save_metadata(f"{args.out}.json", meta)
        print("GKLS simulation done.")

    elif args.cmd == "galton":
        hits = galton_hits(args.rows, args.n, seed=args.seed)
        centers, counts = galton_histogram(hits, args.rows)
        np.savetxt(f"{args.out}_centers.csv", centers, delimiter=",")
        np.savetxt(f"{args.out}_counts.csv", counts, delimiter=",")
        plot_galton(centers, counts, args.rows, f"{args.out}.png")
        meta = {"task":"galton","rows":args.rows,"n":args.n,"seed":args.seed,"env":env_info()}
        save_metadata(f"{args.out}.json", meta)
        print("Galton simulation done.")

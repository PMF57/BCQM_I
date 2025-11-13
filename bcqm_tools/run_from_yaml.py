from __future__ import annotations
import argparse, yaml, numpy as np, os
from .ramsey import compute_D, F_opt_from_D, W_from_D, plot_ramsey
from .gkls import simulate_gkls, plot_bloch_components
from .galton import galton_hits, galton_histogram, plot_galton
from .utils import save_metadata, env_info

def main():
    ap = argparse.ArgumentParser(description="Run bcqm_tools tasks from a YAML config")
    ap.add_argument("--config", required=True, help="Path to YAML config")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    task = cfg.get("task")
    if task == "ramsey":
        gamma = cfg.get("gamma")
        fstar = cfg["fstar"]
        tmax = cfg.get("tmax", 50e-6)
        dt   = cfg.get("dt", 0.2e-6)
        out  = cfg.get("out", "bcqm_tools/examples/ramsey")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        t = np.arange(0.0, tmax + 1e-15, dt)
        D = compute_D(t, gamma=gamma)
        Fbound = F_opt_from_D(D)
        W, idx = W_from_D(t, D, fstar)
        np.savetxt(f"{out}_t.csv", t, delimiter=",")
        np.savetxt(f"{out}_D.csv", D, delimiter=",")
        np.savetxt(f"{out}_Fbound.csv", Fbound, delimiter=",")
        plot_ramsey(t, D, fstar, f"{out}.png")
        meta = {"task":"ramsey","gamma":gamma,"fstar":fstar,"tmax":tmax,"dt":dt,"W":None if W is None else float(W),"env":env_info()}
        save_metadata(f"{out}.json", meta)
        print(f"W result: {W if W is not None else 'not reached within tmax'}")

    elif task == "gkls":
        gamma_phi   = cfg.get("gamma_phi", 0.0)
        gamma_relax = cfg.get("gamma_relax", 0.0)
        tmax = cfg.get("tmax", 50e-6)
        dt   = cfg.get("dt", 0.5e-6)
        out  = cfg.get("out", "bcqm_tools/examples/gkls")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        t = np.arange(0.0, tmax + 1e-15, dt)
        r = simulate_gkls(t, gamma_phi=gamma_phi, gamma_relax=gamma_relax)
        np.savetxt(f"{out}_t.csv", t, delimiter=",")
        np.savetxt(f"{out}_r.csv", r, delimiter=",")
        plot_bloch_components(t, r, f"{out}.png")
        meta = {"task":"gkls","gamma_phi":gamma_phi,"gamma_relax":gamma_relax,"tmax":tmax,"dt":dt,"env":env_info()}
        save_metadata(f"{out}.json", meta)
        print("GKLS simulation done.")

    elif task == "galton":
        rows = cfg.get("rows", 10)
        n    = cfg.get("n", 1000)
        seed = cfg.get("seed", 0)
        out  = cfg.get("out", "bcqm_tools/examples/galton")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        hits = galton_hits(rows, n, seed=seed)
        centers, counts = galton_histogram(hits, rows)
        np.savetxt(f"{out}_centers.csv", centers, delimiter=",")
        np.savetxt(f"{out}_counts.csv", counts, delimiter=",")
        plot_galton(centers, counts, rows, f"{out}.png")
        meta = {"task":"galton","rows":rows,"n":n,"seed":seed,"env":env_info()}
        save_metadata(f"{out}.json", meta)
        print("Galton simulation done.")
    else:
        raise SystemExit(f"Unknown task in YAML: {task!r}")

if __name__ == "__main__":
    main()

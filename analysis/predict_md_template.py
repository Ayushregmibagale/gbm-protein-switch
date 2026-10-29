#!/usr/bin/env python3
"""
Template-guided Boltz-2 prediction for Track 1 caspase geometry check.

Uses the 10 ns MD frame as a structural template so Boltz-2 starts from the
assembled state rather than folding from scratch. If the halves remain close
(<15 Å) in the output, the MD result is more credible. If they separate,
that is a weaker sign.

Three seeds are run; endpoint C419–D420 CA distance is reported for each.
Results written to evidence/template_guided_controls/.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJECT  = Path("/home/ayushregmibagale/gbm-protein-switch")
BOLTZ    = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")
TEMPLATE = PROJECT / "evidence" / "track1_bound_10ns_template.pdb"
BOTH_YAML = PROJECT / "evidence" / "matched_controls" / "inputs" / "BOTH.yaml"
OUT_DIR  = PROJECT / "evidence" / "template_guided_controls"

# Endpoint residue indices in the Boltz output PDB (chain C last res = 250th,
# chain D first res = 1st).  Boltz renumbers per-chain from 1.
# Chain C: 250 residues → last = 250.  Chain D: first = 1.
CA_C_RES = 250   # last residue of chain C (caspase_N terminus)
CA_D_RES = 1     # first residue of chain D (caspase_C N-terminus)
PROXIMITY_CUTOFF = 15.0


def write_yaml(path: Path, template_path: Path) -> None:
    """Write YAML with the same sequences as BOTH.yaml + template block.

    Boltz-2 schema uses 'pdb:' (not 'path:') and 'chain_id:' (not 'ids:').
    """
    original = BOTH_YAML.read_text()
    yaml_text = original.rstrip("\n") + f"""
templates:
- pdb: {template_path}
  chain_id: [A, B, C, D]
"""
    path.write_text(yaml_text)


def _pred_dir(work_dir: Path) -> Path:
    """Boltz-2 names its output dir after the YAML stem, e.g.
    <out_dir>/boltz_results_<stem>/predictions/<stem>/
    Fall back to the old 'input' convention used by other scripts.
    """
    for result_dir in sorted(work_dir.glob("boltz_results_*")):
        preds = result_dir / "predictions"
        if preds.is_dir():
            for sub in preds.iterdir():
                if sub.is_dir() and any(sub.glob("*.pdb")):
                    return sub
            # Return first subdir even if no PDB yet (for existence check)
            subdirs = [s for s in preds.iterdir() if s.is_dir()]
            if subdirs:
                return subdirs[0]
    # Legacy fallback
    nested = work_dir / "boltz_results_input" / "predictions" / "input"
    flat   = work_dir / "predictions" / "input"
    return nested if nested.exists() else flat


def run_boltz(work_dir: Path, yaml_path: Path, seed: int) -> None:
    pd = _pred_dir(work_dir)
    if pd.exists() and any(pd.glob("*.pdb")):
        print(f"    seed{seed}: already done", flush=True)
        return
    work_dir.mkdir(parents=True, exist_ok=True)
    r = subprocess.run(
        [str(BOLTZ), "predict", str(yaml_path),
         "--out_dir", str(work_dir),
         "--recycling_steps", "3",
         "--sampling_steps", "200",
         "--output_format", "pdb",
         "--no_kernels",
         "--seed", str(seed)],
        capture_output=True, timeout=1800,
    )
    if r.returncode != 0:
        print(f"    seed{seed}: Boltz failed\n{r.stderr.decode()[:400]}", flush=True)


def get_iptm(work_dir: Path) -> float:
    pd = _pred_dir(work_dir)
    cfs = list(pd.glob("confidence_*_model_0.json")) if pd.exists() else []
    if not cfs:
        return 0.0
    return float(json.load(cfs[0].open()).get("iptm", 0.0))


def parse_ca(pdb_path: str, chain: str, resnum: int) -> np.ndarray | None:
    """Return the CA coordinate for the given chain/resnum, or None."""
    try:
        with open(pdb_path) as f:
            for line in f:
                if not (line.startswith("ATOM") and line[12:16].strip() == "CA"):
                    continue
                if line[21] != chain:
                    continue
                if int(line[22:26]) == resnum:
                    return np.array([float(line[30:38]),
                                     float(line[38:46]),
                                     float(line[46:54])])
    except Exception:
        pass
    return None


def endpoint_distance(pdb_path: str) -> float | None:
    ca_c = parse_ca(pdb_path, "C", CA_C_RES)
    ca_d = parse_ca(pdb_path, "D", CA_D_RES)
    if ca_c is None or ca_d is None:
        return None
    return float(np.linalg.norm(ca_c - ca_d))


def get_pdb(work_dir: Path) -> str:
    pd = _pred_dir(work_dir)
    pdbs = list(pd.glob("*_model_0.pdb")) if pd.exists() else []
    return str(pdbs[0]) if pdbs else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--status", action="store_true",
                        help="Report results without running new predictions")
    args = parser.parse_args()

    if not TEMPLATE.exists():
        sys.exit(f"Template PDB not found: {TEMPLATE}")
    if not BOTH_YAML.exists():
        sys.exit(f"BOTH.yaml not found: {BOTH_YAML}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    for seed in range(args.seeds):
        wd   = OUT_DIR / f"seed{seed}"
        yaml = OUT_DIR / f"seed{seed}_input.yaml"

        if not args.status:
            write_yaml(yaml, TEMPLATE)
            print(f"  Running seed {seed} ...", flush=True)
            run_boltz(wd, yaml, seed)

        iptm = get_iptm(wd)
        pdb  = get_pdb(wd)
        dist = endpoint_distance(pdb) if pdb else None
        close = dist is not None and dist < PROXIMITY_CUTOFF

        results.append({
            "seed": seed,
            "iptm": iptm,
            "endpoint_dist_A": dist,
            "halves_close": close,
            "pdb": pdb,
        })
        status = f"{dist:.1f} Å {'✓ close' if close else '✗ separated'}" if dist is not None else "no PDB"
        print(f"  seed{seed}: ipTM={iptm:.3f}  endpoint={status}")

    print("\n=== Template-guided prediction summary ===")
    print(f"  Template: {TEMPLATE.name}  (split site 9.54 Å in MD)")
    print(f"  Endpoint: chain-C res {CA_C_RES} CA  ↔  chain-D res {CA_D_RES} CA")
    print(f"  Proximity cutoff: {PROXIMITY_CUTOFF} Å\n")

    n_close = sum(r["halves_close"] for r in results)
    for r in results:
        status = f"{r['endpoint_dist_A']:.1f} Å" if r["endpoint_dist_A"] is not None else "n/a"
        flag   = "✓ close" if r["halves_close"] else "✗ separated"
        print(f"  seed{r['seed']}: ipTM={r['iptm']:.3f}  dist={status}  {flag}")

    print(f"\n  {n_close}/{len(results)} seeds keep halves close (<{PROXIMITY_CUTOFF} Å)")
    if n_close == len(results):
        print("  → Boltz-2 keeps the assembly: MD result is more credible")
    elif n_close == 0:
        print("  → Boltz-2 separates the halves even from the assembled template: concerning")
    else:
        print("  → Mixed: check individual PDBs")

    summary = {
        "template": str(TEMPLATE),
        "template_split_dist_A": 9.54,
        "endpoint_chain_C_res": CA_C_RES,
        "endpoint_chain_D_res": CA_D_RES,
        "cutoff_A": PROXIMITY_CUTOFF,
        "seeds": results,
        "n_close": n_close,
        "interpretation": (
            "assembly retained" if n_close == len(results) else
            "assembly lost" if n_close == 0 else
            "mixed"
        ),
    }
    out_json = OUT_DIR / "summary.json"
    json.dump(summary, out_json.open("w"), indent=2)
    print(f"\nSaved → {out_json}")


if __name__ == "__main__":
    main()

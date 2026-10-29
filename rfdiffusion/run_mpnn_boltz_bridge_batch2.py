#!/usr/bin/env python3
"""
MPNN + Boltz-2 pipeline for RFdiffusion bridge binders — batch 2.
Processes backbones ranked 21–80 by dual-hotspot CA distance (skipping the top 20
already done in batch 1). Outputs to bridge_mpnn_boltz_batch2/, accepted PDBs
also copied to the shared bridge_mpnn_boltz/Accepted/.

Usage:
  python run_mpnn_boltz_bridge_batch2.py [--batch_start N] [--batch_end M] [--mpnn_seqs N] [--temp T]
Defaults: batch_start=21  batch_end=80  mpnn_seqs=10  temp=0.15
"""

import argparse, csv, json, os, shutil, subprocess
from pathlib import Path
import numpy as np

# ── paths ──────────────────────────────────────────────────────────────────
GBM       = Path("/mnt/d/gbm-protein-switch")
RFD_OUT   = GBM / "rfdiffusion" / "bridge_both_chains"
BATCH1    = GBM / "rfdiffusion" / "bridge_mpnn_boltz"
OUT_DIR   = GBM / "rfdiffusion" / "bridge_mpnn_boltz_batch2"
ACCEPT    = BATCH1 / "Accepted"           # shared accepted dir
MPNN_DIR  = Path("/home/ayushregmibagale/projects/ProteinMPNN")
PYTHON_RF = "/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin/python"
BOLTZ     = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

SOX2_HOTSPOTS  = {193, 199}
SOX9_HOTSPOTS  = {466, 505, 508}
BOLTZ_IPTM_THRESH    = 0.40
BOLTZ_CONTACT_THRESH = 0.12


def parse_ca(pdb_path, chain):
    atoms = {}
    with open(pdb_path) as f:
        for l in f:
            if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == chain:
                res = int(l[22:26])
                atoms[res] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return atoms


def min_hotspot_dist(ca_binder, ca_target, hotspots):
    if not ca_binder or not ca_target:
        return 999.0
    c = np.array(list(ca_binder.values()))
    dists = [np.linalg.norm(c - ca_target[r], axis=1).min()
             for r in hotspots if r in ca_target]
    return min(dists) if dists else 999.0


def contact_fraction(ca_binder, ca_target, cutoff=8.0):
    if not ca_binder or not ca_target:
        return 0.0
    b = np.array(list(ca_binder.values()))
    t = np.array(list(ca_target.values()))
    return sum(1 for ba in b if np.linalg.norm(t - ba, axis=1).min() < cutoff) / len(b)


def rank_all_backbones():
    ranked = []
    for pdb in sorted(RFD_OUT.glob("bridge_*.pdb")):
        ca_A = parse_ca(pdb, "A")
        ca_B = parse_ca(pdb, "B")
        ca_C = parse_ca(pdb, "C")
        if not ca_C:
            continue
        s2 = min_hotspot_dist(ca_C, ca_A, SOX2_HOTSPOTS)
        s9 = min_hotspot_dist(ca_C, ca_B, SOX9_HOTSPOTS)
        ranked.append((max(s2, s9), s2, s9, len(ca_C), pdb))
    ranked.sort(key=lambda x: x[0])
    return ranked


def run_mpnn(pdb_path: Path, work_dir: Path, n_seqs: int, temp: str) -> list[str]:
    work_dir.mkdir(parents=True, exist_ok=True)
    parsed   = work_dir / "parsed.jsonl"
    filtered = work_dir / "filtered.jsonl"
    fixed    = work_dir / "fixed.jsonl"
    mpnn_out = work_dir / "mpnn_out"
    mpnn_out.mkdir(exist_ok=True)

    # Copy just this PDB to a temp dir so parse_multiple_chains only sees it
    tmp_pdb_dir = work_dir / "pdb_tmp"
    tmp_pdb_dir.mkdir(exist_ok=True)
    shutil.copy(pdb_path, tmp_pdb_dir / pdb_path.name)

    subprocess.run([PYTHON_RF,
        str(MPNN_DIR / "helper_scripts" / "parse_multiple_chains.py"),
        f"--input_path={tmp_pdb_dir}",
        f"--output_path={parsed}",
    ], check=True, capture_output=True)

    pdb_name = pdb_path.stem
    with open(parsed) as fin, open(filtered, "w") as fout:
        for line in fin:
            if pdb_name in line:
                fout.write(line)
                break

    if filtered.stat().st_size == 0:
        print(f"  WARN: {pdb_name} not in parsed JSONL")
        return []

    subprocess.run([PYTHON_RF,
        str(MPNN_DIR / "helper_scripts" / "assign_fixed_chains.py"),
        f"--input_path={filtered}",
        f"--output_path={fixed}",
        "--chain_list=C",
    ], check=True, capture_output=True)

    subprocess.run([PYTHON_RF,
        str(MPNN_DIR / "protein_mpnn_run.py"),
        f"--jsonl_path={filtered}",
        f"--chain_id_jsonl={fixed}",
        f"--out_folder={mpnn_out}",
        f"--num_seq_per_target={n_seqs}",
        f"--sampling_temp={temp}",
        "--seed=37",
        "--batch_size=1",
    ], check=True, capture_output=True)

    seqs = []
    for fa in (mpnn_out / "seqs").glob("*.fa"):
        with open(fa) as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith(">"):
                seq_full = lines[i+1].strip() if i+1 < len(lines) else ""
                if seq_full and not seq_full.startswith(">"):
                    chain_c = seq_full.split("/")[-1]
                    if chain_c:
                        seqs.append(chain_c)
                i += 2
            else:
                i += 1
    return seqs


def write_yaml(path: Path, binder_seq: str):
    path.write_text(f"""sequences:
  - protein:
      id: A
      sequence: {SOX2_TARGET}
      msa: empty
  - protein:
      id: B
      sequence: {SOX9_TARGET}
      msa: empty
  - protein:
      id: C
      sequence: {binder_seq}
      msa: empty
""")


def run_boltz(yaml_path: Path, out_dir: Path):
    try:
        subprocess.run([
            str(BOLTZ), "predict", str(yaml_path),
            "--out_dir", str(out_dir),
            "--recycling_steps", "3",
            "--sampling_steps", "200",
            "--output_format", "pdb",
            "--no_kernels",
        ], check=True, capture_output=True, timeout=600)
    except Exception as e:
        print(f"    Boltz error: {e}")
        return None, None

    pdb_files  = list(out_dir.rglob("*model_0.pdb"))
    conf_files = list(out_dir.rglob("confidence_*model_0.json"))
    pdb  = pdb_files[0]  if pdb_files  else None
    conf = json.loads(conf_files[0].read_text()) if conf_files else {}
    return conf, pdb


def score_pdb(pdb_path: Path) -> dict:
    ca_A = parse_ca(pdb_path, "A")
    ca_B = parse_ca(pdb_path, "B")
    ca_C = parse_ca(pdb_path, "C")
    return {
        "sox2_contact":      round(contact_fraction(ca_C, ca_A), 3),
        "sox9_contact":      round(contact_fraction(ca_C, ca_B), 3),
        "sox2_hotspot_dist": round(min_hotspot_dist(ca_C, ca_A, SOX2_HOTSPOTS), 2),
        "sox9_hotspot_dist": round(min_hotspot_dist(ca_C, ca_B, SOX9_HOTSPOTS), 2),
        "binder_len":        len(ca_C),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch_start", type=int, default=21)
    ap.add_argument("--batch_end",   type=int, default=80)
    ap.add_argument("--mpnn_seqs",   type=int, default=10)
    ap.add_argument("--temp",        type=str, default="0.15")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ACCEPT.mkdir(parents=True, exist_ok=True)

    print("=== Ranking all 200 backbones ===")
    all_ranked = rank_all_backbones()
    batch = all_ranked[args.batch_start - 1 : args.batch_end]
    print(f"Processing backbones {args.batch_start}–{min(args.batch_end, len(all_ranked))} "
          f"({len(batch)} total, {args.mpnn_seqs} seqs each at temp={args.temp})")
    print(f"Expected: {len(batch) * args.mpnn_seqs} Boltz-2 runs")

    results = []
    total_seqs = 0

    for rank_idx, (worst, s2_bb, s9_bb, blen_bb, pdb) in enumerate(batch):
        backbone  = pdb.stem
        rank_num  = args.batch_start + rank_idx
        work_dir  = OUT_DIR / backbone
        print(f"\n=== [{rank_num}] {backbone}  SOX2={s2_bb:.2f}Å SOX9={s9_bb:.2f}Å len={blen_bb} ===")

        print(f"  Running ProteinMPNN ({args.mpnn_seqs} seqs, temp={args.temp})...")
        seqs = run_mpnn(pdb, work_dir / "mpnn", args.mpnn_seqs, args.temp)
        print(f"  Got {len(seqs)} sequences")
        if not seqs:
            continue

        for si, seq in enumerate(seqs):
            total_seqs += 1
            tag       = f"{backbone}_s{si:02d}"
            boltz_dir = work_dir / f"boltz_{si:02d}"
            yaml_path = boltz_dir / "input.yaml"
            boltz_dir.mkdir(exist_ok=True)
            write_yaml(yaml_path, seq)

            print(f"  [{si+1}/{len(seqs)}] {tag} (len={len(seq)})... ", end="", flush=True)
            conf, pred_pdb = run_boltz(yaml_path, boltz_dir)

            if conf is None:
                print("FAILED")
                results.append({"tag": tag, "backbone": backbone, "seq_idx": si,
                                 "rank": rank_num, "sequence": seq, "status": "FAILED"})
                continue

            iptm  = conf.get("iptm", 0)
            plddt = conf.get("plddt", 0)
            ptm   = conf.get("ptm", 0)
            geo   = score_pdb(pred_pdb) if pred_pdb else {}
            sox2_c = geo.get("sox2_contact", 0)
            sox9_c = geo.get("sox9_contact", 0)
            s2_h   = geo.get("sox2_hotspot_dist", 999)
            s9_h   = geo.get("sox9_hotspot_dist", 999)

            passed = (iptm >= BOLTZ_IPTM_THRESH and
                      sox2_c >= BOLTZ_CONTACT_THRESH and
                      sox9_c >= BOLTZ_CONTACT_THRESH)

            print(f"ipTM={iptm:.3f} pLDDT={plddt:.3f} SOX2_c={sox2_c:.2f} "
                  f"SOX9_c={sox9_c:.2f}  {'✓ PASS' if passed else 'fail'}")

            if passed and pred_pdb:
                dst = ACCEPT / f"{tag}.pdb"
                shutil.copy(pred_pdb, dst)
                print(f"    → Accepted: {dst.name}")

            results.append({
                "tag": tag, "backbone": backbone, "seq_idx": si, "rank": rank_num,
                "sequence": seq, "iptm": round(iptm, 4), "plddt": round(plddt, 4),
                "ptm": round(ptm, 4), "sox2_contact": sox2_c, "sox9_contact": sox9_c,
                "sox2_hotspot_dist": s2_h, "sox9_hotspot_dist": s9_h,
                "bb_sox2_dist": round(s2_bb, 2), "bb_sox9_dist": round(s9_bb, 2),
                "binder_len": geo.get("binder_len", len(seq)),
                "passed": passed, "status": "OK",
            })

    tsv = OUT_DIR / "results.tsv"
    if results:
        with open(tsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0].keys()), delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(results)

    ok     = [r for r in results if r.get("status") == "OK"]
    passed = [r for r in ok if r.get("passed")]
    print(f"\n{'='*60}")
    print(f"Sequences tested: {total_seqs}  Boltz OK: {len(ok)}  Passed filter: {len(passed)}")
    if passed:
        passed.sort(key=lambda x: -(x["sox2_contact"] + x["sox9_contact"]))
        print(f"\n{'tag':<30} {'ipTM':>6} {'SOX2_c':>7} {'SOX9_c':>7} {'len':>5}")
        for r in passed:
            print(f"{r['tag']:<30} {r['iptm']:>6.3f} {r['sox2_contact']:>7.3f} "
                  f"{r['sox9_contact']:>7.3f} {r['binder_len']:>5}")
    print(f"\nResults: {tsv}")
    print(f"Accepted PDBs: {ACCEPT}")
    print("=== DONE ===")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
ProteinMPNN + Boltz-2 pipeline for RFdiffusion bridge binder designs.

Steps:
  1. Select top 20 dual-contact backbones from bridge_both_chains/
  2. ProteinMPNN: design chain C (binder), fix chains A (SOX2) and B (SOX9), 10 seqs/backbone
  3. Boltz-2: 3-chain prediction per sequence — SOX2_HMG + SOX9_HMG + binder
  4. Geometric filter: ipTM + dual-chain contact fractions
  5. Save results TSV + accepted PDBs

Usage:
  python run_mpnn_boltz_bridge.py
"""

import os, json, subprocess, shutil, csv, glob, pickle, tempfile
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

# ── paths ──────────────────────────────────────────────────────────────────
GBM        = Path("/mnt/d/gbm-protein-switch")
RFD_OUT    = GBM / "rfdiffusion" / "bridge_both_chains"
OUT_DIR    = GBM / "rfdiffusion" / "bridge_mpnn_boltz"
MPNN_DIR   = Path("/home/ayushregmibagale/projects/ProteinMPNN")
PYTHON_RF  = "/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin/python"
BOLTZ      = Path("/home/ayushregmibagale/miniconda3/envs/boltz2/bin/boltz")
PYTHON_B2  = "/home/ayushregmibagale/miniconda3/envs/boltz2/bin/python3"

# ── targets for Boltz-2 ────────────────────────────────────────────────────
SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

# ── selection criteria ──────────────────────────────────────────────────────
SOX2_HOTSPOTS  = {193, 199}
SOX9_HOTSPOTS  = {466, 505, 508}
CONTACT_CUTOFF = 8.0    # Å, for selecting top backbones
N_BACKBONES    = 20     # top N to run MPNN on
MPNN_SEQS      = 10     # sequences per backbone
MPNN_TEMP      = "0.1"  # sampling temperature (lower = more conservative)
BOLTZ_IPTM_THRESH    = 0.40
BOLTZ_CONTACT_THRESH = 0.12  # fraction of binder within 8 Å of each target chain

pdb_parser = PDBParser(QUIET=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── helper: parse Cα coords ────────────────────────────────────────────────
def parse_ca(pdb_path, chain):
    atoms = {}
    with open(pdb_path) as f:
        for l in f:
            if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == chain:
                res = int(l[22:26])
                atoms[res] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return atoms


def min_hotspot_dist(ca_binder, ca_target, hotspots):
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


# ── step 1: rank backbones ─────────────────────────────────────────────────
def rank_backbones():
    pdbs = sorted(RFD_OUT.glob("bridge_*.pdb"))
    ranked = []
    for pdb in pdbs:
        ca_A = parse_ca(pdb, "A")
        ca_B = parse_ca(pdb, "B")
        ca_C = parse_ca(pdb, "C")
        if not ca_C:
            continue
        s2 = min_hotspot_dist(ca_C, ca_A, SOX2_HOTSPOTS)
        s9 = min_hotspot_dist(ca_C, ca_B, SOX9_HOTSPOTS)
        ranked.append((max(s2, s9), s2, s9, len(ca_C), pdb))
    ranked.sort(key=lambda x: x[0])
    return ranked[:N_BACKBONES]


# ── step 2: ProteinMPNN ────────────────────────────────────────────────────
def run_mpnn(pdb_path: Path, work_dir: Path) -> list[str]:
    """Run MPNN on one PDB, design chain C, fix A and B. Returns list of sequences."""
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1) parse chains
    parsed = work_dir / "parsed.jsonl"
    subprocess.run([
        PYTHON_RF,
        str(MPNN_DIR / "helper_scripts" / "parse_multiple_chains.py"),
        f"--input_path={pdb_path.parent}",
        f"--output_path={parsed}",
    ], check=True, capture_output=True)

    # filter parsed.jsonl to only this PDB (parse_multiple_chains scans the whole dir)
    pdb_name = pdb_path.stem
    filtered = work_dir / "filtered.jsonl"
    with open(parsed) as fin, open(filtered, "w") as fout:
        for line in fin:
            if pdb_name in line:
                fout.write(line)
                break

    if filtered.stat().st_size == 0:
        print(f"  WARN: {pdb_name} not found in parsed JSONL")
        return []

    # 2) assign chain C as the designed chain (A, B are implicitly fixed)
    fixed = work_dir / "fixed.jsonl"
    subprocess.run([
        PYTHON_RF,
        str(MPNN_DIR / "helper_scripts" / "assign_fixed_chains.py"),
        f"--input_path={filtered}",
        f"--output_path={fixed}",
        '--chain_list=C',
    ], check=True, capture_output=True)

    # 3) run MPNN
    mpnn_out = work_dir / "mpnn_out"
    mpnn_out.mkdir(exist_ok=True)
    subprocess.run([
        PYTHON_RF,
        str(MPNN_DIR / "protein_mpnn_run.py"),
        f"--jsonl_path={filtered}",
        f"--chain_id_jsonl={fixed}",
        f"--out_folder={mpnn_out}",
        f"--num_seq_per_target={MPNN_SEQS}",
        f"--sampling_temp={MPNN_TEMP}",
        "--seed=42",
        "--batch_size=1",
    ], check=True, capture_output=True)

    # 4) parse FASTA output → extract chain C sequences
    seqs = []
    for fa in (mpnn_out / "seqs").glob("*.fa"):
        with open(fa) as f:
            lines = f.readlines()
        # MPNN FASTA: header line has chain info; sequence is interleaved
        i = 0
        while i < len(lines):
            if lines[i].startswith(">"):
                header = lines[i].strip()
                seq_full = lines[i+1].strip() if i+1 < len(lines) else ""
                # When only chain C is designed, MPNN outputs just its sequence (no / separator)
                if seq_full and not seq_full.startswith(">"):
                    chain_c = seq_full.split("/")[-1]  # handles both single-chain and multi-chain output
                    if chain_c:
                        seqs.append(chain_c)
                i += 2
            else:
                i += 1
    return seqs


# ── step 3: Boltz-2 ───────────────────────────────────────────────────────
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


def run_boltz(yaml_path: Path, out_dir: Path) -> tuple[dict | None, Path | None]:
    cmd = [
        str(BOLTZ), "predict", str(yaml_path),
        "--out_dir", str(out_dir),
        "--recycling_steps", "3",
        "--sampling_steps", "200",
        "--output_format", "pdb",
        "--no_kernels",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=600)
    except Exception as e:
        print(f"    Boltz error: {e}")
        return None, None

    pdb_files = list(out_dir.rglob("*model_0.pdb"))
    conf_files = list(out_dir.rglob("confidence_*model_0.json"))
    pdb = pdb_files[0] if pdb_files else None
    conf = json.loads(conf_files[0].read_text()) if conf_files else {}
    return conf, pdb


def score_boltz_pdb(pdb_path: Path) -> dict:
    ca_A = parse_ca(pdb_path, "A")
    ca_B = parse_ca(pdb_path, "B")
    ca_C = parse_ca(pdb_path, "C")
    return {
        "sox2_contact": round(contact_fraction(ca_C, ca_A), 3),
        "sox9_contact": round(contact_fraction(ca_C, ca_B), 3),
        "sox2_hotspot_dist": round(min_hotspot_dist(ca_C, ca_A, SOX2_HOTSPOTS), 2),
        "sox9_hotspot_dist": round(min_hotspot_dist(ca_C, ca_B, SOX9_HOTSPOTS), 2),
        "binder_len": len(ca_C),
    }


# ── main ───────────────────────────────────────────────────────────────────
def main():
    print("=== Step 1: Ranking backbones ===")
    top = rank_backbones()
    print(f"Selected {len(top)} backbones:")
    for worst, s2, s9, bl, pdb in top:
        print(f"  {pdb.stem}: SOX2={s2:.2f}Å SOX9={s9:.2f}Å len={bl}")

    results = []
    accepted_dir = OUT_DIR / "Accepted"
    accepted_dir.mkdir(exist_ok=True)

    total_seqs = 0
    for rank, (worst, s2_bb, s9_bb, blen_bb, pdb) in enumerate(top):
        backbone = pdb.stem
        print(f"\n=== [{rank+1}/{len(top)}] {backbone} (SOX2={s2_bb:.2f} SOX9={s9_bb:.2f}) ===")

        work_dir = OUT_DIR / backbone
        work_dir.mkdir(exist_ok=True)

        # MPNN
        print(f"  Running ProteinMPNN ({MPNN_SEQS} seqs)...")
        seqs = run_mpnn(pdb, work_dir / "mpnn")
        print(f"  Got {len(seqs)} sequences")
        if not seqs:
            continue

        # Boltz-2 for each sequence
        for si, seq in enumerate(seqs):
            total_seqs += 1
            tag = f"{backbone}_s{si:02d}"
            boltz_dir = work_dir / f"boltz_{si:02d}"
            yaml_path = boltz_dir / "input.yaml"
            boltz_dir.mkdir(exist_ok=True)
            write_yaml(yaml_path, seq)

            print(f"  [{si+1}/{len(seqs)}] {tag} (len={len(seq)})... ", end="", flush=True)
            conf, pred_pdb = run_boltz(yaml_path, boltz_dir)

            if conf is None:
                print("FAILED")
                results.append({"tag": tag, "backbone": backbone, "seq_idx": si,
                                 "sequence": seq, "status": "FAILED"})
                continue

            iptm  = conf.get("iptm", 0)
            plddt = conf.get("plddt", 0)
            ptm   = conf.get("ptm", 0)

            geo = score_boltz_pdb(pred_pdb) if pred_pdb else {}
            sox2_c = geo.get("sox2_contact", 0)
            sox9_c = geo.get("sox9_contact", 0)
            s2_h   = geo.get("sox2_hotspot_dist", 999)
            s9_h   = geo.get("sox9_hotspot_dist", 999)

            passed = (iptm >= BOLTZ_IPTM_THRESH and
                      sox2_c >= BOLTZ_CONTACT_THRESH and
                      sox9_c >= BOLTZ_CONTACT_THRESH)

            print(f"ipTM={iptm:.3f} pLDDT={plddt:.3f} SOX2_c={sox2_c:.2f} SOX9_c={sox9_c:.2f} "
                  f"hotspot: S2={s2_h:.1f}Å S9={s9_h:.1f}Å  {'✓ PASS' if passed else 'fail'}")

            if passed and pred_pdb:
                dst = accepted_dir / f"{tag}.pdb"
                shutil.copy(pred_pdb, dst)
                print(f"    → Saved to {dst}")

            results.append({
                "tag": tag, "backbone": backbone, "seq_idx": si, "sequence": seq,
                "iptm": round(iptm, 4), "plddt": round(plddt, 4), "ptm": round(ptm, 4),
                "sox2_contact": sox2_c, "sox9_contact": sox9_c,
                "sox2_hotspot_dist": s2_h, "sox9_hotspot_dist": s9_h,
                "bb_sox2_dist": round(s2_bb, 2), "bb_sox9_dist": round(s9_bb, 2),
                "binder_len": geo.get("binder_len", len(seq)),
                "passed": passed, "status": "OK",
            })

    # Save results
    tsv = OUT_DIR / "mpnn_boltz_results.tsv"
    if results:
        keys = list(results[0].keys())
        with open(tsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(results)

    ok      = [r for r in results if r.get("status") == "OK"]
    passed  = [r for r in ok if r.get("passed")]
    print(f"\n{'='*60}")
    print(f"Total sequences tested: {total_seqs}")
    print(f"Boltz-2 OK: {len(ok)}  |  Passed filter: {len(passed)}")
    if passed:
        passed.sort(key=lambda x: -(x["sox2_contact"] + x["sox9_contact"]))
        print("\nTop accepted designs:")
        print(f"{'tag':<30} {'ipTM':>6} {'SOX2_c':>7} {'SOX9_c':>7} {'S2_hot':>7} {'S9_hot':>7}")
        for r in passed[:10]:
            print(f"{r['tag']:<30} {r['iptm']:>6.3f} {r['sox2_contact']:>7.3f} "
                  f"{r['sox9_contact']:>7.3f} {r['sox2_hotspot_dist']:>7.1f} {r['sox9_hotspot_dist']:>7.1f}")
    print(f"\nResults: {tsv}")
    print(f"Accepted PDBs: {accepted_dir}")
    print("=== DONE ===")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Phase 2: Two-chain linker screen.
Fuse top stage-3-selective SOX2+SOX9 binder pairs into a single chain with GGS linkers.
Run Boltz-2 3-chain prediction: SOX2_HMG (A) + SOX9_HMG (B) + fused_binder (C).
Filter by ipTM and pLDDT; measure which orientation works.

Usage:
  python bridge_linker_screen.py
"""

from selection_support import require_execution, reference_sequence, validate_recorded_candidates, campaign_output, campaign_table
require_execution(__name__)

import os, json, subprocess
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

BASE     = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
OUT_DIR  = campaign_output(BASE / "selection" / "bridge_linker_outputs")
BOLTZ    = Path(os.environ.get("BOLTZ_BIN", "boltz"))

# Top selective pairs from stage 3 (both_selective=True, sorted by pair_ipTM)
PAIRS = [
    {
        "name":     "mpnn17_x_sox9run1",
        "sox2_seq": "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM",
        "sox9_seq": "PPLSPEEQELEEIEEGYRNLTEKSGEALLKALGITSKNEEEREKKMRELALELAREVTRKLSEKTGVRYAEDGTRLP",
    },
    {
        "name":     "mpnn18_x_sox9run4",
        "sox2_seq": "PYEIPEEDREWDEYSFRMVEVHERLMEARRELRLLAPEETEKKRELEEEIARAEKELEEVHKEGEAIERKM",
        "sox9_seq": "SEGIELMKELGLPVKDGKPYIEDLENLPELTNEEKLELGLLYGELNGGEKKAREAVEKAKKEGTEEAIERAFEEAVEAGDEGLKKK",
    },
    {
        "name":     "ph_sox2run18_x_sox9run18",
        "sox2_seq": "SRLEELVDHIVDSMEEANQGGLTDEAREMIRNNLLEDLLDVNPNVEAELDAASDAEFAELLAGVDRIIREGTTSDLLGFLEALAEGGVAGALAELAK",
        "sox9_seq": "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL",
    },
]

LINKERS = ["GGS", "GGSGGS", "GGSGGSGGS"]

pdb_parser = PDBParser(QUIET=True)
results = []


def write_yaml(path: Path, fused_seq: str):
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
      sequence: {fused_seq}
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
        print(f"  Boltz error: {e}")
        return None, None

    pdb_files = list(out_dir.rglob("*model_0.pdb"))
    if not pdb_files:
        return None, None
    pdb_path = pdb_files[0]

    conf_files = list(out_dir.rglob("confidence_*model_0.json"))
    if not conf_files:
        return None, pdb_path
    conf = json.loads(conf_files[0].read_text())
    return conf, pdb_path


def chain_contact_fraction(pdb_path: Path, binder_chain: str, target_chain: str, cutoff: float = 8.0) -> float:
    """Fraction of binder Cα atoms within cutoff Å of any target Cα."""
    try:
        struct = pdb_parser.get_structure("x", pdb_path)[0]
        b_cas = [r["CA"].get_vector().get_array() for r in struct[binder_chain].get_residues() if "CA" in r]
        t_cas = [r["CA"].get_vector().get_array() for r in struct[target_chain].get_residues() if "CA" in r]
        if not b_cas or not t_cas:
            return 0.0
        b = np.array(b_cas)
        t = np.array(t_cas)
        contacts = 0
        for ba in b:
            if np.min(np.linalg.norm(t - ba, axis=1)) < cutoff:
                contacts += 1
        return contacts / len(b)
    except Exception:
        return 0.0


def midpoint_split_dist(pdb_path: Path, sox2_len: int, sox9_len: int, linker_len: int) -> float | None:
    """
    Distance between the Cα centroid of the SOX2-binding half of the fused binder (chain C,
    first sox2_len residues) and the SOX9-binding half (last sox9_len residues).
    """
    try:
        struct = pdb_parser.get_structure("x", pdb_path)[0]
        residues = list(struct["C"].get_residues())
        sox2_half = [r["CA"].get_vector().get_array() for r in residues[:sox2_len] if "CA" in r]
        sox9_half = [r["CA"].get_vector().get_array() for r in residues[sox2_len + linker_len:] if "CA" in r]
        if not sox2_half or not sox9_half:
            return None
        c1 = np.mean(sox2_half, axis=0)
        c2 = np.mean(sox9_half, axis=0)
        return float(np.linalg.norm(c1 - c2))
    except Exception:
        return None


def screen_pair(pair: dict):
    name = pair["name"]
    sox2_seq = pair["sox2_seq"]
    sox9_seq = pair["sox9_seq"]

    for linker in LINKERS:
        for orient, (first, second, first_len, second_len) in [
            ("fwd", (sox2_seq, sox9_seq, len(sox2_seq), len(sox9_seq))),
            ("rev", (sox9_seq, sox2_seq, len(sox9_seq), len(sox2_seq))),
        ]:
            fused = first + linker + second
            tag = f"{name}_{orient}_{linker}"
            run_dir = OUT_DIR / tag
            run_dir.mkdir(parents=True, exist_ok=True)
            yaml_path = run_dir / "input.yaml"

            write_yaml(yaml_path, fused)
            print(f"\n--- {tag} (len={len(fused)}) ---")

            conf, pdb_path = run_boltz(yaml_path, run_dir)
            if conf is None:
                print("  FAILED")
                results.append({"tag": tag, "status": "FAILED"})
                continue

            iptm   = conf.get("iptm", 0)
            plddt  = conf.get("plddt", 0)
            ptm    = conf.get("ptm", 0)

            sox2_contact = chain_contact_fraction(pdb_path, "C", "A") if pdb_path else 0
            sox9_contact = chain_contact_fraction(pdb_path, "C", "B") if pdb_path else 0

            if orient == "fwd":
                half_dist = midpoint_split_dist(pdb_path, first_len, second_len, len(linker)) if pdb_path else None
            else:
                half_dist = midpoint_split_dist(pdb_path, first_len, second_len, len(linker)) if pdb_path else None

            print(f"  ipTM={iptm:.3f}  pLDDT={plddt:.3f}  ptm={ptm:.3f}")
            print(f"  SOX2 contact frac={sox2_contact:.2f}  SOX9 contact frac={sox9_contact:.2f}")
            print(f"  Half-domain dist={half_dist:.1f} Å" if half_dist else "  Half-domain dist=N/A")

            results.append({
                "tag": tag, "pair": name, "orient": orient, "linker": linker,
                "fused_len": len(fused),
                "iptm": round(iptm, 4), "plddt": round(plddt, 4), "ptm": round(ptm, 4),
                "sox2_contact": round(sox2_contact, 3),
                "sox9_contact": round(sox9_contact, 3),
                "half_dist_A": round(half_dist, 2) if half_dist else None,
                "status": "OK",
            })


if __name__ == "__main__":
    import csv
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for pair in PAIRS:
        screen_pair(pair)

    out_tsv = OUT_DIR / "linker_screen_results.tsv"
    if results:
        keys = list(results[0].keys())
        with open(out_tsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(results)
        print(f"\n=== Results written to {out_tsv} ===")

        ok = [r for r in results if r.get("status") == "OK"]
        if ok:
            ok.sort(key=lambda x: (float(x.get("sox2_contact", 0)) + float(x.get("sox9_contact", 0))), reverse=True)
            print("\nTop hits by dual-contact:")
            for r in ok[:5]:
                print(f"  {r['tag']}  ipTM={r['iptm']}  SOX2={r['sox2_contact']}  SOX9={r['sox9_contact']}  dist={r.get('half_dist_A')}")

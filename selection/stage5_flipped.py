#!/usr/bin/env python3
"""
Stage 5 — flipped orientation test.
Root cause (Option 4 diagnosis): SOX9 binder C-terminus points AWAY from SOX2 binder.
Fix (Option 3): flip chain D to N-terminal fusion.

  Original: chain D = BinderSOX9 + GGS + CaspC  (CaspC at far end, pointing away)
  Flipped:  chain D = CaspC + GGS + BinderSOX9   (CaspC at near end, directly adjacent to CaspN)

Reconstitution distance is now simply: C-term(chain C) <-> N-term(chain D) = standard measurement.

Run: top 4 pairs × GGS linker ×1,2,3 × 3 seeds = 36 runs.
Pass: recon dist < 15 Å in ≥2/3 seeds.
"""

from selection_support import require_execution, reference_sequence, validate_recorded_candidates, campaign_output, campaign_table
require_execution(__name__)

import os, json, subprocess
import numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser

SOX2_TARGET = "DRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPG"
SOX9_TARGET = "PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

BASE      = Path(os.environ.get("GBM_ROOT", Path(__file__).resolve().parents[1]))
OUT_DIR   = campaign_output(BASE / "selection" / "stage5_flipped_outputs")
BOLTZ_BIN = Path(os.environ.get("BOLTZ_BIN", "boltz"))

PAIRS = [
    ("SOX2_mini_l71_s878170_mpnn18_model2",
     "PYEIPEEDREWDEYSFRMVEVHERLMEARRELRLLAPEETEKKRELEEEIARAEKELEEVHKEGEAIERKM",
     "ph_sox9_run18",
     "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"),
    ("SOX2_mini_l71_s878170_mpnn17_model2",
     "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM",
     "ph_sox9_run18",
     "MEEVRELLERARELEEKVRELEDDPSKREERERLLEEAERLREEALDLLVEQIREPLLQWMREEGLL"),
    ("SOX2_mini_l71_s878170_mpnn17_model2",
     "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM",
     "ph_sox9_run30",
     "SEEEKLFEEYERLLAERERDGDLPDLEEIESLLSEEQWQELLKKADEILKKEEEEREEKKKKEEEERKKELQKLLEE"),
    ("SOX2_mini_l71_s878170_mpnn17_model2",
     "PYEIPEEDKEWDEYSFAMVEVHERLMEAERQLRLLSPEEEEKKRELEEEIERARKELEEVHKRGEAIERKM",
     "ph_sox9_run1",
     "PPLSPEEQELEEIEEGYRNLTEKSGEALLKALGITSKNEEEREKKMRELALELAREVTRKLSEKTGVRYAEDGTRLP"),
]

LINKER_NS = [1, 2, 3]
SEEDS     = [0, 1, 2]

pdb_parser = PDBParser(QUIET=True)


def get_casp_halves() -> tuple[str, str]:
    cache = Path(__file__).resolve().parents[1] / "inputs/reference/P55211.fasta"
    seq = "".join(cache.read_text().strip().split("\n")[1:])
    return seq[139:315], seq[315:416]


def recon_dist_flipped(pdb_path: Path) -> float | None:
    """With flipped chain D, CaspC is at the N-terminus of D.
    Reconstitution = C-term(chain C) <-> N-term(chain D)."""
    try:
        s = pdb_parser.get_structure("x", pdb_path)
        m = s[0]
        c_res = list(m["C"].get_residues())
        d_res = list(m["D"].get_residues())
        ca_c = c_res[-1]["CA"].get_vector().get_array()
        ca_d = d_res[0]["CA"].get_vector().get_array()
        return float(np.linalg.norm(ca_c - ca_d))
    except Exception as e:
        print(f"  WARN: {e}")
        return None


def run_boltz(yaml_path: Path, out_dir: Path, seed: int) -> tuple[float | None, float | None]:
    cmd = [str(BOLTZ_BIN), "predict", str(yaml_path),
           "--out_dir", str(out_dir),
           "--recycling_steps", "3", "--sampling_steps", "200",
           "--output_format", "pdb", "--no_kernels", "--seed", str(seed)]
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/usr/lib/wsl/lib:" + env.get("LD_LIBRARY_PATH", "")
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(f"  ERROR: {r.stderr[-300:]}")
        return None, None
    conf_files = list(out_dir.rglob("confidence_*.json"))
    pdb_files  = list(out_dir.rglob("*.pdb"))
    if not conf_files: return None, None
    with open(conf_files[0]) as f: c = json.load(f)
    iptm = c.get("iptm", c.get("protein_iptm"))
    dist = recon_dist_flipped(pdb_files[0]) if pdb_files else None
    return iptm, dist


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    casp_n, casp_c = get_casp_halves()
    print(f"=== Stage 5 flipped orientation ===")
    print(f"Chain C: BinderSOX2 + GGS(n) + CaspN  [unchanged]")
    print(f"Chain D: CaspC + GGS(n) + BinderSOX9  [FLIPPED — CaspC at N-terminus]")

    rows = []
    total = len(PAIRS) * len(LINKER_NS) * len(SEEDS)
    done = 0

    for s2_name, s2_seq, s9_name, s9_seq in PAIRS:
        for n in LINKER_NS:
            linker = "GGS" * n
            chain_c = s2_seq + linker + casp_n          # unchanged
            chain_d = casp_c + linker + s9_seq          # FLIPPED

            for seed in SEEDS:
                done += 1
                name = f"flip_{s2_name[-10:]}_{s9_name}_GGS{n}_s{seed}"
                run_dir = OUT_DIR / name
                run_dir.mkdir(exist_ok=True)

                print(f"\n[{done}/{total}] {s2_name[-20:]} x {s9_name} | GGS×{n} | seed {seed}")

                cached_conf = list(run_dir.rglob("confidence_*.json"))
                if cached_conf:
                    cached_pdb = list(run_dir.rglob("*.pdb"))
                    with open(cached_conf[0]) as f: c = json.load(f)
                    iptm = c.get("iptm", c.get("protein_iptm"))
                    dist = recon_dist_flipped(cached_pdb[0]) if cached_pdb else None
                    print(f"  (cached) recon={dist:.1f}Å  iptm={iptm:.3f}" if dist else "  cached")
                    rows.append({"sox2": s2_name, "sox9": s9_name, "linker_n": n, "seed": seed,
                                 "recon_dist_A": dist, "iptm": iptm, "pass": dist is not None and dist < 15.0})
                    continue

                yaml = f"""sequences:
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
      sequence: {chain_c}
      msa: empty
  - protein:
      id: D
      sequence: {chain_d}
      msa: empty
"""
                (run_dir / "input.yaml").write_text(yaml)
                iptm, dist = run_boltz(run_dir / "input.yaml", run_dir, seed)
                passed = dist is not None and dist < 15.0
                print(f"  recon={dist:.1f}Å  iptm={iptm:.3f}  pass={passed}" if dist else "  FAILED")
                rows.append({"sox2": s2_name, "sox9": s9_name, "linker_n": n, "seed": seed,
                             "recon_dist_A": dist, "iptm": iptm, "pass": passed})

    df = pd.DataFrame(rows)
    out_path = campaign_table(BASE / "selection" / "stage5_flipped_results.tsv")
    df.to_csv(out_path, sep="\t", index=False)

    print(f"\n=== Flipped orientation results ===")
    summary = df.groupby(["sox2","sox9","linker_n"]).agg(
        mean_dist=("recon_dist_A","mean"), min_dist=("recon_dist_A","min"),
        n_pass=("pass","sum"), n_seeds=("pass","count")
    ).reset_index().sort_values("min_dist")
    print(summary.to_string(index=False))
    print(f"\nTotal passing runs: {df['pass'].sum()}/{len(df)}")


if __name__ == "__main__":
    main()

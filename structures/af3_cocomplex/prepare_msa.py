#!/usr/bin/env python3
"""Convert Boltz-2 CSV MSAs to A3M format for AlphaFold3."""
import csv
import json
from pathlib import Path

MSA_DIR = Path("/mnt/d/gbm-protein-switch/structures/boltz2_cocomplex/output/boltz_results_sox2_sox9_cocomplex/msa")
OUT_DIR = Path("/mnt/d/gbm-protein-switch/structures/af3_cocomplex/msa")
OUT_DIR.mkdir(exist_ok=True)

def csv_to_a3m(csv_path: Path, query_name: str) -> Path:
    a3m_path = OUT_DIR / (csv_path.stem + ".a3m")
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    with open(a3m_path, "w") as out:
        for i, row in enumerate(rows):
            seq = row["sequence"]
            if i == 0:
                out.write(f">{query_name}\n{seq}\n")
            else:
                out.write(f">seq_{i}\n{seq}\n")
    print(f"Wrote {len(rows)} sequences to {a3m_path}")
    return a3m_path

sox2_a3m = csv_to_a3m(MSA_DIR / "sox2_sox9_cocomplex_0.csv", "SOX2_P48431")
sox9_a3m = csv_to_a3m(MSA_DIR / "sox2_sox9_cocomplex_1.csv", "SOX9_P48436")

# Update AF3 JSON with MSA paths and empty templates
json_path = Path("/mnt/d/gbm-protein-switch/structures/af3_cocomplex/sox2_sox9_cocomplex_af3.json")
with open(json_path) as f:
    data = json.load(f)

for entry in data["sequences"]:
    if "protein" in entry:
        chain_id = entry["protein"]["id"]
        if chain_id == "A":
            entry["protein"]["unpairedMsaPath"] = str(sox2_a3m)
            entry["protein"]["pairedMsa"] = ""
            entry["protein"]["templates"] = []
        elif chain_id == "B":
            entry["protein"]["unpairedMsaPath"] = str(sox9_a3m)
            entry["protein"]["pairedMsa"] = ""
            entry["protein"]["templates"] = []

out_json = Path("/mnt/d/gbm-protein-switch/structures/af3_cocomplex/sox2_sox9_cocomplex_af3_withmsa.json")
with open(out_json, "w") as f:
    json.dump(data, f, indent=2)
print(f"Written updated JSON to {out_json}")

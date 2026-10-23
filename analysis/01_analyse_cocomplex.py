#!/usr/bin/env python3
"""
Analyse Boltz-2 SOX2:SOX9:DNA co-complex model_0.
Outputs:
  - DNA-contacting residues on SOX2 and SOX9
  - SOX2-SOX9 protein-protein interface residues
  - Exposed epitopes (good targets for miniprotein binder design)
  - PyMOL colouring script
"""
import numpy as np
from pathlib import Path
from Bio import PDB
from Bio.PDB import SASA

PDB_PATH = Path("/mnt/d/gbm-protein-switch/structures/boltz2_cocomplex/output/"
                "boltz_results_sox2_sox9_cocomplex/predictions/sox2_sox9_cocomplex/"
                "sox2_sox9_cocomplex_model_0.pdb")
PLDDT_PATH = Path("/mnt/d/gbm-protein-switch/structures/boltz2_cocomplex/output/"
                  "boltz_results_sox2_sox9_cocomplex/predictions/sox2_sox9_cocomplex/"
                  "plddt_sox2_sox9_cocomplex_model_0.npz")
OUT_DIR = Path("/mnt/d/gbm-protein-switch/analysis")
OUT_DIR.mkdir(exist_ok=True)

DNA_CONTACT_CUTOFF  = 4.5   # Å — residue atom to DNA atom
PPI_CONTACT_CUTOFF  = 5.0   # Å — SOX2 atom to SOX9 atom
SASA_EXPOSED_CUTOFF = 15.0  # Å² — residue considered surface-exposed
PLDDT_CUTOFF        = 50.0  # only report confident residues

# ── Load structure ────────────────────────────────────────────────────────────
parser = PDB.PDBParser(QUIET=True)
struct = parser.get_structure("sox_complex", PDB_PATH)
model  = struct[0]

chain_sox2 = model["A"]   # SOX2 (317 aa)
chain_sox9 = model["B"]   # SOX9 (509 aa)
chain_dna1 = model["C"]   # forward strand
chain_dna2 = model["D"]   # reverse strand

# ── pLDDT (stored in B-factor) ────────────────────────────────────────────────
def get_residue_plddt(chain):
    plddts = {}
    for res in chain.get_residues():
        if res.id[0] != " ":
            continue
        bvals = [a.get_bfactor() for a in res.get_atoms()]
        plddts[res.id[1]] = np.mean(bvals)
    return plddts

sox2_plddt = get_residue_plddt(chain_sox2)
sox9_plddt = get_residue_plddt(chain_sox9)

# ── SASA ──────────────────────────────────────────────────────────────────────
sr = SASA.ShrakeRupley()
sr.compute(struct, level="R")

def get_sasa(chain):
    return {res.id[1]: res.sasa for res in chain.get_residues() if res.id[0] == " "}

sox2_sasa = get_sasa(chain_sox2)
sox9_sasa = get_sasa(chain_sox9)

# ── Helper: atom list for a chain ────────────────────────────────────────────
def atoms(chain):
    return [a for r in chain.get_residues() if r.id[0]==" " for a in r.get_atoms()]

def atoms_of_chain_list(chains):
    return [a for ch in chains for a in atoms(ch)]

# ── DNA-contacting residues ───────────────────────────────────────────────────
dna_atoms = atoms_of_chain_list([chain_dna1, chain_dna2])

def find_dna_contacts(chain, cutoff=DNA_CONTACT_CUTOFF):
    contacts = set()
    for res in chain.get_residues():
        if res.id[0] != " ":
            continue
        for a in res.get_atoms():
            for da in dna_atoms:
                if a - da < cutoff:
                    contacts.add(res.id[1])
                    break
            if res.id[1] in contacts:
                break
    return contacts

print("Computing DNA contacts...")
sox2_dna = find_dna_contacts(chain_sox2)
sox9_dna = find_dna_contacts(chain_sox9)

# ── SOX2-SOX9 PPI interface ───────────────────────────────────────────────────
sox9_atoms = atoms(chain_sox9)
sox2_atoms = atoms(chain_sox2)

def find_ppi_contacts(chain_a_residues, other_atoms, cutoff=PPI_CONTACT_CUTOFF):
    contacts = set()
    for res in chain_a_residues.get_residues():
        if res.id[0] != " ":
            continue
        for a in res.get_atoms():
            for oa in other_atoms:
                if a - oa < cutoff:
                    contacts.add(res.id[1])
                    break
            if res.id[1] in contacts:
                break
    return contacts

print("Computing SOX2-SOX9 PPI interface...")
sox2_ppi = find_ppi_contacts(chain_sox2, sox9_atoms)
sox9_ppi = find_ppi_contacts(chain_sox9, sox2_atoms)

# ── Exposed epitopes (surface, not DNA/PPI, confident) ───────────────────────
def find_epitopes(chain, sasa_dict, plddt_dict, dna_contacts, ppi_contacts):
    epitopes = []
    for res in chain.get_residues():
        if res.id[0] != " ":
            continue
        rnum = res.id[1]
        sasa  = sasa_dict.get(rnum, 0)
        pldd  = plddt_dict.get(rnum, 0)
        if (sasa  >= SASA_EXPOSED_CUTOFF and
            pldd  >= PLDDT_CUTOFF and
            rnum not in dna_contacts and
            rnum not in ppi_contacts):
            epitopes.append({
                "resnum": rnum,
                "resname": res.resname,
                "sasa": round(sasa, 1),
                "plddt": round(pldd, 1),
            })
    return epitopes

sox2_epitopes = find_epitopes(chain_sox2, sox2_sasa, sox2_plddt, sox2_dna, sox2_ppi)
sox9_epitopes = find_epitopes(chain_sox9, sox9_sasa, sox9_plddt, sox9_dna, sox9_ppi)

# ── Summarise ─────────────────────────────────────────────────────────────────
def rnum_list(s): return sorted(s)

print("\n" + "="*60)
print("SOX2 (chain A) — 317 aa")
print(f"  DNA contacts ({DNA_CONTACT_CUTOFF}Å): {sorted(sox2_dna)}")
print(f"  PPI contacts ({PPI_CONTACT_CUTOFF}Å): {sorted(sox2_ppi)}")
print(f"  Exposed epitopes (pLDDT≥{PLDDT_CUTOFF}, SASA≥{SASA_EXPOSED_CUTOFF}): "
      f"{[e['resnum'] for e in sox2_epitopes]}")

print("\nSOX9 (chain B) — 509 aa")
print(f"  DNA contacts ({DNA_CONTACT_CUTOFF}Å): {sorted(sox9_dna)}")
print(f"  PPI contacts ({PPI_CONTACT_CUTOFF}Å): {sorted(sox9_ppi)}")
print(f"  Exposed epitopes (pLDDT≥{PLDDT_CUTOFF}, SASA≥{SASA_EXPOSED_CUTOFF}): "
      f"{[e['resnum'] for e in sox9_epitopes]}")

# ── Per-residue table for top epitope residues ────────────────────────────────
print("\n── SOX2 top epitope residues (sorted by SASA) ──")
for e in sorted(sox2_epitopes, key=lambda x: -x["sasa"])[:30]:
    print(f"  A:{e['resnum']:4d} {e['resname']}  SASA={e['sasa']:6.1f}  pLDDT={e['plddt']:.0f}")

print("\n── SOX9 top epitope residues (sorted by SASA) ──")
for e in sorted(sox9_epitopes, key=lambda x: -x["sasa"])[:30]:
    print(f"  B:{e['resnum']:4d} {e['resname']}  SASA={e['sasa']:6.1f}  pLDDT={e['plddt']:.0f}")

# ── Identify contiguous epitope patches (for BindCraft hotspot selection) ─────
def find_patches(epitopes, gap=8):
    """Group consecutive residues (within gap) into patches."""
    if not epitopes:
        return []
    nums = sorted(e["resnum"] for e in epitopes)
    patches, cur = [[nums[0]]], [nums[0]]
    for n in nums[1:]:
        if n - cur[-1] <= gap:
            cur.append(n)
        else:
            patches.append(cur)
            cur = [n]
    patches.append(cur)
    return patches

sox2_patches = find_patches(sox2_epitopes)
sox9_patches = find_patches(sox9_epitopes)

print("\n── SOX2 epitope patches ──")
for p in sox2_patches:
    avg_sasa  = np.mean([e["sasa"]  for e in sox2_epitopes if e["resnum"] in p])
    avg_plddt = np.mean([e["plddt"] for e in sox2_epitopes if e["resnum"] in p])
    print(f"  Residues {p[0]}-{p[-1]} (n={len(p)})  avg_SASA={avg_sasa:.0f}  avg_pLDDT={avg_plddt:.0f}")

print("\n── SOX9 epitope patches ──")
for p in sox9_patches:
    avg_sasa  = np.mean([e["sasa"]  for e in sox9_epitopes if e["resnum"] in p])
    avg_plddt = np.mean([e["plddt"] for e in sox9_epitopes if e["resnum"] in p])
    print(f"  Residues {p[0]}-{p[-1]} (n={len(p)})  avg_SASA={avg_sasa:.0f}  avg_pLDDT={avg_plddt:.0f}")

# ── Write PyMOL colouring script ──────────────────────────────────────────────
pml = OUT_DIR / "colour_cocomplex.pml"
with open(pml, "w") as f:
    f.write(f"load {PDB_PATH}, sox_complex\n")
    f.write("hide everything\n")
    f.write("show cartoon, chain A+B\n")
    f.write("show sticks, chain C+D\n")
    f.write("color grey80, chain A\n")
    f.write("color grey60, chain B\n")
    f.write("color yellow, chain C+D\n")
    # DNA contacts
    if sox2_dna:
        sel = "+".join(str(r) for r in sorted(sox2_dna))
        f.write(f"select sox2_dna, chain A and resi {sel}\n")
        f.write("color orange, sox2_dna\n")
    if sox9_dna:
        sel = "+".join(str(r) for r in sorted(sox9_dna))
        f.write(f"select sox9_dna, chain B and resi {sel}\n")
        f.write("color orange, sox9_dna\n")
    # PPI interface
    if sox2_ppi:
        sel = "+".join(str(r) for r in sorted(sox2_ppi))
        f.write(f"select sox2_ppi, chain A and resi {sel}\n")
        f.write("color red, sox2_ppi\n")
    if sox9_ppi:
        sel = "+".join(str(r) for r in sorted(sox9_ppi))
        f.write(f"select sox9_ppi, chain B and resi {sel}\n")
        f.write("color red, sox9_ppi\n")
    # Epitopes
    if sox2_epitopes:
        sel = "+".join(str(e["resnum"]) for e in sox2_epitopes)
        f.write(f"select sox2_epitopes, chain A and resi {sel}\n")
        f.write("color cyan, sox2_epitopes\n")
    if sox9_epitopes:
        sel = "+".join(str(e["resnum"]) for e in sox9_epitopes)
        f.write(f"select sox9_epitopes, chain B and resi {sel}\n")
        f.write("color cyan, sox9_epitopes\n")
    f.write("zoom\n")
    f.write("# Colour key:\n")
    f.write("# grey   = protein backbone\n")
    f.write("# yellow = DNA\n")
    f.write("# orange = DNA-contacting residues\n")
    f.write("# red    = SOX2-SOX9 interface residues\n")
    f.write("# cyan   = exposed epitopes (binder targets)\n")

print(f"\nPyMOL script written to {pml}")

# ── Save epitope residue lists for BindCraft ──────────────────────────────────
import json
out = {
    "model": "sox2_sox9_cocomplex_model_0",
    "sox2_dna_contacts":  sorted(sox2_dna),
    "sox2_ppi_contacts":  sorted(sox2_ppi),
    "sox2_epitopes": [{"chain": "A", **e} for e in sox2_epitopes],
    "sox2_patches": sox2_patches,
    "sox9_dna_contacts":  sorted(sox9_dna),
    "sox9_ppi_contacts":  sorted(sox9_ppi),
    "sox9_epitopes": [{"chain": "B", **e} for e in sox9_epitopes],
    "sox9_patches": sox9_patches,
}
with open(OUT_DIR / "epitope_analysis.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"Epitope data saved to {OUT_DIR}/epitope_analysis.json")

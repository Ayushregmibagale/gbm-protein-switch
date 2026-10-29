# Caspase Geometry Modelling Plan

## Goal
Verify that the binder-linker-caspase fusion constructs adopt a geometry where:
- Caspase half is not buried against the SOX2/SOX9 surface
- For Option 1 (split caspase): the two halves are within ~15 Å when both binders are co-bound
- For Option 2 (interface binder + full caspase): caspase is accessible and not sterically clashing

## Why not pRosettaC
pRosettaC is built specifically for PROTAC ternary complexes (small molecule bridge). This system is all protein-protein. Use AF3/Boltz-2 + RosettaDock + GROMACS instead.

## Workflow

### Step 1 — Full assembly prediction (AF3 or Boltz-2)
Model entire assembly as separate chains with linker encoded in sequence:

**Option 1 (split caspase):**
```
Chain A: SOX2 HMG (6WX8 chain B)
Chain B: SOX9 HMG (4EUW chain A)  
Chain C: BinderSOX2 + GGSGGSGGSGGS + Caspase9_N-half (aa 140-315)
Chain D: BinderSOX9 + GGSGGSGGSGGS + Caspase9_C-half (aa 316-416)
```

**Option 2 (interface binder + full caspase):**
```
Chain A: SOX2 HMG
Chain B: SOX9 HMG
Chain C: InterfaceBinder + GGSGGSGGSGGS + Caspase9_full (CARD-deleted, aa 140-416)
```

Caspase-9 sequences (CARD-deleted):
- N-half (aa 140-315): from UniProt P55211
- C-half (aa 316-416): from UniProt P55211
- Full CARD-deleted (aa 140-416): from UniProt P55211

### Step 2 — Check geometry from model
Two criteria:
1. **Caspase not buried**: solvent-accessible surface of caspase catalytic residues (C287, H237)
2. **Split reconstitution distance**: C-terminus of CaspA_half ↔ N-terminus of CaspB_half < 15 Å

### Step 3 — Linker length sweep
Run AF3/Boltz-2 with GGS × 1, 2, 3, 4 repeats (5, 10, 15, 20 aa linkers).
Geometry is sensitive to linker length — need to find the right window.

### Step 4 — GROMACS MD
20–50 ns on best pose per linker length.
Monitor:
- Caspase half RMSD (stable fold?)
- Inter-half distance (does it converge toward active geometry?)
- Binder-target contact persistence

## Tools available
- Boltz-2: `boltz2` conda env, `~/miniconda3/envs/boltz2/`
- AF3: `alphafold3` conda env, weights at `~/alphafold3/weights/af3.bin.zst`
- GROMACS: full PLUMED+GROMACS 2024.3 build
- RosettaDock: via PyRosetta (protac env)

## Trigger
Start AFTER binder design campaigns finish and post-design ColabFold quaternary filter is done.
Priority order: quaternary ColabFold filter first → linker sweep → MD

## Status
PLANNED — not started

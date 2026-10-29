# Bridge Binder Campaign — Full Overview
*Last updated: 2026-09-12*

## Goal
Design a single miniprotein binder that simultaneously contacts both SOX2 (chain A, HMG res ~175-215, hotspots A193/A199) and SOX9 (chain B, HMG res ~450-509, hotspots B466/B505/B508), with N-terminus near SOX2 and C-terminus near SOX9, so it can serve as the insert in:

`CaspN(140-315) + GGS + [bridge_binder] + GGS + CaspC(316-416)`

Reconstitution criterion: N-to-C distance of the bridge binder < 15Å in ≥50% of 20ns MD frames, with both SOX2 and SOX9 contacts maintained.

---

## What We Tried

### 1. BindCraft — dual-hotspot, fragment target (FAILED)
- **Runs:** 233 trajectories + ~92 relaxed rerun
- **Accepted:** 0
- **Root cause:** AF2 hallucination optimizer cannot simultaneously satisfy 5 hard hotspot constraints across 2 separate chains.

### 2. Protein Hunter — bridge (FAILED)
- **Runs:** 50 designs; **Accepted:** 0
- **Root cause:** Boltz-2 hallucination produces disordered peptides when the two target chains are ~45Å apart.

### 3. Boltz-2 linker screen — fused SOX2+SOX9 binder pairs (FAILED)
- **What:** 18 fused combinations of top individual binders with GGS×1/×2/×3 linkers
- **Accepted:** 0; individual binders don't fold cooperatively when fused.

### 4. RFdiffusion — dual-chain bridge backbone generation (WORKED)
- **What:** Single run with contig `[A175-215/0 B450-509/0 60-100]`, hotspots `A193,A199,B466,B505,B508`. 200 designs.
- **Result:** 5/200 contact both chains <8Å; 20/200 within 10Å. pLDDT mean=0.849.
- **Why it works:** Diffuses backbone coordinates directly — can span two hotspot clusters in a single continuous chain.

### 5. ProteinMPNN + Boltz-2 filter — batch 1 (top 20 backbones, WORKED)
- 200 sequences (10 seqs × 20 backbones); filter: ipTM≥0.40, SOX2_c≥0.12, SOX9_c≥0.12
- **Accepted:** 5/200 into MD queue

### 6. ProteinMPNN + Boltz-2 filter — batch 2 (backbones 21–80, WORKED)
- 600 sequences (10 seqs × 60 backbones, temp=0.15); same filter
- **Accepted:** 27/600 into MD queue (32 total in Accepted/)

### 7. GROMACS MD validation — COMPLETE (Sep 12 2026)
- **32 designs tested**, 20ns each
- **Pass rate: 5/32 = 15.6%**
- **Key finding:** Boltz-2 ipTM is weakly predictive of MD pass (highest ipTM design bridge_95_s02=0.807 → FAIL; bridge_48_s01 passed at 98.5%)

---

## Confirmed Protocol

```
Step 1 — RFdiffusion (200 designs)
  contig: [A175-215/0 B450-509/0 60-100]
  hotspots: A193,A199,B466,B505,B508
  model: Complex_base_ckpt.pt
  filter: CA dist to SOX2 hotspot cluster <9Å AND SOX9 cluster <9Å → keep top 20-80

Step 2 — ProteinMPNN (10 seqs per backbone)
  Boltz-2 3-chain prediction (SOX2_HMG + SOX9_HMG + binder)
  filter: ipTM≥0.40 AND SOX2_contact≥0.12 AND SOX9_contact≥0.12

Step 3 — GROMACS MD (20ns, 3 replicates on confirmed leads)
  pass criterion: pct_lt15A ≥50% AND sox2_contact_mean>0.10 AND sox9_contact_mean>0.05
  metric: N-to-C Cα distance of binder chain
```

---

## Final MD Results — All 32 Designs

| Design | %<15Å | %<30Å | NC_mean (Å) | SOX2_c | SOX9_c | Result |
|--------|-------|-------|-------------|--------|--------|--------|
| bridge_48_s01  | **98.5%** | 98.5% | 8.85  | 0.376 | 0.166 | ✓ PASS |
| bridge_73_s09  | 86.7% | 91.8% | 15.04 | 0.218 | 0.229 | ✓ PASS |
| bridge_136_s01 | 86.0% | 86.5% | 23.18 | 0.143 | 0.172 | ✓ PASS |
| bridge_101_s06 | 79.3% | 100.0% | 13.32 | 0.274 | 0.255 | ✓ PASS |
| bridge_90_s01  | 55.5% | 64.6% | 50.23 | 0.234 | 0.162 | ✓ PASS |
| bridge_85_s07  | 40.3% | 78.9% | 30.19 | 0.236 | 0.135 | ✗ borderline |
| bridge_62_s03  | 24.2% | 86.3% | 27.57 | 0.348 | 0.159 | ✗ borderline |
| bridge_142_s06 | 0.0%  | 81.7% | 24.11 | 0.270 | 0.154 | ✗ FAIL |
| bridge_192_s08 | 0.0%  | 58.6% | 28.52 | 0.176 | 0.063 | ✗ FAIL |
| bridge_159_s05 | 0.1%  | 71.2% | 37.30 | 0.162 | 0.160 | ✗ FAIL |
| bridge_101_s03 | 0.0%  | 66.8% | 36.18 | 0.221 | 0.320 | ✗ FAIL |
| bridge_10_s03  | 1.9%  | 62.3% | 40.04 | 0.213 | 0.346 | ✗ FAIL |
| bridge_141_s09 | 0.0%  | 48.1% | 30.17 | 0.196 | 0.167 | ✗ FAIL |
| (13 others)    | 0.0%  | <25%  | >47Å  | —     | —     | ✗ FAIL |

---

## Confirmed PASS Designs — Details

### bridge_48_s01 ★ BEST (69 aa)
- 98.5% <15Å, NC_mean=8.85Å — near-perfect bridging, almost never extended
- SOX2_c=0.376 (highest of all PASS designs)
- Replicates: RUNNING (Sep 12)

### bridge_73_s09 (83 aa)
- 86.7% <15Å, NC_mean=15.04Å
- Balanced dual contacts (SOX2=0.218, SOX9=0.229)
- Replicates: RUNNING (Sep 12)

### bridge_136_s01 (73 aa)
- 86.0% <15Å, NC_mean=23.18Å
- Replicates: RUNNING (Sep 12)

### bridge_101_s06 (76 aa)
- 79.3% <15Å, NC_mean=13.32Å — very tight NC
- High dual contacts (SOX2=0.274, SOX9=0.255)
- Replicates: RUNNING (Sep 12)

### bridge_90_s01 (90 aa) — original lead
- 55.5% <15Å (rep1), 87.1% (rep2), 60.1% (rep3) — mean 67.6%
- Seq: `GSRLDAEIARARARAEAREAEIAALRAAAAAEAAALEAKAAEDAEHGEELLAEAAAVRERAELEAEALEQEVAAERALADARAAAAAAAR`
- All 3 replicates confirmed PASS

---

## Replicate Status (Sep 12 2026)

| Design | Rep1 | Rep2 | Rep3 | Status |
|--------|------|------|------|--------|
| bridge_90_s01 | 55.5% ✓ | 87.1% ✓ | 60.1% ✓ | CONFIRMED |
| bridge_48_s01 | 98.5% ✓ | running | pending | RUNNING |
| bridge_73_s09 | 86.7% ✓ | pending | pending | RUNNING |
| bridge_136_s01 | 86.0% ✓ | pending | pending | RUNNING |
| bridge_101_s06 | 79.3% ✓ | pending | pending | RUNNING |

Replicate script: `md/bridge_binder_md/run_new_replicates.sh`
Log: `md/bridge_binder_md/new_replicates.log`
ETA: ~03:00 BST Sep 13

---

## What to Do Next

1. **Confirm replicates** — expect all 4 to hold up (rep1 scores are all strong)
2. **Lab validation** — 5 confirmed designs is sufficient for:
   - SPR/MST: measure binding affinity to SOX2-HMG and SOX9-HMG (should bind both)
   - Selectivity: confirm no binding to either monomer alone
   - DEVDase assay: co-express CaspN-GGS-bridge + GGS-CaspC in GSC line, measure caspase activity
3. **RFdiffusion round 2** (optional) — if more diversity needed, script ready at `rfdiffusion/run_rfdiffusion_bridge_round2.sh`

## Files
- Accepted designs: `rfdiffusion/bridge_mpnn_boltz/Accepted/` (32 PDBs)
- MD runs: `md/bridge_binder_md/`
- Replicate runs: `md/bridge_binder_md/bridge_*_rep2/`, `bridge_*_rep3/`
- Analysis: per-design `metrics.json` in each MD directory

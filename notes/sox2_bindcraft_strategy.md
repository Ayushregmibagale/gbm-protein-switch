# SOX2 BindCraft Strategy & Filter Rationale

## Problem
SOX2 acceptance rate: ~0.5% (7/1430 trajectories) vs SOX9 v3: ~45% (40/88 trajectories).
Root cause: SOX2 HMG wing-2 face (H2-H3 helices, residues 62-81) is flat — no groove or pocket.
Flat surfaces score poorly on ShapeComplementarity even when binders are geometrically valid.

## SOX2 crash (Aug 29)
DSSP segfaults at `ffffffffff600400` (WSL2 vsyscall page incompatibility) propagated into Python.
Fix: watchdog wrapper scripts auto-restart on crash until 20 accepted designs found.
Scripts: `run_sox2_watchdog.sh`, `run_sox9_watchdog.sh`

## Filter relaxation (Campaign 1 — current)
Created `~/projects/BindCraft/settings_filters/sox2_relaxed_filters.json`:

| Filter                | Default | Relaxed |
|-----------------------|---------|---------|
| ShapeComplementarity (avg) | ≥ 0.60 | ≥ 0.45 |
| ShapeComplementarity (1/2) | ≥ 0.55 | ≥ 0.40 |
| n_InterfaceHbonds (avg/1/2) | ≥ 3   | ≥ 2    |
| i_pAE (avg/1/2)       | ≤ 0.35  | ≤ 0.42 |

Rationale: flat surfaces intrinsically score ~0.4-0.5 on SC; ColabFold re-validates all
designs anyway so relaxing here doesn't reduce final quality.

Hotspots unchanged: **B62, B70, B74, B81** (wing-2 face, confirmed correct AND-gate geometry)
Design dir: `designs/sox2_binder_xtal/` (7 accepted already kept)

## Groove binding site (Campaign 2 — ready, not launched)
H1/L1/H2 junction concavity (auth ~59-64) identified as shallow groove on same face as
SOX9 hotspots. Better ShapeComplementarity expected vs flat helix face.

Groove hotspot geometry (from AF3 co-complex, best model seed-3_sample-3):

| Auth | AA  | Angle to SOX9 epi | DNA dist | Notes |
|------|-----|-------------------|----------|-------|
| 59   | Met | 84.5°             | 8.1 Å    | groove entrance |
| 61   | Gln | 87.0°             | 6.9 Å    | groove |
| 62   | Glu | 80.9°             | 9.9 Å    | shared with campaign 1 |
| 63   | Asn | 70.8°             | 9.0 Å    | groove |
| 70   | Glu | 54.6°             | 6.1 Å    | shared with campaign 1 |

All angles < 90° → same face as SOX9 hotspots → AND-gate geometry preserved.
None at SOX2-SOX9 PPI interface. All DNA dist > 6 Å.

Settings: `~/projects/BindCraft/settings_target/SOX2_groove_miniprotein.json`
Design dir: `designs/sox2_binder_groove/`
Filter: `sox2_relaxed_filters.json`
Script: `run_sox2_groove_watchdog.sh` (NOT launched — GPU busy with campaign 1)

## GPU constraint
RTX 5090 VRAM: 32.6 GB. One BindCraft run uses ~21 GB → cannot run two simultaneously.
Launch groove campaign AFTER campaign 1 finishes or pipeline (PH/BoltzGen) completes.

## Pipeline launcher
v6 PID 4146779, log: `pipeline_launcher_v6.log`
Sequence: SOX2 (campaign 1) done → PH SOX2 → PH SOX9 → BoltzGen SOX2 → BoltzGen SOX9
SOX9 v3: COMPLETE (40 designs in `sox9_binder_xtal_v3/Accepted/`)

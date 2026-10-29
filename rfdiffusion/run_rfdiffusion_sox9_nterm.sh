#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# RFdiffusion: SOX9 binders with hotspots biased to the N-terminal half of the HMG fragment.
# Goal: produce binders whose N-terminus is on the SOX2-proximal face,
# enabling CaspC + GGS + BinderSOX9 N-terminal caspase fusion.
#
# Target: sox2_sox9_bridge_fragment.pdb, chain B only (SOX9 HMG res 450-509)
# Hotspots: B453,B457,B463,B468 — N-terminal half of SOX9 fragment (res 450-470 face)
# Contig: B450-509/0 55-90  (binder length 55-90 aa)
# 200 designs → filter by N-terminal direction test (Stage5 terminus diagnosis)
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

TARGET=/mnt/d/gbm-protein-switch/bindcraft/targets/sox2_sox9_bridge_fragment.pdb
OUTDIR=/mnt/d/gbm-protein-switch/rfdiffusion/sox9_nterm
RFDIR=/home/ayushregmibagale/projects/RFdiffusion
MODEL=$RFDIR/models/Complex_base_ckpt.pt

mkdir -p "$OUTDIR"

echo "=== RFdiffusion SOX9 N-term binders — $(date) ==="
echo "Hotspots biased to B450-B470 face to produce N-terminal-proximal binding mode"

python $RFDIR/scripts/run_inference.py \
    inference.input_pdb=$TARGET \
    inference.output_prefix=$OUTDIR/sox9_nterm \
    inference.model_directory_path=$RFDIR/models \
    inference.ckpt_override_path=$MODEL \
    'contigmap.contigs=[B450-509/0 55-90]' \
    'ppi.hotspot_res=[B453,B457,B463,B468]' \
    inference.num_designs=200

echo "=== RFdiffusion SOX9 N-term done — $(date) ==="
echo "Next: ProteinMPNN → Boltz-2 Stage1 filter → Stage5 terminus direction diagnosis"
echo "Select designs where N-terminus direction vector points TOWARD SOX2 (d_toward > 0)"

#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# RFdiffusion: single binder targeting BOTH SOX2 (chain A) and SOX9 (chain B) simultaneously.
# Target: sox2_sox9_bridge_fragment.pdb
#   Chain A: SOX2 HMG res 175-215 (41 aa)
#   Chain B: SOX9 HMG res 450-509 (60 aa)
# Hotspots: A193,A199 (SOX2) + B466,B505,B508 (SOX9)
# Contig: A175-215/0 B450-509/0 60-100  (binder length 60-100 aa)
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

TARGET=/mnt/d/gbm-protein-switch/bindcraft/targets/sox2_sox9_bridge_fragment.pdb
OUTDIR=/mnt/d/gbm-protein-switch/rfdiffusion/bridge_both_chains
RFDIR=/home/ayushregmibagale/projects/RFdiffusion
MODEL=$RFDIR/models/Complex_base_ckpt.pt

mkdir -p "$OUTDIR"

echo "=== RFdiffusion bridge both chains — $(date) ==="
echo "Target: $TARGET"
echo "Output: $OUTDIR"

python $RFDIR/scripts/run_inference.py \
    inference.input_pdb=$TARGET \
    inference.output_prefix=$OUTDIR/bridge \
    inference.model_directory_path=$RFDIR/models \
    inference.ckpt_override_path=$MODEL \
    'contigmap.contigs=[A175-215/0 B450-509/0 60-100]' \
    'ppi.hotspot_res=[A193,A199,B466,B505,B508]' \
    inference.num_designs=200

echo "=== RFdiffusion bridge both chains done — $(date) ==="
echo "Outputs in: $OUTDIR"
echo "Next: run ProteinMPNN + Boltz-2 filter on the 200 backbones"

#!/bin/bash
# RFdiffusion round 2: 500 more bridge binder backbones.
# Two sub-runs to sample different length regimes:
#   run A: 60-100 aa (same as round 1) — 250 designs
#   run B: 50-80 aa  (shorter/more rigid) — 250 designs
# Outputs to rfdiffusion/bridge_round2/
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH:-}

TARGET=/mnt/d/gbm-protein-switch/bindcraft/targets/sox2_sox9_bridge_fragment.pdb
OUTDIR=/mnt/d/gbm-protein-switch/rfdiffusion/bridge_round2
RFDIR=/home/ayushregmibagale/projects/RFdiffusion
MODEL=$RFDIR/models/Complex_base_ckpt.pt

mkdir -p "$OUTDIR/long" "$OUTDIR/short"

echo "=== RFdiffusion bridge round 2 — $(date) ==="

echo "--- Run A: 60-100 aa (250 designs) ---"
python $RFDIR/scripts/run_inference.py \
    inference.input_pdb=$TARGET \
    inference.output_prefix=$OUTDIR/long/bridge \
    inference.model_directory_path=$RFDIR/models \
    inference.ckpt_override_path=$MODEL \
    'contigmap.contigs=[A175-215/0 B450-509/0 60-100]' \
    'ppi.hotspot_res=[A193,A199,B466,B505,B508]' \
    inference.num_designs=250

echo "--- Run B: 50-80 aa (250 designs) ---"
python $RFDIR/scripts/run_inference.py \
    inference.input_pdb=$TARGET \
    inference.output_prefix=$OUTDIR/short/bridge \
    inference.model_directory_path=$RFDIR/models \
    inference.ckpt_override_path=$MODEL \
    'contigmap.contigs=[A175-215/0 B450-509/0 50-80]' \
    'ppi.hotspot_res=[A193,A199,B466,B505,B508]' \
    inference.num_designs=250

echo "=== RFdiffusion round 2 done — $(date) ==="
echo "Next: run run_mpnn_boltz_bridge_round2.py on these backbones"

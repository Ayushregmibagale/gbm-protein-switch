#!/bin/bash
# RFdiffusion binder design for SOX9 epitope (B:141-172)
# Hotspots: B148, B152, B159, B166, B171 — HMG C-tail patch
# 200 diffusion trajectories
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

TARGET=/mnt/d/gbm-protein-switch/bindcraft/targets/sox2_sox9_cocomplex.pdb
OUTDIR=/mnt/d/gbm-protein-switch/rfdiffusion/sox9_binder
RFDIR=/home/ayushregmibagale/projects/RFdiffusion
MODEL=$RFDIR/models/Complex_base_ckpt.pt

echo "=== RFdiffusion SOX9 binder — $(date) ==="

python $RFDIR/scripts/run_inference.py \
    inference.input_pdb=$TARGET \
    inference.output_prefix=$OUTDIR/sox9_binder \
    inference.model_directory_path=$RFDIR/models \
    inference.ckpt_override_path=$MODEL \
    'contigmap.contigs=[B1-509/0 60-100]' \
    'ppi.hotspot_res=[B148,B152,B159,B166,B171]' \
    inference.num_designs=200 \
    inference.noise_scale_ca=1 \
    inference.noise_scale_frame=1 \
    denoiser.noise_scale_ca=1 \
    denoiser.noise_scale_frame=1

echo "=== RFdiffusion SOX9 done — $(date) ==="

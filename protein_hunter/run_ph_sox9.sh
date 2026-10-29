#!/bin/bash
# Protein Hunter — SOX9 binder design
# Target: SOX9 HMG domain only (72 aa, auth 103-174), hotspots 66,67,69,70,71 (= auth 168,169,171,172,173)
# Using HMG-only to avoid VRAM pressure from full-length 509 aa sequence
# 50 designs, 7 cycles, ipTM threshold 0.7
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/boltz2/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export MPLBACKEND=Agg

SOX9_SEQ="PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

echo "=== Protein Hunter SOX9 binder — $(date) ==="

cd /home/ayushregmibagale/Protein-Hunter

python /home/ayushregmibagale/Protein-Hunter/boltz_ph/design.py \
    --name SOX9_binder \
    --protein_seqs "$SOX9_SEQ" \
    --contact_residues "66,67,69,70,71" \
    --min_protein_length 60 \
    --max_protein_length 100 \
    --num_designs 50 \
    --num_cycles 7 \
    --msa_mode single \
    --gpu_id 0 \
    --high_iptm_threshold 0.7 \
    --high_plddt_threshold 0.7 \
    --alanine_bias \
    --plot \
    --save_dir /mnt/d/gbm-protein-switch/protein_hunter/sox9_binder_hmg

echo "=== Protein Hunter SOX9 done — $(date) ==="

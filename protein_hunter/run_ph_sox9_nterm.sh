#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# Protein Hunter — SOX9 N-terminal face binder design.
# contact_residues 5,10,15,19,25 = label_seq_id for auth 107,112,117,121,127.
# These are the N-terminal half of the SOX9 HMG domain.
# Goal: binders whose N-terminus is on the SOX2-proximal face for CaspC+GGS+BinderSOX9 fusion.
# 50 designs, 7 cycles, same protocol as original sox9 campaign.
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/boltz2/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export MPLBACKEND=Agg

SOX9_SEQ="PHVKRPMNAFMVWAQAARRKLADQYPHLHNAELSKTLGKLWRLLNESEKRPFVEEAERLRVQHKKDHPDYKY"

echo "=== Protein Hunter SOX9 N-term binders — $(date) ==="
echo "Hotspots: label_seq 5,10,15,19,25 (auth 107,112,117,121,127) — N-terminal half of HMG"

cd /home/ayushregmibagale/Protein-Hunter

python /home/ayushregmibagale/Protein-Hunter/boltz_ph/design.py \
    --name SOX9_nterm \
    --protein_seqs "$SOX9_SEQ" \
    --contact_residues "5,10,15,19,25" \
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
    --save_dir /mnt/d/gbm-protein-switch/protein_hunter/sox9_nterm

echo "=== Protein Hunter SOX9 N-term done — $(date) ==="

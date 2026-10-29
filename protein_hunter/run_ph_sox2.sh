#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# Protein Hunter — SOX2 binder design
# Target: SOX2 (317 aa), hotspots A62,A70,A74,A81 (HMG wing-2 face)
# 50 designs, 7 cycles, ipTM threshold 0.7
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/boltz2/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export MPLBACKEND=Agg

SOX2_SEQ="MYNMMETELKPPGPQQTSGGGGGNSTAAAAGGNQKNSPDRVKRPMNAFMVWSRGQRRKMAQENPKMHNSEISKRLGAEWKLLSETEKRPFIDEAKRLRALHMKEHPDYKYRPRRKTKTLMKKDKYTLPGGLLAPGGNSMASGVGVGAGLGAGVNQRMDSYAHMNGWSNGSYSMMQDQLGYPQHPGLNAHGAAQMQPMHRYDVSALQYNSMTSSQTYMNGSPTYSMSYSQQGTPGMALGSMGSVVKSEASSSPPVVTSSSHSRAPCQAGDLRDMISMYLPGAEVPEPAAPSRLHMSQHYQSGPVPGTAINGTLPLSHM"

echo "=== Protein Hunter SOX2 binder — $(date) ==="

cd /home/ayushregmibagale/Protein-Hunter

python /home/ayushregmibagale/Protein-Hunter/boltz_ph/design.py \
    --name SOX2_binder \
    --protein_seqs "$SOX2_SEQ" \
    --contact_residues "62,70,74,81" \
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
    --save_dir /mnt/d/gbm-protein-switch/protein_hunter/sox2_binder

echo "=== Protein Hunter SOX2 done — $(date) ==="

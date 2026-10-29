#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# Protein Hunter — SOX2/SOX9 interface bridging binder (truncated fragments)
# Target: SOX2 aa175-215 (41aa) chain B + SOX9 aa450-509 (60aa) chain C
# Hotspots in fragment coords: SOX2 pos 19,25 | SOX9 pos 17,56,59
# Binder (chain A) must contact BOTH chains simultaneously
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/boltz2/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export MPLBACKEND=Agg

SOX2_FRAG="QDQLGYPQHPGLNAHGAAQMQPMHRYDVSALQYNSMTSSQT"
SOX9_FRAG="SYYSHAAGQGTGLYSTFTYMNPAQRPMYTPIADTSGVPSIPQTHSPQHWEQPVYTQLTRP"

echo "=== Protein Hunter bridge binder (truncated) — $(date) ==="

cd /home/ayushregmibagale/Protein-Hunter

python /home/ayushregmibagale/Protein-Hunter/boltz_ph/design.py \
    --name bridge_SOX2_SOX9 \
    --protein_seqs "${SOX2_FRAG}:${SOX9_FRAG}" \
    --contact_residues "19,25|17,56,59" \
    --min_protein_length 60 \
    --max_protein_length 100 \
    --num_designs 50 \
    --num_cycles 7 \
    --msa_mode single \
    --gpu_id 0 \
    --high_iptm_threshold 0.65 \
    --high_plddt_threshold 0.65 \
    --alanine_bias \
    --plot \
    --save_dir /mnt/d/gbm-protein-switch/protein_hunter/bridge_binder

echo "=== Protein Hunter bridge done — $(date) ==="

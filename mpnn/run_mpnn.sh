#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# ProteinMPNN sequence design on RFdiffusion binder backbones
# Run for both SOX2 and SOX9 backbones
# Generates 8 sequences per backbone → ~1600 sequences per target
set -e

export PATH="/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin:$PATH"
MPNN=/home/ayushregmibagale/projects/ProteinMPNN

for TARGET in sox2 sox9; do
    BACKBONES=/mnt/d/gbm-protein-switch/rfdiffusion/${TARGET}_binder
    OUTDIR=/mnt/d/gbm-protein-switch/mpnn/${TARGET}_binder
    mkdir -p $OUTDIR/parsed $OUTDIR/assigned $OUTDIR/seqs

    echo "=== ProteinMPNN $TARGET — $(date) ==="

    # Parse PDBs into MPNN format
    python $MPNN/helper_scripts/parse_multiple_chains.py \
        --input_path $BACKBONES \
        --output_path $OUTDIR/parsed/parsed.jsonl

    # Fix target chains (A for sox2, B for sox9), design binder chain only
    FIXED_CHAIN=$([ "$TARGET" = "sox2" ] && echo "A" || echo "B")
    python $MPNN/helper_scripts/assign_fixed_chains.py \
        --input_path $OUTDIR/parsed/parsed.jsonl \
        --output_path $OUTDIR/assigned/assigned.jsonl \
        --chain_list "$FIXED_CHAIN"

    # Design sequences (8 per backbone, temp=0.1 for stringency)
    python $MPNN/protein_mpnn_run.py \
        --jsonl_path $OUTDIR/parsed/parsed.jsonl \
        --fixed_positions_jsonl $OUTDIR/assigned/assigned.jsonl \
        --out_folder $OUTDIR/seqs \
        --num_seq_per_target 8 \
        --sampling_temp "0.1" \
        --score_only 0 \
        --use_soluble_model

    echo "=== ProteinMPNN $TARGET done — $(date) ==="
done

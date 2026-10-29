#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# BoltzGen SOX9 binder design
# Target: SOX9 HMG domain (4EUW chain A, 2.77 Å) + DNA context
# Hotspots: label_seq_id 46,50,57,64,69 (auth 148,152,159,166,171)
set -e

export PATH="$HOME/miniconda3/envs/boltzgen/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

SPEC=~/gbm-protein-switch/boltzgen/specs/sox9_binder.yaml
OUT=~/gbm-protein-switch/boltzgen/output/sox9

echo "=== BoltzGen SOX9 — $(date) ==="
echo "Checking spec first..."
boltzgen check "$SPEC" && echo "Spec OK"

echo "Running BoltzGen SOX9 (num_designs=1000, budget=20)..."
boltzgen run "$SPEC" \
  --output "$OUT" \
  --protocol protein-anything \
  --num_designs 1000 \
  --budget 20

echo "=== BoltzGen SOX9 done — $(date) ==="

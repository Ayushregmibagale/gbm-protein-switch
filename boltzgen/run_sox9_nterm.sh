#!/bin/bash
# BoltzGen SOX9 N-terminal face binder design.
# Hotspots on the N-terminal half of HMG: label_seq 5,10,15,19,25 (auth 107,112,117,121,127).
# Goal: binders whose N-terminus is on the SOX2-proximal face for CaspC+GGS+BinderSOX9 fusion.
set -e

export PATH="$HOME/miniconda3/envs/boltzgen/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

SPEC=~/gbm-protein-switch/boltzgen/specs/sox9_nterm_binder.yaml
OUT=~/gbm-protein-switch/boltzgen/output/sox9_nterm

echo "=== BoltzGen SOX9 N-term — $(date) ==="
boltzgen check "$SPEC" && echo "Spec OK"

boltzgen run "$SPEC" \
  --output "$OUT" \
  --protocol protein-anything \
  --num_designs 1000 \
  --budget 20

echo "=== BoltzGen SOX9 N-term done — $(date) ==="

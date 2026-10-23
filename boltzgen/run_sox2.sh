#!/bin/bash
# BoltzGen SOX2 binder design
# Target: SOX2 HMG domain (6WX8 chain B, 2.3 Å)
# Hotspots: label_seq_id 24,32,36,43 (auth 62,70,74,81)
set -e

export PATH="$HOME/miniconda3/envs/boltzgen/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

SPEC=~/gbm-protein-switch/boltzgen/specs/sox2_binder.yaml
OUT=~/gbm-protein-switch/boltzgen/output/sox2

echo "=== BoltzGen SOX2 — $(date) ==="
echo "Checking spec first..."
boltzgen check "$SPEC" && echo "Spec OK"

echo "Running BoltzGen SOX2 (num_designs=1000, budget=20)..."
boltzgen run "$SPEC" \
  --output "$OUT" \
  --protocol protein-anything \
  --num_designs 1000 \
  --budget 20

echo "=== BoltzGen SOX2 done — $(date) ==="

#!/bin/bash
# BoltzGen — interface bridging binder (SOX2+SOX9 co-complex)
# Binder must contact A193,A199 (SOX2) AND B466,B505,B508 (SOX9) simultaneously
set -e

export PATH="$HOME/miniconda3/envs/boltzgen/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

SPEC=~/gbm-protein-switch/boltzgen/specs/bridge_binder.yaml
OUT=~/gbm-protein-switch/boltzgen/output/bridge

echo "=== BoltzGen bridge binder — $(date) ==="
echo "Checking spec first..."
boltzgen check "$SPEC" && echo "Spec OK"

echo "Running BoltzGen bridge (num_designs=1000, budget=20)..."
boltzgen run "$SPEC" \
  --output "$OUT" \
  --protocol protein-anything \
  --num_designs 1000 \
  --budget 20

echo "=== BoltzGen bridge done — $(date) ==="

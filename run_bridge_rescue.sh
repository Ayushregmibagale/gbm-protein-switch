#!/bin/bash
# Master launcher for bridge binder rescue.
# GPU tasks are run sequentially; linker screen runs in parallel on CPU.
# Phase 1: RFdiffusion bridge both chains (~30-60 min)
# Phase 2: Boltz-2 linker screen (parallel, CPU)
# Phase 3: BindCraft no-hotspot (overnight)
# Phase 4: RFdiffusion SOX9 N-term binders (~30 min)
set -euo pipefail

LOG=/mnt/d/gbm-protein-switch/bridge_rescue.log
GBM=/mnt/d/gbm-protein-switch
PYTHON=/home/ayushregmibagale/miniconda3/envs/boltz2/bin/python3

exec > >(tee -a "$LOG") 2>&1
echo "=============================="
echo "Bridge rescue started: $(date)"
echo "=============================="

# Phase 2: linker screen in background (CPU, no GPU needed)
echo "[$(date)] Launching linker screen in background..."
$PYTHON $GBM/selection/bridge_linker_screen.py &
LINKER_PID=$!
echo "Linker screen PID: $LINKER_PID"

# Phase 1: RFdiffusion dual-chain bridge
echo "[$(date)] === Phase 1: RFdiffusion dual-chain bridge ==="
bash $GBM/rfdiffusion/run_rfdiffusion_bridge_both.sh
echo "[$(date)] Phase 1 DONE"

# Phase 3: BindCraft no-hotspot (runs overnight)
echo "[$(date)] === Phase 3: BindCraft no-hotspot ==="
bash $GBM/bindcraft/run_bridge_nohot.sh
echo "[$(date)] Phase 3 DONE"

# Phase 4: RFdiffusion SOX9 N-term binders
echo "[$(date)] === Phase 4: RFdiffusion SOX9 N-term ==="
bash $GBM/rfdiffusion/run_rfdiffusion_sox9_nterm.sh
echo "[$(date)] Phase 4 DONE"

# Wait for linker screen
echo "[$(date)] Waiting for linker screen (PID $LINKER_PID)..."
wait $LINKER_PID && echo "[$(date)] Linker screen DONE" || echo "[$(date)] Linker screen FAILED"

echo "=============================="
echo "Bridge rescue complete: $(date)"
echo "=============================="
echo "Results:"
echo "  RFdiffusion bridge:  $GBM/rfdiffusion/bridge_both_chains/ (200 backbones)"
echo "  Linker screen:       $GBM/selection/bridge_linker_outputs/linker_screen_results.tsv"
echo "  BindCraft nohot:     $GBM/bindcraft/designs/sox2_sox9_bridge_nohot/"
echo "  RFdiffusion sox9_nt: $GBM/rfdiffusion/sox9_nterm/ (200 backbones)"

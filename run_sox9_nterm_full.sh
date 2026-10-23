#!/bin/bash
# SOX9 N-terminal binder full campaign: RFdiffusion + Protein Hunter + BoltzGen + BindCraft.
# Waits for the bridge_rescue master (PID passed as $1, or auto-detected) to finish first.
# GPU tasks run sequentially; BoltzGen runs in background on GPU in parallel with PH.
set -euo pipefail

GBM=/mnt/d/gbm-protein-switch
LOG=$GBM/sox9_nterm_full.log
PYTHON_BOLTZ=/home/ayushregmibagale/miniconda3/envs/boltz2/bin/python3

exec > >(tee -a "$LOG") 2>&1

BRIDGE_PID=${1:-$(pgrep -f "run_bridge_rescue" | head -1 || true)}

if [[ -n "$BRIDGE_PID" ]] && kill -0 "$BRIDGE_PID" 2>/dev/null; then
    echo "[$(date)] Waiting for bridge_rescue (PID $BRIDGE_PID) to finish..."
    while kill -0 "$BRIDGE_PID" 2>/dev/null; do sleep 60; done
    echo "[$(date)] bridge_rescue finished."
else
    echo "[$(date)] No bridge_rescue running — starting N-term campaign immediately."
fi

echo "=============================="
echo "SOX9 N-term full campaign started: $(date)"
echo "=============================="

# Phase A: Protein Hunter (GPU — boltz2 env)
echo "[$(date)] === Phase A: Protein Hunter SOX9 N-term (50 runs × 7 cycles) ==="
bash $GBM/protein_hunter/run_ph_sox9_nterm.sh
echo "[$(date)] Phase A DONE"

# Phase B: BoltzGen (GPU — boltzgen env)
echo "[$(date)] === Phase B: BoltzGen SOX9 N-term (1000 designs) ==="
bash $GBM/boltzgen/run_sox9_nterm.sh
echo "[$(date)] Phase B DONE"

# Phase C: BindCraft (GPU — BindCraft env, overnight)
echo "[$(date)] === Phase C: BindCraft SOX9 N-term ==="
bash $GBM/bindcraft/run_sox9_nterm_bc.sh
echo "[$(date)] Phase C DONE"

echo "=============================="
echo "SOX9 N-term full campaign COMPLETE: $(date)"
echo "=============================="
echo "Results:"
echo "  PH:        $GBM/protein_hunter/sox9_nterm/summary_all_runs.csv"
echo "  BoltzGen:  $GBM/boltzgen/output/sox9_nterm/"
echo "  BindCraft: $GBM/bindcraft/designs/sox9_nterm/"
echo "  RFdiffusion: $GBM/rfdiffusion/sox9_nterm/ (from bridge_rescue Phase 4)"
echo ""
echo "Next: run terminus direction diagnosis on all designs"
echo "  Select those with d_toward > 0 (N-terminus facing SOX2)"
echo "  Then: Stage5 flipped caspase geometry → MD"

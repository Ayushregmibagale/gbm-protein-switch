#!/bin/bash
# ProteinMPNN + Boltz-2 pipeline for RFdiffusion bridge binder designs.
# Waits for BindCraft no-hotspot (PID passed as $1) to free the GPU first.
set -euo pipefail

LOG=/mnt/d/gbm-protein-switch/rfdiffusion/mpnn_boltz_bridge.log
exec > >(tee -a "$LOG") 2>&1

BINDCRAFT_PID=${1:-$(pgrep -f "bindcraft.py.*bridge_nohot" | head -1 || true)}

if [[ -n "$BINDCRAFT_PID" ]] && kill -0 "$BINDCRAFT_PID" 2>/dev/null; then
    echo "[$(date)] Waiting for BindCraft no-hotspot (PID $BINDCRAFT_PID)..."
    while kill -0 "$BINDCRAFT_PID" 2>/dev/null; do sleep 60; done
    echo "[$(date)] BindCraft finished."
else
    echo "[$(date)] BindCraft not running — starting immediately."
fi

echo "[$(date)] === Running ProteinMPNN + Boltz-2 on top 20 bridge backbones ==="
export PATH="/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin:$PATH"
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}

/home/ayushregmibagale/miniconda3/envs/rf_diffusion/bin/python \
    /mnt/d/gbm-protein-switch/rfdiffusion/run_mpnn_boltz_bridge.py

echo "[$(date)] === Pipeline complete ==="

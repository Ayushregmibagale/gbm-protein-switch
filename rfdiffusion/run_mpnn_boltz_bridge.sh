#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
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

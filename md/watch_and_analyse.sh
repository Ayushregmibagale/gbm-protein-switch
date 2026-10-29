#!/usr/bin/env bash
# Wait for production MD (PID stored in $1 or auto-detect) to finish, then run analysis.
set -euo pipefail

MD_DIR=/mnt/d/gbm-protein-switch/md/lead_construct
SCRIPT=/mnt/d/gbm-protein-switch/md/run_lead_md.py
GMX=/home/ayushregmibagale/.conda/envs/gromacs/bin/gmx
PYTHON=/home/ayushregmibagale/.conda/envs/gromacs/bin/python3

PROD_PID=${1:-$(pgrep -f "deffnm prod" 2>/dev/null | head -1 || true)}

if [[ -z "$PROD_PID" ]]; then
    echo "No running prod mdrun found — checking state..."
else
    echo "Watching prod mdrun PID $PROD_PID..."
    while kill -0 "$PROD_PID" 2>/dev/null; do
        sleep 60
    done
    echo "mdrun finished."
fi

echo "Setting state to PROD_DONE..."
echo "PROD_DONE" > "$MD_DIR/state.txt"

echo "Running analysis..."
cd /mnt/d/gbm-protein-switch
$PYTHON md/run_lead_md.py 2>&1 | tee -a "$MD_DIR/analysis.log"

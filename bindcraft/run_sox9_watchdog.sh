#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# Watchdog: restart SOX9 BindCraft v3 on crash until 20 accepted designs found
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export CUDA_HOME="$HOME/miniconda3"
export PATH="$HOME/miniconda3/bin:$PATH"
export MPLBACKEND=Agg

DESIGN_DIR=/mnt/d/gbm-protein-switch/bindcraft/designs/sox9_binder_xtal_v3
LOG=~/gbm-protein-switch/bindcraft/logs/sox9_bindcraft_v3.log
TARGET=20

count_accepted() { find "$DESIGN_DIR/Accepted" -name "*.pdb" 2>/dev/null | wc -l; }

while [ "$(count_accepted)" -lt "$TARGET" ]; do
    echo "=== SOX9 BindCraft v3 start/restart $(date) — $(count_accepted)/$TARGET accepted ===" >> "$LOG"
    cd /home/ayushregmibagale/projects/BindCraft
    /home/ayushregmibagale/miniconda3/envs/BindCraft/bin/python bindcraft.py \
        -s settings_target/SOX9_miniprotein.json \
        -f settings_filters/default_filters.json \
        -a settings_advanced/default_4stage_multimer_noanim.json >> "$LOG" 2>&1
    echo "=== SOX9 BindCraft v3 exited $(date) — $(count_accepted)/$TARGET accepted ===" >> "$LOG"
    [ "$(count_accepted)" -lt "$TARGET" ] && sleep 5
done
echo "=== SOX9 BindCraft v3 COMPLETE: $(count_accepted) accepted designs ===" >> "$LOG"

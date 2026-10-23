#!/bin/bash
# Watchdog: restart SOX2 BindCraft on crash until 20 accepted designs found
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export CUDA_HOME="$HOME/miniconda3"
export PATH="$HOME/miniconda3/bin:$PATH"
export MPLBACKEND=Agg

DESIGN_DIR=/mnt/d/gbm-protein-switch/bindcraft/designs/sox2_binder_xtal
LOG=~/gbm-protein-switch/bindcraft/sox2_xtal_bindcraft.log
TARGET=13

count_accepted() { find "$DESIGN_DIR/Accepted" -name "*.pdb" 2>/dev/null | wc -l; }

while [ "$(count_accepted)" -lt "$TARGET" ]; do
    echo "=== SOX2 BindCraft start/restart $(date) — $(count_accepted)/$TARGET accepted ===" >> "$LOG"
    cd /home/ayushregmibagale/projects/BindCraft
    /home/ayushregmibagale/miniconda3/envs/BindCraft/bin/python bindcraft.py \
        -s settings_target/SOX2_miniprotein.json \
        -f settings_filters/sox2_relaxed_filters.json \
        -a settings_advanced/default_4stage_multimer_noanim.json >> "$LOG" 2>&1
    echo "=== SOX2 BindCraft exited $(date) — $(count_accepted)/$TARGET accepted ===" >> "$LOG"
    [ "$(count_accepted)" -lt "$TARGET" ] && sleep 5
done
echo "=== SOX2 BindCraft COMPLETE: $(count_accepted) accepted designs ===" >> "$LOG"

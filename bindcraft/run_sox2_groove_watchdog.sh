#!/bin/bash
# Watchdog: SOX2 groove binding site (H1/L1/H2 junction, hotspots B59,B61,B62,B63,B70)
# Run after current SOX2 campaign completes or GPU is free
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export CUDA_HOME="$HOME/miniconda3"
export PATH="$HOME/miniconda3/bin:$PATH"
export MPLBACKEND=Agg

DESIGN_DIR=/mnt/d/gbm-protein-switch/bindcraft/designs/sox2_binder_groove
LOG=~/gbm-protein-switch/bindcraft/logs/sox2_groove_bindcraft.log
TARGET=20

count_accepted() { find "$DESIGN_DIR/Accepted" -name "*.pdb" 2>/dev/null | wc -l; }

while [ "$(count_accepted)" -lt "$TARGET" ]; do
    echo "=== SOX2 groove start/restart $(date) — $(count_accepted)/$TARGET accepted ===" >> "$LOG"
    cd /home/ayushregmibagale/projects/BindCraft
    /home/ayushregmibagale/miniconda3/envs/BindCraft/bin/python bindcraft.py \
        -s settings_target/SOX2_groove_miniprotein.json \
        -f settings_filters/sox2_relaxed_filters.json \
        -a settings_advanced/default_4stage_multimer_noanim.json >> "$LOG" 2>&1
    echo "=== SOX2 groove exited $(date) — $(count_accepted)/$TARGET accepted ===" >> "$LOG"
    [ "$(count_accepted)" -lt "$TARGET" ] && sleep 5
done
echo "=== SOX2 groove COMPLETE: $(count_accepted) accepted designs ===" >> "$LOG"

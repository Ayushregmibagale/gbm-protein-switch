#!/bin/bash
# Historical campaign launcher: explicit opt-in is required.
if [[ "${GBM_ALLOW_LEGACY_CAMPAIGNS:-0}" != "1" ]]; then
    echo "Archived campaign launcher. See WORKFLOW.md for supported entry points." >&2
    echo "For deliberate legacy execution set GBM_ALLOW_LEGACY_CAMPAIGNS=1." >&2
    exit 2
fi
# Poll until both crystal-structure BindCraft runs finish, then run Protein Hunter
# Auto-cancels if GPU temp exceeds 70°C
SOX2_PID=4146777
SOX9_PID=1711852
SOX2_DESIGN_DIR=/mnt/d/gbm-protein-switch/bindcraft/designs/sox2_binder_xtal
SOX9_DESIGN_DIR=/mnt/d/gbm-protein-switch/bindcraft/designs/sox9_binder_xtal_v3
GPU_TEMP_LIMIT=70

gpu_temp() {
    nvidia-smi --query-gpu=temperature.gpu --format=csv,noheader,nounits 2>/dev/null | head -1
}

count_accepted() {
    local dir=$1
    find "$dir/Accepted" -name "*.pdb" 2>/dev/null | wc -l
}

echo "=== Pipeline launcher started $(date) ==="
echo "SOX2 PID=$SOX2_PID (crystal 6WX8)  SOX9 PID=$SOX9_PID (crystal 4EUW)  GPU limit=${GPU_TEMP_LIMIT}°C"

while kill -0 $SOX2_PID 2>/dev/null || kill -0 $SOX9_PID 2>/dev/null; do
    TEMP=$(gpu_temp)
    SOX2_STATUS=$(kill -0 $SOX2_PID 2>/dev/null && echo "running($(count_accepted $SOX2_DESIGN_DIR)acc)" || echo "done($(count_accepted $SOX2_DESIGN_DIR)acc)")
    SOX9_STATUS=$(kill -0 $SOX9_PID 2>/dev/null && echo "running($(count_accepted $SOX9_DESIGN_DIR)acc)" || echo "done($(count_accepted $SOX9_DESIGN_DIR)acc)")
    echo "$(date '+%H:%M') SOX2=$SOX2_STATUS  SOX9=$SOX9_STATUS  GPU=${TEMP}°C"

    if [ -n "$TEMP" ] && [ "$TEMP" -ge "$GPU_TEMP_LIMIT" ]; then
        echo "!!! GPU temp ${TEMP}°C >= ${GPU_TEMP_LIMIT}°C — killing all jobs $(date) !!!"
        kill $SOX2_PID $SOX9_PID 2>/dev/null
        echo "Jobs killed. Pipeline aborted."
        exit 1
    fi

    sleep 300
done

echo "=== Both BindCraft done — $(date) ==="

# Safety check: don't proceed with empty design folders
SOX2_ACC=$(count_accepted $SOX2_DESIGN_DIR)
SOX9_ACC=$(count_accepted $SOX9_DESIGN_DIR)
echo "Accepted designs: SOX2=$SOX2_ACC  SOX9=$SOX9_ACC"

if [ "$SOX2_ACC" -eq 0 ] && [ "$SOX9_ACC" -eq 0 ]; then
    echo "ERROR: No accepted designs from either BindCraft run — aborting Protein Hunter."
    exit 1
fi

if [ "$SOX2_ACC" -gt 0 ]; then
    echo "--- Protein Hunter SOX2 starting $(date) ---"
    bash /mnt/d/gbm-protein-switch/protein_hunter/run_ph_sox2.sh
else
    echo "SKIP: SOX2 BindCraft produced 0 accepted designs — skipping Protein Hunter SOX2"
fi

if [ "$SOX9_ACC" -gt 0 ]; then
    echo "--- Protein Hunter SOX9 starting $(date) ---"
    bash /mnt/d/gbm-protein-switch/protein_hunter/run_ph_sox9.sh
else
    echo "SKIP: SOX9 BindCraft produced 0 accepted designs — skipping Protein Hunter SOX9"
fi

echo "=== Protein Hunter complete — starting BoltzGen $(date) ==="

echo "--- BoltzGen check SOX2 ---"
export PATH="$HOME/miniconda3/envs/boltzgen/bin:$PATH"
boltzgen check /mnt/d/gbm-protein-switch/boltzgen/specs/sox2_binder.yaml \
    && echo "SOX2 spec OK" || { echo "SOX2 spec FAILED — skipping BoltzGen SOX2"; }

echo "--- BoltzGen SOX2 starting $(date) ---"
bash /mnt/d/gbm-protein-switch/boltzgen/run_sox2.sh \
    && echo "BoltzGen SOX2 done" || echo "BoltzGen SOX2 FAILED"

echo "--- BoltzGen SOX9 starting $(date) ---"
bash /mnt/d/gbm-protein-switch/boltzgen/run_sox9.sh \
    && echo "BoltzGen SOX9 done" || echo "BoltzGen SOX9 FAILED"

echo "=== Full pipeline complete $(date) ==="

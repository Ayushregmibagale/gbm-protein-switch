#!/bin/bash
# Chain: SOX9 BindCraft → RFdiffusion SOX2 → RFdiffusion SOX9 → ProteinMPNN
# Run this after SOX2 BindCraft finishes (or run separately)
set -e

LOG=/mnt/d/gbm-protein-switch/bindcraft/pipeline.log
exec >> $LOG 2>&1

echo "=== Pipeline started $(date) ==="

# Step 1: SOX9 BindCraft
echo "--- SOX9 BindCraft starting $(date) ---"
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export CUDA_HOME="$HOME/miniconda3"
export PATH="$HOME/miniconda3/bin:$PATH"
export MPLBACKEND=Agg

cd /home/ayushregmibagale/projects/BindCraft
/home/ayushregmibagale/miniconda3/envs/BindCraft/bin/python bindcraft.py \
  -s settings_target/SOX9_miniprotein.json \
  -f settings_filters/default_filters.json \
  -a settings_advanced/default_4stage_multimer_noanim.json
echo "--- SOX9 BindCraft done $(date) ---"

# Step 2: RFdiffusion SOX2
echo "--- RFdiffusion SOX2 starting $(date) ---"
bash /mnt/d/gbm-protein-switch/rfdiffusion/run_rfdiffusion_sox2.sh
echo "--- RFdiffusion SOX2 done $(date) ---"

# Step 3: RFdiffusion SOX9
echo "--- RFdiffusion SOX9 starting $(date) ---"
bash /mnt/d/gbm-protein-switch/rfdiffusion/run_rfdiffusion_sox9.sh
echo "--- RFdiffusion SOX9 done $(date) ---"

# Step 4: ProteinMPNN on all backbones
echo "--- ProteinMPNN starting $(date) ---"
bash /mnt/d/gbm-protein-switch/mpnn/run_mpnn.sh
echo "--- ProteinMPNN done $(date) ---"

echo "=== Full pipeline complete $(date) ==="

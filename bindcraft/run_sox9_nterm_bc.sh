#!/bin/bash
# BindCraft — SOX9 N-terminal face binder design.
# Hotspots: A107,A112,A117,A121,A127 (N-terminal half of 4EUW HMG domain).
# Goal: binders whose N-terminus is on the SOX2-proximal face for CaspC+GGS+BinderSOX9 fusion.
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export CUDA_HOME="$HOME/miniconda3"
export PATH="$HOME/miniconda3/bin:$PATH"
export MPLBACKEND=Agg

cd /home/ayushregmibagale/projects/BindCraft

echo "=== BindCraft SOX9 N-term binders — $(date) ==="

/home/ayushregmibagale/miniconda3/envs/BindCraft/bin/python bindcraft.py \
  -s settings_target/SOX9_nterm_miniprotein.json \
  -f settings_filters/default_filters.json \
  -a settings_advanced/default_4stage_multimer_noanim.json

echo "=== BindCraft SOX9 N-term done — $(date) ==="

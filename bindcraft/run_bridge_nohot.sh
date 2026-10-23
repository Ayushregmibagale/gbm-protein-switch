#!/bin/bash
# BindCraft bridge binder — no hotspot constraint, very relaxed filters.
# Rationale: with pinned hotspots (A193,A199,B466,B505,B508) the AF2 optimizer
# achieves i_pTM ~0.12-0.16 and never passes the base AF2 filter (233 trajectories, 0 accepted).
# Removing hotspots lets BindCraft find natural bridging contacts freely.
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export CUDA_HOME="$HOME/miniconda3"
export PATH="$HOME/miniconda3/bin:$PATH"
export MPLBACKEND=Agg

cd /home/ayushregmibagale/projects/BindCraft

echo "=== BindCraft bridge nohot — $(date) ==="

/home/ayushregmibagale/miniconda3/envs/BindCraft/bin/python bindcraft.py \
  -s settings_target/SOX2_SOX9_bridge_nohot.json \
  -f settings_filters/very_relaxed_bridge_filters.json \
  -a settings_advanced/default_4stage_multimer_noanim.json

echo "=== BindCraft bridge nohot done — $(date) ==="

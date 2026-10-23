#!/bin/bash
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:${LD_LIBRARY_PATH}
export CUDA_HOME="$HOME/miniconda3"
export PATH="$HOME/miniconda3/bin:$PATH"
export MPLBACKEND=Agg

cd /home/ayushregmibagale/projects/BindCraft

/home/ayushregmibagale/miniconda3/envs/BindCraft/bin/python bindcraft.py \
  -s settings_target/SOX2_SOX9_bridge_fragment.json \
  -f settings_filters/default_filters.json \
  -a settings_advanced/default_4stage_multimer_noanim.json

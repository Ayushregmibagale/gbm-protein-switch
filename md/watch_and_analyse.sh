#!/usr/bin/env bash
# Compatibility entry point: validate explicit completed trajectory inputs.
# Process exit and state-file contents do not establish complete production data.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec "${PYTHON_BIN:-python3}" "$SCRIPT_DIR/../analysis/audit_control.py" "$@"

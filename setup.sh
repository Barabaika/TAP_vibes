#!/usr/bin/env bash
# Create an isolated environment and fetch checksum-pinned model assets.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
prefix=".tap-env"
assets="assets"
device="cpu"
metrics_only=false
while (($#)); do
  case "$1" in
    --prefix) prefix="${2:?Missing environment directory}"; shift 2 ;;
    --assets) assets="${2:?Missing assets directory}"; shift 2 ;;
    --device) device="${2:?Choose cpu or cuda}"; shift 2 ;;
    --metrics-only) metrics_only=true; shift ;;
    --help|-h)
      echo 'Usage: bash setup.sh [--device cpu|cuda] [--metrics-only] [--prefix .tap-env] [--assets assets]'
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || {
  echo 'Linux x86_64 is required by the pinned PSA executable. Use WSL2 on Windows.' >&2; exit 2;
}
[[ "$device" == cpu || "$device" == cuda ]] || { echo 'Device must be cpu or cuda' >&2; exit 2; }
if command -v micromamba >/dev/null 2>&1; then
  manager=micromamba
elif command -v conda >/dev/null 2>&1; then
  manager=conda
else
  echo 'Install Miniforge (conda) or micromamba, then rerun setup.sh. No sudo is needed.' >&2
  exit 2
fi
action=create
[[ ! -d "$prefix/conda-meta" ]] || action=install
"$manager" "$action" --yes --prefix "$prefix" --override-channels \
  --channel conda-forge --channel bioconda python=3.10 pip hmmer=3.4 libstdcxx-ng
run_env() { "$manager" run --prefix "$prefix" "$@"; }
run_env python -m pip install --upgrade pip
if "$metrics_only"; then
  run_env python -m pip install -e '.[test]'
  run_env tap-assets --psa-only --output "$assets"
else
  index=cpu
  [[ "$device" != cuda ]] || index=cu121
  run_env python -m pip install "torch==2.1.0+$index" --index-url "https://download.pytorch.org/whl/$index"
  run_env python -m pip install -e '.[fold,test]'
  run_env tap-assets --output "$assets"
fi
run_env python -m pip check
activation_prefix="$prefix"
[[ "$prefix" == /* ]] || activation_prefix="./$prefix"
echo "Setup complete. Activate with: $manager activate $activation_prefix"
echo "Or run: $manager run --prefix $prefix python -m pytest -q --run-structure --psa $assets/psa"

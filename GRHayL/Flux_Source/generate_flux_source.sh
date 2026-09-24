#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 OUTPUT_DIRECTORY" >&2
  exit 2
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Pin the NRPy 2 revision used to generate the checked-in C kernels.
nrpy_commit=46654c2e22682b2ffe2757022868c2293cd45640
checkout="$(mktemp -d)"
trap 'rm -rf -- "$checkout"' EXIT

# Fetch the exact revision, including when the upstream branch moves.
git clone --quiet --filter=blob:none --no-checkout --depth=1 \
  https://github.com/nrpy/nrpy.git "$checkout/nrpy"
git -C "$checkout/nrpy" fetch --quiet --depth=1 origin "$nrpy_commit"
git -C "$checkout/nrpy" checkout --quiet --detach "$nrpy_commit"
NRPY_ROOT="$checkout/nrpy" PYTHONPATH="$checkout/nrpy${PYTHONPATH:+:$PYTHONPATH}" \
  python3 "$script_dir/generate_flux_source.py" "$1"

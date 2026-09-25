#!/usr/bin/env bash
# Run the scoped Radiation M1 inventory.
#
# Stored THC_M1 values are read only from the repository-local fixture package.
# This runner never invokes a fixture generator, THC_M1, or Verification.
# Keep the runner compatible with the macOS system Bash 3.2: use indexed arrays.
set -euo pipefail

m1_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
m1_fixture_dir="$m1_root/Unit_Tests/data/m1_thcm1"
m1_jthick_fixture="$m1_root/Unit_Tests/data/jthick_thcm1.fixture"
m1_manifest="$m1_root/Unit_Tests/make.code.defn"
m1_fixture_audit="$m1_fixture_dir/audit_package.py"
m1_fixture_archive="$m1_fixture_dir/payloads.tar.gz"
cd "$m1_root"

if [[ ! -s "$m1_manifest" ]]; then
  echo "Missing scoped M1 inventory: $m1_manifest" >&2
  exit 1
fi

# The diffusion executable loads this paired Jthick fixture from a fixed
# repository-relative path.  Check it before compiling so a clean checkout
# cannot report a build failure that is actually a missing validation input.
if [[ ! -s "$m1_jthick_fixture" ]]; then
  echo "Missing repository-local M1 fixture: $m1_jthick_fixture" >&2
  exit 1
fi

m1_sources=()
while IFS= read -r m1_source; do
  m1_sources+=("$m1_source")
done < <(
  awk '
    /^SRCS[[:space:]]*=/ {
      active = 1
      sub(/^[^=]*=[[:space:]]*/, "")
    }
    !active { next }
    {
      continued = ($0 ~ /\\[[:space:]]*$/)
      sub(/[[:space:]]*\\[[:space:]]*$/, "")
      for(i = 1; i <= NF; ++i)
        if($i ~ /^unit_test_.*\.c$/) print $i
      if(!continued) exit
    }
  ' "$m1_manifest"
)

if (( ${#m1_sources[@]} == 0 )); then
  echo "Scoped M1 inventory has no SRCS entries: $m1_manifest" >&2
  exit 1
fi

# The manifest is the build/run input, so it cannot also be the oracle: a
# deleted SRCS entry would silently shrink coverage while every existence and
# duplicate check below still passed.  This independent inventory is compared
# against it in both directions, so an omission and an unowned addition are
# each a hard failure.
m1_expected=(
  unit_test_m1_closure_fallback.c
  unit_test_m1_diffusion_flux.c
  unit_test_m1_error_handling.c
  unit_test_m1_fd_jacobian.c
  unit_test_m1_neutrino_rusanov_flux.c
  unit_test_m1_neutrino_seeded_invariants.c
  unit_test_m1_neutrino_source_update.c
  unit_test_m1_rate_provider.c
  unit_test_m1_thcm1_blended_rusanov.c
  unit_test_rusanov_flux.c
)

m1_inventory_mismatch=0
for m1_want in "${m1_expected[@]}"; do
  m1_found=0
  for m1_have in "${m1_sources[@]}"; do
    if [[ "$m1_have" == "$m1_want" ]]; then m1_found=1; break; fi
  done
  if (( m1_found == 0 )); then
    echo "Scoped M1 inventory is missing an expected test: $m1_want" >&2
    m1_inventory_mismatch=1
  fi
done
for m1_have in "${m1_sources[@]}"; do
  m1_found=0
  for m1_want in "${m1_expected[@]}"; do
    if [[ "$m1_have" == "$m1_want" ]]; then m1_found=1; break; fi
  done
  if (( m1_found == 0 )); then
    echo "Scoped M1 inventory lists an unexpected test: $m1_have" >&2
    echo "Add it to m1_expected in $0 so the runner owns it." >&2
    m1_inventory_mismatch=1
  fi
done
if (( m1_inventory_mismatch != 0 )); then
  echo "Scoped M1 inventory does not match the runner's expected set." >&2
  exit 1
fi

if [[ ! -s "$m1_fixture_audit" ]]; then
  echo "Missing M1 fixture-package audit: $m1_fixture_audit" >&2
  exit 1
fi

# Clean checkouts carry the retained payloads in a compressed archive. Expand
# only when no plaintext payload is present; the audit must catch incomplete or
# modified local packages rather than silently replacing their contents.
m1_payload_present=0
for m1_payload in "$m1_fixture_dir"/*.dat "$m1_fixture_dir"/*.m1; do
  if [[ -e "$m1_payload" || -L "$m1_payload" ]]; then
    m1_payload_present=1
    break
  fi
done
if (( m1_payload_present == 0 )); then
  if [[ ! -s "$m1_fixture_archive" ]]; then
    echo "Missing repository-local M1 fixture archive: $m1_fixture_archive" >&2
    exit 1
  fi
  tar -xzf "$m1_fixture_archive" -C "$m1_fixture_dir"
fi
python3 "$m1_fixture_audit" --root "$m1_fixture_dir"

m1_tests=()
for m1_source in "${m1_sources[@]}"; do
  m1_test="${m1_source%.c}"
  for ((m1_seen_index = 0; m1_seen_index < ${#m1_tests[@]}; ++m1_seen_index)); do
    if [[ "${m1_tests[$m1_seen_index]}" == "$m1_test" ]]; then
      echo "Duplicate scoped M1 test: $m1_test" >&2
      exit 1
    fi
  done
  if [[ ! -f "$m1_root/Unit_Tests/$m1_source" ]]; then
    echo "Scoped M1 source is missing: Unit_Tests/$m1_source" >&2
    exit 1
  fi
  m1_tests+=("$m1_test")
done

# These are the checked-in inputs for the four stored-reference consumers.
# Their expected outputs are retained THC_M1 producer values from the checked-in
# fixture package; no local generator may replace a missing file during a
# normal run.
m1_fixture_files=(
  pointwise_closure_moments.m1
  stress_energy.m1
  m1_thcm1_instantaneous_sources.m1
  transport_four_point_d0.dat
  transport_four_point_d1.dat
  transport_four_point_d2.dat
  transport_four_point_varying_d0.dat
  transport_four_point_varying_d1.dat
  transport_four_point_varying_d2.dat
  transport_four_point_varying_controls.dat
  rusanov_neutrino.dat
  rusanov_neutrino_current.dat
  rusanov_generic.dat
)

case "${1:-}" in
  --build|--build-only)
    m1_targets=()
    for m1_test in "${m1_tests[@]}"; do
      m1_targets+=("test/$m1_test")
    done
    make "${m1_targets[@]}"
    if [[ "$1" == --build-only ]]; then exit 0; fi
    shift
    ;;
esac
if (( $# != 0 )); then
  echo "Usage: bash Unit_Tests/run_m1_tests.sh [--build|--build-only]" >&2
  exit 1
fi

for m1_fixture in "${m1_fixture_files[@]}"; do
  if [[ ! -s "$m1_fixture_dir/$m1_fixture" ]]; then
    echo "Missing repository-local M1 fixture: $m1_fixture_dir/$m1_fixture" >&2
    exit 1
  fi
done

# This runner uses the configured default build directory, not an installed libghl.
export LD_LIBRARY_PATH="$m1_root/build/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="$m1_root/build/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
# Distinct LLVM profiles prevent later executables from overwriting earlier ones.
export LLVM_PROFILE_FILE="${LLVM_PROFILE_FILE:-$m1_root/test/m1-%m-%p.profraw}"
for m1_test in "${m1_tests[@]}"; do
  echo "Running $m1_test"
  case "$m1_test" in
    unit_test_m1_neutrino_seeded_invariants|unit_test_m1_neutrino_rusanov_flux|unit_test_m1_thcm1_blended_rusanov|unit_test_rusanov_flux)
      "$m1_root/test/$m1_test" --fixture-dir "$m1_fixture_dir"
      ;;
    unit_test_m1_rate_provider)
      # This creates a temporary EOS input for provider coverage.  It is not a
      # trusted M1 output generator and is unrelated to the retired explicit-RHS
      # datagen path.
      "$m1_root/test/$m1_test" --generated-fixture
      ;;
    *) "$m1_root/test/$m1_test" ;;
  esac
done

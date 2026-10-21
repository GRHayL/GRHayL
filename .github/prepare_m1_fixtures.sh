#!/bin/bash
# Prepare a raw M1 fixture directory from raw or gzip-compressed payloads.
#
# Standalone usage:
#   .github/prepare_m1_fixtures.sh SOURCE_DIR DEST_DIR
#
# Sourcing usage (no caller set-options or traps are changed):
#   source .github/prepare_m1_fixtures.sh
#   prepare_m1_fixtures "$SOURCE_DIR" "$DEST_DIR"
#
# SOURCE_DIR supplies each expected member as either "name.bin" (raw) or
# "name.bin.gz" (gzip-compressed). Raw members take precedence when both
# forms are present. DEST_DIR must exist, be a directory, and be empty; it
# receives raw "name.bin" copies or decompressed payloads upon success.
# SOURCE_DIR is never modified. A missing member, a failed gzip
# decompression/integrity check, or a failed copy returns nonzero and names
# the member. Prepared members in DEST_DIR are removed on failure.

# Expected M1 fixture member inventory, shared by the standalone invocation
# and the .github/run_tests.sh M1 replay path.
m1_fixture_members=(
  pointwise_closure_moments.bin
  m1_thcm1_instantaneous_sources.bin
  stress_energy.bin
  rusanov_neutrino.bin
  rusanov_neutrino_current.bin
  rusanov_generic.bin
  transport_four_point_d0.bin
  transport_four_point_d1.bin
  transport_four_point_d2.bin
  transport_four_point_varying_d0.bin
  transport_four_point_varying_d1.bin
  transport_four_point_varying_d2.bin
  transport_four_point_varying_controls.bin
  jthick_thcm1.bin
)

prepare_m1_fixtures() {
  local source_dir="$1"
  local dest_dir="$2"
  # Failure paths are explicit checks only, so this function changes no
  # set-options and sourcing callers keep their own options/traps. The
  # ${arr[@]+"${arr[@]}"} expansions keep set -u callers safe while
  # prepared_members is still empty.
  local prepared_members=()

  if [[ ! -d "$source_dir" ]]; then
    echo "M1 fixture source is not a directory: $source_dir" >&2
    return 1
  fi
  if [[ ! -d "$dest_dir" ]]; then
    echo "M1 fixture destination is not a directory: $dest_dir" >&2
    return 1
  fi
  if [[ -n $(find "$dest_dir" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
    echo "M1 fixture destination is not empty: $dest_dir" >&2
    return 1
  fi

  local member source_member dest_member
  for member in "${m1_fixture_members[@]}"; do
    dest_member="$dest_dir/$member"
    if [[ -f "$source_dir/$member" ]]; then
      # Raw members take precedence over a supplied .gz counterpart.
      source_member="$source_dir/$member"
      if ! cp -- "$source_member" "$dest_member"; then
        echo "M1 fixture member failed to copy: $source_member" >&2
        rm -f -- "$dest_member" ${prepared_members[@]+"${prepared_members[@]}"} 2>/dev/null
        return 1
      fi
    elif [[ -f "$source_dir/$member.gz" ]]; then
      source_member="$source_dir/$member.gz"
      # Successful gzip exit includes the decompression integrity check.
      if ! gzip -dc -- "$source_member" >"$dest_member"; then
        echo "M1 fixture member failed gzip decompression: $source_member" >&2
        rm -f -- "$dest_member" ${prepared_members[@]+"${prepared_members[@]}"} 2>/dev/null
        return 1
      fi
    else
      echo "M1 fixture member is missing: $source_dir/$member (or $member.gz)" >&2
      rm -f -- ${prepared_members[@]+"${prepared_members[@]}"} 2>/dev/null
      return 1
    fi
    prepared_members+=("$dest_member")
  done

  return 0
}

# Standalone invocation only; sourcing executes nothing and changes no
# caller set-options or traps.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  set -e
  if [[ $# -ne 2 ]]; then
    echo "Usage: $0 SOURCE_DIR DEST_DIR" >&2
    exit 2
  fi
  prepare_m1_fixtures "$1" "$2"
fi

#!/bin/sh
# The gate for formal/, run locally and by CI: build the library with every warning an error, run
# the environment linters, refuse an unfinished proof or a new assumption in the sources, and audit
# through the environment what every declaration rests on.
#
#     ./formal/check.sh        # about 15 s once Mathlib is in place (`lake exe cache get`)
set -eu
cd "$(dirname "$0")"

lake build --wfail Array16
lake lint

if grep -nwE 'sorry|axiom|admit|native_decide' Array16.lean Array16/*.lean Array16/*/*.lean; then
  echo "an unfinished proof or a new assumption: the lines above" >&2
  exit 1
fi

audit=$(lake env lean Audit.lean 2>&1) || { printf '%s\n' "$audit" >&2; exit 1; }
printf '%s\n' "$audit"
printf '%s\n' "$audit" | grep -q 'declarations audited' || { echo "the audit reported no count" >&2; exit 1; }

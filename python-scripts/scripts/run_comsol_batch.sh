#!/usr/bin/env bash
# run_comsol_batch.sh -- solve a study of a COMSOL template with `comsol batch`
# (no GUI, no MPh, no LiveLink) and extract the eigenfrequencies.
#
# A licence is checked out only while the COMSOL processes run and is returned
# when they exit. Nothing is held before, between runs, or afterwards.
#
# Usage:
#   scripts/run_comsol_batch.sh                       # rib template, opt TE at X
#   scripts/run_comsol_batch.sh --study std2 --study std4       # opt TE + opt TM
#   scripts/run_comsol_batch.sh --study "opt TM" -p a=600[nm] -p t=400[nm]
#   scripts/run_comsol_batch.sh --outdir results/my_run --watch-license
#   scripts/run_comsol_batch.sh --check-license       # what the template needs
#   scripts/run_comsol_batch.sh --help                # every option
#
# Defaults: template comsol/omc_unitcell_rib.mph, study std2 ('opt TE'), all
# parameters as saved in the template (kF = pi/a, the X point), output in a
# new results/batch_<time stamp>/ folder:
#   solved_<study>.mph   the solved model (open it to look at the modes)
#   solve_<study>.log    COMSOL batch log of the solve
#   eig_<study>.csv      one row per mode: f_real_Hz, f_imag_Hz, integrals
#   summary.json         the same plus light line, Q, energy fractions
#   license_watch.json   only with --watch-license
#   OmcEigExtract.java/.class, extract.log, extract_stdout.txt   extraction step
#
# The rib studies (comsol/build_omc_unitcell_rib.m):
#   std1 'mech sym'   std3 'mech antisym'   std2 'opt TE'   std4 'opt TM'
#
# Parameter overrides (-p NAME=VALUE) are passed with COMSOL's -pname/-plist.
# A value is a COMSOL expression without commas or spaces, e.g. 600[nm].
#
# LICENCE FEATURES A RUN CHECKS OUT
# Measured on 2026-11-23 with --watch-license (COMSOL 6.3 build 420, server
# 27005@license.seas.harvard.edu, rib template, opt TE and opt TM at X):
#
#   solve (about 3 min; default --license batch = COMSOL's -usebatchlic):
#       COMSOLBATCH  WAVEOPTICSBATCH  MEMSBATCH  CADIMPORTBATCH  CLUSTERNODE
#       (+ COMSOLUSER, CADIMPORTUSER)
#     No interactive COMSOL, WAVEOPTICS or RF seat is used, so the solve runs
#     even when WAVEOPTICS is 3/3. MEMSBATCH (2 seats) and CADIMPORTBATCH are
#     taken although only an optical study runs: the template lists COMSOL,
#     WAVEOPTICS, MEMS and CADIMPORT as its licences (see --check-license).
#
#   extract (about 45 s, seats held for the last ~30 s of it):
#       COMSOL  WAVEOPTICS  MEMS  CADIMPORT   (+ COMSOLUSER, CADIMPORTUSER)
#     These are the ORDINARY features. COMSOL refuses -usebatchlic for a class
#     file ("Input file does not support batch licenses"), so this step needs
#     an interactive WAVEOPTICS and MEMS seat for a few seconds. If none is
#     free the solved .mph is kept; get the numbers later with
#       scripts/wait_for_comsol_seat.sh -f WAVEOPTICS -f MEMS --all -r 3 -- \
#           scripts/run_comsol_batch.sh --extract-only --outdir <folder> --study <tag>
#
#   Everything was released within 3 s of each process exiting.
#
# Batch mode by itself does NOT select the batch pool: -usebatchlic does.
# With --license interactive the flag is left out and, per the COMSOL
# documentation, the ordinary features are used for the solve as well (not
# measured here). Whether the solve falls back to RFBATCH when WAVEOPTICSBATCH
# is full was not measured either.
#
# The work is done by scripts/comsol_batch.py (standard library only); this
# wrapper only picks a Python. The COMSOL options used are the documented
# ones of COMSOL 6.3 (`comsolbatch -help`): -inputfile -outputfile -study
# -pname -plist -batchlog -nosave -usebatchlic -checklicense.
#
# Environment (every path is configurable; nothing is hard-coded):
#   PYTHON          Python 3.8+ to use (default: python3, then python)
#   COMSOL_ROOT     the .../COMSOL63/Multiphysics folder (or --comsol-root).
#                   Default: the newest COMSOL under the usual locations
#                   (/mnt/c/Program Files/COMSOL, C:/Program Files/COMSOL,
#                   /usr/local, /opt, /Applications).
#   COMSOL_BATCH    the batch executable, if it is not <root>/bin/win64/
#                   comsolbatch.exe or <root>/bin/comsol
#   COMSOL_JAVAC    javac, if not the one under <root>/java
#   COMSOL_BATCH_LICENSE   batch (default) | interactive, as --license
#   COMSOL_LICENSE_SERVER, LMUTIL, COMSOL_LM_USER   for --watch-license
#
# Exit status: 0 ok, 1 COMSOL or extraction failed, 2 template not found,
#              3 failed with a licence error (no seat), 64 usage, 130 Ctrl-C.
# On a failure the end of the COMSOL log is printed to stderr, so
#   scripts/wait_for_comsol_seat.sh --batch -r 3 -- scripts/run_comsol_batch.sh ...
# recognises a lost race for the seat and retries.

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "${1:-}" in
    -h|--help)
        awk 'NR == 1 { next } !/^#/ { exit } { sub(/^# ?/, ""); print }' "$0"
        echo
        echo "---- options (scripts/comsol_batch.py) ----"
        ;;
esac

PY="${PYTHON:-}"
if [ -z "$PY" ]; then
    for c in python3 python; do
        if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
    done
fi
if [ -z "$PY" ]; then
    echo "run_comsol_batch: no python3/python on PATH; set PYTHON=/path/to/python" >&2
    exit 1
fi

# A Windows python started from WSL does not understand /mnt/c/... script paths.
SCRIPT="$HERE/comsol_batch.py"
case "$PY" in
    *.exe) command -v wslpath >/dev/null 2>&1 && SCRIPT="$(wslpath -w "$SCRIPT")" ;;
esac

exec "$PY" "$SCRIPT" "$@"

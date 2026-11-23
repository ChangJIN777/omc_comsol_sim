#!/usr/bin/env bash
# check_comsol_license.sh -- show free COMSOL licence seats on the Harvard
# FlexLM server, WITHOUT checking a licence out.
#
# check_comsol.py checks the licence by starting COMSOL via MPh, which holds a
# seat while it runs. This script only asks the server (lmutil lmstat), so it
# is safe to run at any time, e.g. before an MPh / LiveLink job.
#
# Usage:
#   scripts/check_comsol_license.sh                 # features these models use
#   scripts/check_comsol_license.sh -a              # every feature on the server
#   scripts/check_comsol_license.sh -f WAVEOPTICS   # one feature (repeatable)
#   scripts/check_comsol_license.sh -u              # also list who holds seats
#   scripts/check_comsol_license.sh -s 27005@host   # another licence server
#
# Environment overrides:
#   COMSOL_LICENSE_SERVER  default 27005@license.seas.harvard.edu
#                          (from COMSOL63/Multiphysics/license/license.dat)
#   LMUTIL                 path to lmutil(.exe) if it is not found automatically
#
# The server is on the Harvard network: off campus, connect to the VPN first.
#
# Exit status: 0 = all listed features have a free seat,
#              2 = at least one listed feature is fully in use,
#              1 = lmutil not found or server unreachable.

set -u

SERVER="${COMSOL_LICENSE_SERVER:-27005@license.seas.harvard.edu}"
# Features relevant to the OMC models: the core, LiveLink for MATLAB, the
# mechanical studies, and the optical (ewfd) studies, which can run on a
# Wave Optics or an RF seat.
DEFAULT_FEATURES=(COMSOL COMSOLGUI LLMATLAB STRUCTURALMECHANICS MEMS ACOUSTICS WAVEOPTICS RF)
FEATURES=()
SHOW_ALL=0
SHOW_USERS=0

usage() { sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts ":af:s:uh" opt; do
    case "$opt" in
        a) SHOW_ALL=1 ;;
        f) FEATURES+=("$(printf '%s' "$OPTARG" | tr '[:lower:]' '[:upper:]')") ;;
        s) SERVER="$OPTARG" ;;
        u) SHOW_USERS=1 ;;
        h) usage 0 ;;
        :) echo "option -$OPTARG needs a value" >&2; usage 1 ;;
        *) echo "unknown option -$OPTARG" >&2; usage 1 ;;
    esac
done
[ ${#FEATURES[@]} -eq 0 ] && FEATURES=("${DEFAULT_FEATURES[@]}")

# ---- locate lmutil (WSL / Windows install, Linux, macOS, or PATH) ----------
find_lmutil() {
    if [ -n "${LMUTIL:-}" ]; then printf '%s\n' "$LMUTIL"; return; fi
    local c
    for c in \
        /mnt/c/Program\ Files/COMSOL/COMSOL*/Multiphysics/license/win64/lmutil.exe \
        /c/Program\ Files/COMSOL/COMSOL*/Multiphysics/license/win64/lmutil.exe \
        /usr/local/comsol*/multiphysics/license/glnxa64/lmutil \
        /opt/comsol*/multiphysics/license/glnxa64/lmutil \
        /Applications/COMSOL*/Multiphysics/license/mac*64/lmutil; do
        [ -x "$c" ] && { printf '%s\n' "$c"; return; }
    done
    command -v lmutil 2>/dev/null || command -v lmutil.exe 2>/dev/null
}

LM="$(find_lmutil | tail -n 1)"
if [ -z "$LM" ] || [ ! -x "$LM" ]; then
    echo "lmutil not found. It ships with COMSOL under" >&2
    echo "  <COMSOL>/Multiphysics/license/<arch>/lmutil[.exe]; set LMUTIL=/path/to/lmutil" >&2
    exit 1
fi

# ---- query the server (strip Windows CRs) ----------------------------------
OUT="$("$LM" lmstat -a -c "$SERVER" 2>&1 | tr -d '\r')"

if ! printf '%s\n' "$OUT" | grep -q "license server UP"; then
    echo "Licence server $SERVER is not reachable or not up." >&2
    echo "Off campus? Connect to the Harvard VPN. lmutil said:" >&2
    printf '%s\n' "$OUT" | grep -iE "error|cannot|down|refused|timed out" | head -5 >&2
    exit 1
fi

printf 'Server: %s  (%s)\n' "$SERVER" \
    "$(printf '%s\n' "$OUT" | grep -m1 "license server UP" | sed 's/^ *//')"
printf 'Time:   %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S')"
printf '%-24s %6s %6s %6s  %s\n' FEATURE TOTAL USED FREE STATUS
printf '%-24s %6s %6s %6s  %s\n' -------  ----- ---- ---- ------

# Lines look like:
#   Users of COMSOL:  (Total of 13 licenses issued;  Total of 5 licenses in use)
#   Users of SERIAL:  (Uncounted, node-locked)
any_full=0
print_row() {
    local line="$1" feat total used free status
    feat="$(printf '%s' "$line" | sed -E 's/^Users of ([^:]+):.*/\1/')"
    if printf '%s' "$line" | grep -q "Uncounted"; then
        printf '%-24s %6s %6s %6s  %s\n' "$feat" - - - "uncounted"
        return
    fi
    total="$(printf '%s' "$line" | sed -E 's/.*Total of ([0-9]+) licenses? issued.*/\1/')"
    used="$(printf '%s' "$line" | sed -E 's/.*Total of ([0-9]+) licenses? in use.*/\1/')"
    free=$((total - used))
    if [ "$free" -le 0 ]; then status="FULL"; any_full=1; else status="ok"; fi
    printf '%-24s %6s %6s %6s  %s\n' "$feat" "$total" "$used" "$free" "$status"
}

if [ "$SHOW_ALL" -eq 1 ]; then
    while IFS= read -r line; do print_row "$line"; done \
        < <(printf '%s\n' "$OUT" | grep "^Users of ")
else
    for f in "${FEATURES[@]}"; do
        line="$(printf '%s\n' "$OUT" | grep -m1 "^Users of ${f}:")"
        if [ -z "$line" ]; then
            printf '%-24s %6s %6s %6s  %s\n' "$f" - - - "not on server"
        else
            print_row "$line"
        fi
    done
fi

# ---- optional: who holds the seats -----------------------------------------
if [ "$SHOW_USERS" -eq 1 ]; then
    echo
    echo "Seat holders (user host display version (server) start-time):"
    if [ "$SHOW_ALL" -eq 1 ]; then pat="^Users of "; else
        pat="^Users of ($(IFS='|'; echo "${FEATURES[*]}")):"; fi
    printf '%s\n' "$OUT" | awk -v pat="$pat" '
        /^Users of / { show = ($0 ~ pat); feat = $3; next }
        show && /^ +[^ ]+ [^ ]+ .*start/ { sub(/^ +/, ""); print "  " feat " " $0 }'
fi

if [ "$SHOW_ALL" -ne 1 ]; then
    echo
    echo "Note: the optical studies (ewfd) can use a WAVEOPTICS or an RF seat."
fi
exit $((any_full * 2))

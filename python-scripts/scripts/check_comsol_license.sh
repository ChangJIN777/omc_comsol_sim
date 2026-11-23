#!/usr/bin/env bash
# check_comsol_license.sh -- show free COMSOL licence seats on the Harvard
# FlexLM server and who is holding the seats in use, WITHOUT checking a
# licence out.
#
# check_comsol.py checks the licence by starting COMSOL via MPh, which holds a
# seat while it runs. This script only asks the server (lmutil lmstat), so it
# is safe to run at any time, e.g. before an MPh / LiveLink job.
#
# Usage:
#   scripts/check_comsol_license.sh                 # features these models use
#   scripts/check_comsol_license.sh -a              # every feature on the server
#   scripts/check_comsol_license.sh -f WAVEOPTICS   # one feature (repeatable)
#   scripts/check_comsol_license.sh -q              # table only, no seat holders
#   scripts/check_comsol_license.sh -s 27005@host   # another licence server
#
# Environment overrides:
#   COMSOL_LICENSE_SERVER  default 27005@license.seas.harvard.edu
#                          (from COMSOL63/Multiphysics/license/license.dat)
#   LMUTIL                 path to lmutil(.exe) if it is not found automatically
#
# Under each feature with seats in use, the holders are listed as
#   user  host  client-version  start-time  [xN if one checkout holds N seats]
# -q hides them (for scripting). -u, which used to turn the list on, is still
# accepted and does nothing.
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
SHOW_USERS=1

# Print the header comment: every '#' line after the shebang, up to the first
# line that is not a comment (so the header can change length freely).
usage() {
    awk 'NR == 1 { next } !/^#/ { exit } { sub(/^# ?/, ""); print }' "$0"
    exit "${1:-0}"
}

while getopts ":af:s:quh" opt; do
    case "$opt" in
        a) SHOW_ALL=1 ;;
        f) FEATURES+=("$(printf '%s' "$OPTARG" | tr '[:lower:]' '[:upper:]')") ;;
        s) SERVER="$OPTARG" ;;
        q) SHOW_USERS=0 ;;
        u) ;;   # holders are shown by default now; kept for old command lines
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

# ---- table, with the seat holders under each feature -----------------------
# lmstat prints one block per feature:
#   Users of COMSOL:  (Total of 13 licenses issued;  Total of 5 licenses in use)
#   Users of SERIAL:  (Uncounted, node-locked)
# and, only when seats are in use, after it:
#     "COMSOL" v6.4, vendor: LMCOMSOL, expiry: permanent(no expiration date)
#     floating license
#
#       user host display (v6.3) (10.1.2.3/27005 7815), start Mon 11/23 14:43
#       user host display (v6.4) (10.1.2.3/27005 115), start Mon 11/23 10:06, 2 licenses
# A holder line is recognised by its "(vX) (server/port handle), start" tail,
# so the quoted / vendor / "floating license" lines are never taken for one.
# Holders belong to the most recent "Users of" line. Plain POSIX awk (macOS).
printf '%s\n' "$OUT" | awk -v all="$SHOW_ALL" -v users="$SHOW_USERS" \
                           -v want="${FEATURES[*]}" '
    function row(f, t, u, fr, st) {
        printf "%-24s %6s %6s %6s  %s\n", f, t, u, fr, st
    }
    # Left-justify s in w columns (printf "%-*s" is not portable awk).
    function pad(s, w) { while (length(s) < w) s = s " "; return s }
    # Leading number of the first "Total of N license(s) <what>" in the line.
    function total_of(line, what,    s) {
        if (!match(line, "Total of [0-9]+ licenses? " what)) return ""
        s = substr(line, RSTART + 9, RLENGTH - 9)
        sub(/ .*/, "", s)
        return s + 0
    }
    # "host display" -> host. FlexLM does not delimit the two, so with spaces
    # in either the split is a guess: drop a tty/X display, or a display that
    # repeats the host (the usual Windows case); otherwise keep the lot.
    function host_of(rest,    n, a, h, i, half) {
        n = split(rest, a, " ")
        if (n <= 1) return rest
        if (n == 2) return a[1]
        if (a[n] ~ /^\// || a[n] ~ /^:[0-9]/ || a[n] ~ /:[0-9]+(\.[0-9]+)?$/) {
            h = a[1]
            for (i = 2; i < n; i++) h = h " " a[i]
            return h
        }
        if (n % 2 == 0) {
            half = substr(rest, 1, (length(rest) - 1) / 2)
            if (rest == half " " half) return half
        }
        return rest
    }
    function show(f,    i, k, fr, st) {
        if (!(f in seen)) { row(f, "-", "-", "-", "not on server"); return }
        if (uncounted[f]) row(f, "-", "-", "-", "uncounted")
        else {
            fr = total[f] - used[f]
            if (fr <= 0) { st = "FULL"; any_full = 1 } else st = "ok"
            row(f, total[f], used[f], fr, st)
        }
        if (!users) return
        for (i = 1; i <= nh[f]; i++) {
            k = f SUBSEP i
            printf "    %s  %s  %s  %s%s\n", pad(hu[k], wu), pad(hh[k], wh), \
                pad(hv[k], wv), "since " hs[k], (hn[k] > 1 ? "  x" hn[k] : "")
        }
    }

    /^Users of / {
        feat = $0
        sub(/^Users of /, "", feat); sub(/:.*/, "", feat)
        if (feat in seen) { feat = ""; next }       # keep the first block
        seen[feat] = 1; order[++nfeat] = feat; nh[feat] = 0
        if ($0 ~ /Uncounted/) uncounted[feat] = 1
        else { total[feat] = total_of($0, "issued") + 0
               used[feat]  = total_of($0, "in use") + 0 }
        next
    }
    feat != "" && /^[ \t]/ && match($0, / \(v[^)]*\) \([^)]*\), start /) {
        who  = substr($0, 1, RSTART - 1)
        tail = substr($0, RSTART + 1)               # (v6.3) (srv/port n), start ...
        sub(/^[ \t]+/, "", who)
        ver = tail; sub(/^\(/, "", ver); sub(/\).*/, "", ver)
        st = tail;  sub(/^.*\), start /, "", st)
        n = 1
        if (match(st, /, [0-9]+ licenses?$/)) {
            n = substr(st, RSTART + 2) + 0
            st = substr(st, 1, RSTART - 1)
        }
        user = who; sub(/ .*/, "", user)
        rest = who; sub(/^[^ ]+ */, "", rest)
        k = feat SUBSEP (++nh[feat])
        hu[k] = user; hh[k] = host_of(rest); hv[k] = ver; hs[k] = st; hn[k] = n
    }

    END {
        if (all) nw = nfeat
        else     nw = split(want, wanted, " ")
        # common column widths over the holders that will be printed
        wu = 4; wh = 4; wv = 4
        for (j = 1; j <= nw; j++) {
            f = all ? order[j] : wanted[j]
            for (i = 1; i <= nh[f]; i++) {
                k = f SUBSEP i
                if (length(hu[k]) > wu) wu = length(hu[k])
                if (length(hh[k]) > wh) wh = length(hh[k])
                if (length(hv[k]) > wv) wv = length(hv[k])
            }
        }
        row("FEATURE", "TOTAL", "USED", "FREE", "STATUS")
        row("-------", "-----", "----", "----", "------")
        for (j = 1; j <= nw; j++) show(all ? order[j] : wanted[j])
        exit (any_full ? 2 : 0)
    }'
rc=$?

if [ "$SHOW_ALL" -ne 1 ]; then
    echo
    echo "Note: the optical studies (ewfd) can use a WAVEOPTICS or an RF seat."
fi
exit "$rc"

#!/usr/bin/env bash
# wait_for_comsol_seat.sh -- wait until a COMSOL licence seat is free, then
# run a command. Nothing is checked out while waiting: the wait only asks the
# licence server (check_comsol_license.sh -p, i.e. lmutil lmstat). A seat is
# taken by the command itself, and only for as long as the command runs.
#
# Usage:
#   scripts/wait_for_comsol_seat.sh [options] -- <command> [args...]
#   scripts/wait_for_comsol_seat.sh [options]         # only wait, then exit 0
#
# Options:
#   -f FEATURE        feature to watch (repeatable). Default: WAVEOPTICS RF,
#                     either of which can run the optical (ewfd) studies.
#   --batch           preset: watch WAVEOPTICSBATCH RFBATCH (for jobs run with
#                     run_comsol_batch.sh, which uses the batch licence pool)
#                     The rib template also takes MEMS[BATCH] and
#                     CADIMPORT[BATCH]; to wait for the whole set use e.g.
#                     -f WAVEOPTICSBATCH -f MEMSBATCH -f CADIMPORTBATCH --all
#   --any             start when ANY watched feature has a free seat (default)
#   --all             start only when EVERY watched feature has a free seat
#                     (e.g. -f STRUCTURALMECHANICS -f WAVEOPTICS --all)
#   -i SECONDS        poll interval, default 60, minimum 30 (the licence
#                     server is shared; do not poll it faster)
#   -t SECONDS        give up after waiting this long in total; 0 = never
#                     (default). Time spent running the command is not counted.
#   -r N              if the command fails AND its output shows a licence
#                     error (the seat was taken between the check and the
#                     checkout), go back to waiting and try again, at most N
#                     times. Default 0 = run the command once.
#   --retry-pattern RE  extended regex (case-insensitive) that marks a licence
#                     error in the command output; replaces the default
#   --retry-log FILE  also search FILE for the pattern (for commands that log
#                     to a file instead of the terminal)
#   -s PORT@HOST      licence server (passed on to check_comsol_license.sh)
#   -h, --help        this text
#
# Examples:
#   scripts/wait_for_comsol_seat.sh -r 3 -- python scripts/run_one.py --config configs/run_one_rib.yaml
#   scripts/wait_for_comsol_seat.sh --batch -r 2 -- scripts/run_comsol_batch.sh --study std2
#   scripts/wait_for_comsol_seat.sh -f WAVEOPTICS -t 3600 && echo "seat free"
#
# Status lines go to stderr with a time stamp; the command's own stdout and
# stderr are passed through unchanged (with -r they are also copied to a
# temporary file, which is searched for the licence error and then deleted).
#
# Ctrl-C while waiting stops at once. Ctrl-C while the command runs goes to
# the command (it is in the foreground); the wrapper then exits too, without
# retrying.
#
# Exit status: the command's exit status, or
#   0    no command given and a seat is free
#   124  timed out waiting (-t)
#   130  interrupted while waiting
#   1    the licence server could not be queried 5 times in a row
#   64   usage error
#
# Environment:
#   COMSOL_LICENSE_CHECK    the query script (default: check_comsol_license.sh
#                           next to this file). It is called with -p and must
#                           print "FEATURE TOTAL USED FREE STATUS" lines.
#   COMSOL_LICENSE_SERVER, LMUTIL   see check_comsol_license.sh
#   WAIT_SEAT_MIN_INTERVAL  lowest accepted -i (default 30; for tests only)

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHECK="${COMSOL_LICENSE_CHECK:-$HERE/check_comsol_license.sh}"
MIN_INTERVAL="${WAIT_SEAT_MIN_INTERVAL:-30}"
MAX_CHECK_FAILURES=5
# FlexLM / COMSOL wording for "no seat": error -4 is "Licensed number of users
# already reached". A bare "-4" would match far too much, so it is only
# accepted in the FlexLM forms "(-4,132)" and "error:-4,132".
DEFAULT_PATTERN='licen[sc]e error|licensed number of users already reached|flexnet licensing error|\(-4,[0-9]+|error: ?-4,|could not (obtain|check ?out) licen[sc]e|licen[sc]e.*(not available|unavailable)'

FEATURES=()
PRESET_BATCH=0
MODE=any
INTERVAL=60
TIMEOUT=0
RETRIES=0
PATTERN="$DEFAULT_PATTERN"
RETRY_LOG=""
SERVER_ARGS=()
SLEEP_PID=""
CAPTURE=""

usage() {
    awk 'NR == 1 { next } !/^#/ { exit } { sub(/^# ?/, ""); print }' "$0"
    exit "${1:-0}"
}

log() { printf '[%s] wait_for_comsol_seat: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }

die_usage() { echo "wait_for_comsol_seat: $*" >&2; echo "Try --help." >&2; exit 64; }

is_uint() { case "$1" in ''|*[!0-9]*) return 1 ;; *) return 0 ;; esac; }

need_value() { [ $# -ge 2 ] || die_usage "option $1 needs a value"; }

while [ $# -gt 0 ]; do
    case "$1" in
        -f) need_value "$@"; FEATURES+=("$(printf '%s' "$2" | tr '[:lower:]' '[:upper:]')"); shift 2 ;;
        --batch) PRESET_BATCH=1; shift ;;
        --any) MODE=any; shift ;;
        --all) MODE=all; shift ;;
        -i) need_value "$@"; INTERVAL="$2"; shift 2 ;;
        -t) need_value "$@"; TIMEOUT="$2"; shift 2 ;;
        -r) need_value "$@"; RETRIES="$2"; shift 2 ;;
        --retry-pattern) need_value "$@"; PATTERN="$2"; shift 2 ;;
        --retry-log) need_value "$@"; RETRY_LOG="$2"; shift 2 ;;
        -s) need_value "$@"; SERVER_ARGS=(-s "$2"); shift 2 ;;
        -h|--help) usage 0 ;;
        --) shift; break ;;
        -*) die_usage "unknown option $1" ;;
        *) die_usage "unexpected argument '$1' (put the command after --)" ;;
    esac
done

is_uint "$INTERVAL" || die_usage "-i needs a whole number of seconds"
is_uint "$TIMEOUT"  || die_usage "-t needs a whole number of seconds"
is_uint "$RETRIES"  || die_usage "-r needs a whole number"
is_uint "$MIN_INTERVAL" || MIN_INTERVAL=30
if [ "$INTERVAL" -lt "$MIN_INTERVAL" ]; then
    log "interval $INTERVAL s is below the minimum; using $MIN_INTERVAL s"
    INTERVAL="$MIN_INTERVAL"
fi
if [ ${#FEATURES[@]} -eq 0 ]; then
    if [ "$PRESET_BATCH" -eq 1 ]; then FEATURES=(WAVEOPTICSBATCH RFBATCH)
    else FEATURES=(WAVEOPTICS RF); fi
fi
[ -f "$CHECK" ] || { echo "wait_for_comsol_seat: query script not found: $CHECK" >&2; exit 1; }

cleanup() {
    [ -n "$SLEEP_PID" ] && kill "$SLEEP_PID" 2>/dev/null
    [ -n "$CAPTURE" ] && rm -f "$CAPTURE"
}
on_signal() {
    log "interrupted; nothing was checked out by this script"
    cleanup
    exit 130
}
trap on_signal INT TERM
trap cleanup EXIT

# Sleep that a signal can interrupt at once (a plain foreground sleep would
# delay the trap until it returns).
nap() {
    sleep "$1" &
    SLEEP_PID=$!
    wait "$SLEEP_PID" 2>/dev/null
    SLEEP_PID=""
}

# Query the server once. Sets SUMMARY ("WAVEOPTICS 3/3, RF 0/3") and returns
#   0 = start now, 1 = no seat yet, 2 = the query failed.
seat_free() {
    local out rc feat total used free status n_free=0 n=0 parts=""
    local check_args=(-p)
    for feat in "${FEATURES[@]}"; do check_args+=(-f "$feat"); done
    [ ${#SERVER_ARGS[@]} -gt 0 ] && check_args+=("${SERVER_ARGS[@]}")
    out="$(bash "$CHECK" "${check_args[@]}" 2>&1)"
    rc=$?
    if [ "$rc" -ne 0 ] && [ "$rc" -ne 2 ]; then
        SUMMARY="$(printf '%s' "$out" | head -n 2 | tr '\n' ' ')"
        return 2
    fi
    while read -r feat total used free status; do
        [ -n "$feat" ] || continue
        n=$((n + 1))
        case "$status" in
            ok)        n_free=$((n_free + 1)); parts+="$feat $used/$total, " ;;
            uncounted) n_free=$((n_free + 1)); parts+="$feat uncounted, " ;;
            FULL)      parts+="$feat $used/$total, " ;;
            missing)   parts+="$feat not on server, " ;;
            *)         SUMMARY="unexpected line from $CHECK: $feat $total $used $free $status"
                       return 2 ;;
        esac
    done <<EOF
$out
EOF
    SUMMARY="${parts%, }"
    if [ "$n" -eq 0 ]; then SUMMARY="no feature lines from $CHECK"; return 2; fi
    if [ "$MODE" = all ]; then [ "$n_free" -eq "$n" ]; else [ "$n_free" -gt 0 ]; fi
}

WAITED=0            # seconds spent waiting so far (all rounds)

# Poll until a seat is free. Returns 0, or exits on timeout / query failure.
wait_for_seat() {
    local failures=0 rc
    while :; do
        seat_free; rc=$?
        if [ "$rc" -eq 0 ]; then
            log "seat free ($SUMMARY)"
            return 0
        elif [ "$rc" -eq 2 ]; then
            failures=$((failures + 1))
            log "licence query failed ($failures/$MAX_CHECK_FAILURES): $SUMMARY"
            if [ "$failures" -ge "$MAX_CHECK_FAILURES" ]; then
                log "giving up: the licence server cannot be queried"
                exit 1
            fi
        else
            failures=0
            log "no free seat ($SUMMARY); mode $MODE; next check in $INTERVAL s"
        fi
        if [ "$TIMEOUT" -gt 0 ] && [ $((WAITED + INTERVAL)) -gt "$TIMEOUT" ]; then
            log "timed out after waiting $WAITED s (limit $TIMEOUT s)"
            exit 124
        fi
        nap "$INTERVAL"
        WAITED=$((WAITED + INTERVAL))
    done
}

log "watching ${FEATURES[*]} (mode $MODE, every $INTERVAL s, timeout ${TIMEOUT} s, retries $RETRIES)"

if [ $# -eq 0 ]; then
    wait_for_seat
    exit 0
fi

attempt=0
while :; do
    wait_for_seat
    attempt=$((attempt + 1))
    log "running (attempt $attempt): $*"
    if [ "$RETRIES" -gt 0 ]; then
        CAPTURE="$(mktemp "${TMPDIR:-/tmp}/wait_seat.XXXXXX")"
        "$@" > >(tee -a "$CAPTURE") 2> >(tee -a "$CAPTURE" >&2)
        rc=$?
        sleep 0.3           # let the two tee processes flush
    else
        "$@"
        rc=$?
    fi
    if [ "$rc" -eq 0 ]; then
        log "command finished (exit 0)"
        exit 0
    fi
    licence_error=0
    if [ "$RETRIES" -gt 0 ]; then
        grep -Eiq -- "$PATTERN" "$CAPTURE" 2>/dev/null && licence_error=1
        if [ -n "$RETRY_LOG" ] && [ -f "$RETRY_LOG" ]; then
            grep -Eiq -- "$PATTERN" "$RETRY_LOG" 2>/dev/null && licence_error=1
        fi
        rm -f "$CAPTURE"; CAPTURE=""
    fi
    if [ "$rc" -eq 130 ] || [ "$rc" -eq 143 ]; then
        log "command was interrupted (exit $rc); not retrying"
        exit "$rc"
    fi
    if [ "$licence_error" -eq 1 ] && [ "$attempt" -le "$RETRIES" ]; then
        log "command failed (exit $rc) with a licence error: the seat was taken first. Retry $attempt/$RETRIES after $INTERVAL s"
        nap "$INTERVAL"
        WAITED=$((WAITED + INTERVAL))
        continue
    fi
    if [ "$licence_error" -eq 1 ]; then
        log "command failed (exit $rc) with a licence error; no retries left"
    else
        log "command failed (exit $rc)"
    fi
    exit "$rc"
done

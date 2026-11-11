#!/usr/bin/env bash
# Background worker: works through queue/jobs.txt on a few low-priority cores, forever.
# Start it with `mise run queue:install` (launchd agent) or `mise run queue:run` (foreground).
# Settings: copy queue/config.example to queue/config.local.
set -uo pipefail
cd "$(dirname "$0")/.."

THREADS=2          # cores to use
BACKGROUND=1       # 1: macOS background QoS (efficiency cores, quiet fans) via taskpolicy; else nice
PULL=1             # git pull before each job, so jobs pushed to the repo arrive by themselves
PUSH=0             # 1: also push finished runs
POLL=60            # seconds between checks for new results while a job runs
IDLE_SLEEP=900     # seconds to wait when the queue is empty
NOTIFY=1           # 0: no desktop notifications
[ -f queue/config.local ] && . queue/config.local

export NOTIFY
mkdir -p queue/state
log() { echo "$(date '+%F %T') $*"; }
py=python3; [ -x .venv/bin/python ] && py=.venv/bin/python
if [ "$BACKGROUND" = 1 ] && command -v taskpolicy >/dev/null 2>&1; then runner=(taskpolicy -b); else runner=(nice -n 19); fi

pid=""
trap '[ -n "$pid" ] && kill "$pid" 2>/dev/null; rm -f queue/state/current; log "worker stopped"; exit 0' TERM INT

sync_repo() {
  [ "$PULL" = 1 ] || return 0
  git pull --rebase --autostash -q origin main 2>>queue/state/git.log || log "git pull failed (see queue/state/git.log)"
}

build() {
  if command -v cargo >/dev/null 2>&1; then cargo build --release -q
  elif command -v mise >/dev/null 2>&1; then mise exec -- cargo build --release -q
  else log "no cargo on PATH"; return 1; fi
}

finished() { grep -q '^done ' "runs/$1/run.log" 2>/dev/null; }

# first job in the queue that is not finished
next_job() {
  local line name
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line%%#*}"
    [ -z "${line//[[:space:]]/}" ] && continue
    # shellcheck disable=SC2086
    name="$(scripts/runname.sh $line)"
    if ! finished "$name"; then echo "$line"; return 0; fi
  done < queue/jobs.txt
  return 1
}

log "worker started: threads=$THREADS runner=${runner[*]} pull=$PULL push=$PUSH"
while true; do
  sync_repo
  if job="$(next_job)"; then
    if ! build; then log "build failed"; sleep "$IDLE_SLEEP"; continue; fi
    # shellcheck disable=SC2086
    name="$(scripts/runname.sh $job)"
    # shellcheck disable=SC2086
    set -- $job
    depth="$1"; any="$2"; shift 2
    log "start $name ($job)"
    echo "$name" > queue/state/current
    "${runner[@]}" target/release/webb-toroids --depth "$depth" --any "$any" "$@" --threads "$THREADS" --resume --out "runs/$name" >> queue/state/search.out 2>&1 &
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
      sleep "$POLL" & wait $! 2>/dev/null
      "$py" scripts/queue_notify.py "runs/$name" 2>>queue/state/notify.log || true
    done
    wait "$pid"; rc=$?; pid=""
    rm -f queue/state/current
    if finished "$name"; then
      "$py" scripts/queue_notify.py "runs/$name" --done 2>>queue/state/notify.log || true
      n=$(($(wc -l < "runs/$name/results.tsv") - 1))
      log "finished $name: $n results"
      git add "runs/$name" 2>>queue/state/git.log && git commit -q -m "queue: $name complete ($n results)" 2>>queue/state/git.log
      if [ "$PUSH" = 1 ]; then sync_repo; git push -q origin main 2>>queue/state/git.log || log "git push failed (see queue/state/git.log)"; fi
    else
      log "$name stopped early (exit $rc); will resume it in 10 minutes"
      sleep 600 & wait $! 2>/dev/null
    fi
  else
    sleep "$IDLE_SLEEP" & wait $! 2>/dev/null
  fi
done

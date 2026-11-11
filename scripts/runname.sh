#!/usr/bin/env bash
# Print the run directory name for a set of search arguments: DEPTH ANY [flags].
# Used by `mise run search` and by the queue worker, so both agree on where a job's output is.
set -euo pipefail
depth="${1:-2}"
any="${2:-0}"
if [ "$#" -ge 2 ]; then shift 2; else shift "$#"; fi
name="d${depth}a${any}"
# flags that change what is searched become part of the directory name
prev=""
for a in "$@"; do
  case "$prev" in
    --min-group) name="${name}-g${a}" ;;
    --near) name="${name}-near${a}" ;;
    --links) if [ "$a" = struts ]; then name="${name}-struts"; else name="${name}-l$(printf '%s' "$a" | cksum | cut -d' ' -f1)"; fi ;;
    --max-faces) name="${name}-f${a}" ;;
    --max-radius) name="${name}-r${a}" ;;
    --max-pieces) name="${name}-p${a}" ;;
    --only|--seeds|--seed-range) name="${name}-$(printf '%s' "$a" | cksum | cut -d' ' -f1)" ;;
  esac
  prev="$a"
done
echo "$name"

#!/bin/bash
set -euo pipefail

# Delete RSS/Base *.res files referenced in RELAX/IT*/group_* lists.
# Usage examples:
#   ./cleanup_rss_base.sh                 # dry-run, shows targets (up to 200 entries)
#   ./cleanup_rss_base.sh --apply         # delete after confirmation
#   ./cleanup_rss_base.sh --apply --yes   # delete without prompt
#   ./cleanup_rss_base.sh --save-list list.txt  # keep full target list for review

script_dir=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1
  pwd
)
repo_root=$(cd "$script_dir/.." && pwd)
rss_base="$repo_root/RSS/Base"
dry_run=true
assume_yes=false
save_list=""

usage() {
  cat <<EOF
Usage: $(basename "$0") [--apply] [--yes] [--save-list <path>]
  --apply   perform deletion (defaults to dry-run)
  --yes     skip confirmation prompt (implies --apply)
  --save-list <path>  write full target list to <path>
  --help    show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply)
      dry_run=false
      shift
      ;;
    --yes)
      dry_run=false
      assume_yes=true
      shift
      ;;
    --save-list)
      if [[ $# -lt 2 ]]; then
        echo "--save-list requires a path argument" >&2
        usage
        exit 1
      fi
      save_list=$2
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ ! -d "$rss_base" ]]; then
  echo "RSS/Base directory not found at: $rss_base" >&2
  exit 1
fi

relax_dir="$repo_root/RELAX"
first_group=$(find "$relax_dir" -type f -name 'group_*' -print -quit 2>/dev/null || true)
if [[ -z "$first_group" ]]; then
  echo "No group_* files found under RELAX/. Nothing to do."
  exit 0
fi

targets_file=$(mktemp)
trap 'rm -f "$targets_file"' EXIT

python3 - "$rss_base" "$relax_dir" <<'PY' >"$targets_file"
import os
import sys

rss_base = os.path.realpath(sys.argv[1])
relax_dir = os.path.realpath(sys.argv[2])

seen = set()
warnings = 0

def handle_entry(entry, origin):
    entry = entry.rstrip("\r\n")
    if not entry:
        return
    if entry.startswith("/"):
        cand = entry
    else:
        cand = os.path.join(rss_base, entry)
    real = os.path.realpath(cand)
    if not os.path.isfile(real):
        print(f"Warning: file missing -> {cand} (referenced in {origin})", file=sys.stderr)
        return
    if real == rss_base or real.startswith(rss_base + os.sep):
        if real not in seen:
            seen.add(real)
            print(real)
    else:
        print(f"Warning: skipping path outside RSS/Base -> {real} (from {origin})", file=sys.stderr)

for root, _, files in os.walk(relax_dir):
    for name in files:
        if name.startswith("group_"):
            group_path = os.path.join(root, name)
            try:
                with open(group_path, "r", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        handle_entry(line, group_path)
            except FileNotFoundError:
                print(f"Warning: group file missing -> {group_path}", file=sys.stderr)
            except Exception as exc:
                print(f"Warning: failed to read {group_path}: {exc}", file=sys.stderr)
PY

if [[ ! -s "$targets_file" ]]; then
  echo "No existing RSS/Base files referenced in group_* lists."
  exit 0
fi

if [[ -n "$save_list" ]]; then
  cp "$targets_file" "$save_list"
  echo "Full target list saved to: $save_list"
fi

target_count=$(wc -l < "$targets_file" | tr -d '[:space:]')

echo "Targets ($target_count):"
display_limit=200
if (( target_count <= display_limit )); then
  sed 's/^/  /' "$targets_file"
else
  head -n "$display_limit" "$targets_file" | sed 's/^/  /'
  echo "  ..."
  echo "  (list truncated, use --save-list <path> to keep the full list)"
fi

if $dry_run; then
  echo "Dry-run mode: no files deleted. Re-run with --apply to delete."
  exit 0
fi

if ! $assume_yes; then
  read -rp "Delete all listed files from RSS/Base? [y/N] " reply
  case "$reply" in
    y|Y|yes|YES) ;;
    *)
      echo "Aborted."
      exit 0
      ;;
  esac
fi

removed=0
while IFS= read -r file; do
  if rm -f -- "$file"; then
    ((removed++))
  fi
done < "$targets_file"

echo "Removed $removed file(s) from RSS/Base."

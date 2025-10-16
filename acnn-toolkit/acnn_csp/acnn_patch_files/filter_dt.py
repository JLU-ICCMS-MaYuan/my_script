#!/usr/bin/env python3
"""Filter structures based on emodel error report and move outliers.

Typical workflow in POT/IT*/:
    emodel model-restart/model-100000 DT/ > ss 2>&1
    python filter_dt.py 0.20 0.50 0.60

The script reads the generated `ss` file, compares per-structure e/f/v errors
with supplied thresholds, and moves offending structures from `DT/` into
`trash_bin/` so they can be inspected separately.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Move DT structures whose e/f/v errors exceed thresholds."
    )
    parser.add_argument(
        "e_threshold",
        type=float,
        help="Maximum allowed per-atom energy error.",
    )
    parser.add_argument(
        "f_threshold",
        type=float,
        help="Maximum allowed per-atom force error.",
    )
    parser.add_argument(
        "v_threshold",
        type=float,
        help="Maximum allowed per-atom virial error.",
    )
    parser.add_argument(
        "--ss-file",
        default="ss",
        type=Path,
        help="Path to the emodel output file (default: ss).",
    )
    parser.add_argument(
        "--dt-dir",
        default=Path("DT"),
        type=Path,
        help="Path to the DT directory containing structures.",
    )
    parser.add_argument(
        "--trash-bin",
        default=Path("trash_bin"),
        type=Path,
        help="Destination for filtered structures (default: trash_bin).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only report moves without touching the filesystem.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print diagnostics about resolved paths.",
    )
    return parser.parse_args()


def parse_ss_lines(lines: Iterable[str]) -> Iterator[Tuple[str, float, float, float]]:
    """Yield (structure_name, e_err, f_err, v_err) from an emodel report."""
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 5:
            continue
        name = parts[0]
        # Skip obvious summary lines.
        if name.lower().startswith("rmse") or name.lower().startswith("analysis"):
            continue
        try:
            e_err = float(parts[-3])
            f_err = float(parts[-2])
            v_err = float(parts[-1])
        except ValueError:
            continue
        yield name, e_err, f_err, v_err


def build_name_index(dt_dir: Path) -> Dict[str, List[Path]]:
    """Build a lookup of filename -> paths inside DT."""
    index: Dict[str, List[Path]] = {}
    for path in dt_dir.rglob("*"):
        if path == dt_dir:
            continue
        key_candidates = {path.name}
        stem = path.stem
        if stem:
            key_candidates.add(stem)
        for key in key_candidates:
            index.setdefault(key, []).append(path)
    return index


def resolve_structure_path(
    dt_dir: Path,
    name: str,
    index_cache: Dict[str, Optional[Dict[str, List[Path]]]],
    verbose: bool = False,
) -> Optional[Path]:
    """Return the best-guess path for a structure inside DT."""
    direct = dt_dir / name
    if direct.exists():
        return direct

    stem_candidate = dt_dir / Path(name).stem
    if stem_candidate.exists():
        return stem_candidate

    if index_cache.get("value") is None:
        if verbose:
            print(f"Indexing DT directory: {dt_dir}")
        index_cache["value"] = build_name_index(dt_dir)

    index = index_cache["value"] or {}

    matches: Sequence[Path] = index.get(name) or index.get(Path(name).stem) or []
    if not matches:
        return None

    if verbose and len(matches) > 1:
        print(
            f"Multiple matches for {name}: choosing {matches[0]!s}; "
            "please adjust manually if needed."
        )
    return matches[0]


def ensure_unique_destination(dest: Path) -> Path:
    """Append a numeric suffix if the destination already exists."""
    if not dest.exists():
        return dest
    counter = 1
    while True:
        candidate = dest.with_name(f"{dest.name}__{counter}")
        if not candidate.exists():
            return candidate
        counter += 1


def move_structure(
    src: Path, dt_root: Path, trash_root: Path, dry_run: bool = False
) -> Path:
    """Move a structure relative to the DT root into the trash bin."""
    try:
        rel_path = src.relative_to(dt_root)
    except ValueError:
        rel_path = Path(src.name)
    destination = ensure_unique_destination(trash_root / rel_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return destination
    shutil.move(str(src), str(destination))
    return destination


def main() -> None:
    args = parse_args()

    ss_path: Path = args.ss_file
    dt_dir: Path = args.dt_dir
    trash_dir: Path = args.trash_bin

    if not ss_path.exists():
        raise SystemExit(f"ss file not found: {ss_path}")
    if not dt_dir.exists():
        raise SystemExit(f"DT directory not found: {dt_dir}")

    trash_dir.mkdir(parents=True, exist_ok=True)
    ss_path = ss_path.resolve()
    dt_dir = dt_dir.resolve()
    trash_dir = trash_dir.resolve()

    exceeded: List[Tuple[str, float, float, float]] = []
    seen_names = set()

    with ss_path.open() as handle:
        for name, e_err, f_err, v_err in parse_ss_lines(handle):
            if (
                e_err > args.e_threshold
                or f_err > args.f_threshold
                or v_err > args.v_threshold
            ):
                if name in seen_names:
                    continue
                exceeded.append((name, e_err, f_err, v_err))
                seen_names.add(name)

    if not exceeded:
        print("No structures exceed the specified thresholds.")
        return

    index_cache: Dict[str, Optional[Dict[str, List[Path]]]] = {"value": None}
    moved = []
    missing = []

    for idx, (name, e_err, f_err, v_err) in enumerate(exceeded, start=1):
        resolved = resolve_structure_path(
            dt_dir, name, index_cache, verbose=args.verbose
        )
        if not resolved:
            print(
                f"[{idx}/{len(exceeded)}] {name}: exceeds thresholds "
                f"(e={e_err:.4e}, f={f_err:.4e}, v={v_err:.4e}) but not found in {dt_dir}"
            )
            missing.append(name)
            continue

        dest = move_structure(resolved, dt_dir, trash_dir, dry_run=args.dry_run)
        status = "would move to" if args.dry_run else "moved to"
        print(
            f"[{idx}/{len(exceeded)}] {name}: {status} {dest} "
            f"(e={e_err:.4e}, f={f_err:.4e}, v={v_err:.4e})"
        )
        moved.append(dest)

    if args.dry_run:
        print(f"\nDry run complete. {len(moved)} structures would be relocated.")
    else:
        print(f"\nFinished. {len(moved)} structures moved to {trash_dir}.")

    if missing:
        print(
            f"\nCould not locate {len(missing)} structures inside {dt_dir}. "
            "Please verify their paths manually."
        )


if __name__ == "__main__":
    main()

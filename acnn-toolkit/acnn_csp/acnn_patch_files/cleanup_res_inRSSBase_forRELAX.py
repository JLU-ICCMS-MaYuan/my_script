#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

BAR_WIDTH = 40


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Remove RSS/Base files referenced by RELAX group_* lists."
    )
    parser.add_argument("--apply", action="store_true", help="delete the files")
    parser.add_argument("--yes", action="store_true", help="skip confirmation prompt")
    parser.add_argument(
        "--save-list", metavar="PATH", help="write the list of targets to a file"
    )
    parser.add_argument(
        "--keep-it",
        metavar="N",
        type=int,
        action="append",
        default=[],
        help="skip deleting RSS/Base IT{N}/res entries (may be used multiple times)",
    )
    parser.add_argument(
        "--project-root",
        metavar="PATH",
        help="root directory that contains RELAX and RSS/Base subdirectories",
    )
    parser.add_argument(
        "--relax-dir",
        metavar="PATH",
        help="explicit RELAX directory path (overrides project root detection)",
    )
    parser.add_argument(
        "--rss-base",
        metavar="PATH",
        help="explicit RSS/Base directory path (overrides project root detection)",
    )
    return parser.parse_args(argv)


def first_existing_directory(candidates):
    seen = set()
    for candidate in candidates:
        if candidate is None:
            continue
        path = Path(candidate)
        try:
            resolved = path.resolve()
        except OSError:
            continue
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_dir():
            return resolved
    return None


def resolve_directories(args):
    script_dir = Path(__file__).resolve().parent
    cwd = Path.cwd()
    repo_root = None
    if args.project_root:
        candidate = Path(args.project_root).expanduser()
        if not candidate.is_dir():
            raise ValueError(f"Project root not found: {candidate}")
        repo_root = candidate.resolve()

    relax_dir = None
    if args.relax_dir:
        candidate = Path(args.relax_dir).expanduser()
        if not candidate.is_dir():
            raise ValueError(f"RELAX directory not found: {candidate}")
        relax_dir = candidate.resolve()

    rss_base = None
    if args.rss_base:
        candidate = Path(args.rss_base).expanduser()
        if not candidate.is_dir():
            raise ValueError(f"RSS/Base directory not found: {candidate}")
        rss_base = candidate.resolve()

    if repo_root:
        if relax_dir is None:
            candidate = repo_root / "RELAX"
            if candidate.is_dir():
                relax_dir = candidate.resolve()
        if rss_base is None:
            candidate = repo_root / "RSS" / "Base"
            if candidate.is_dir():
                rss_base = candidate.resolve()

    if relax_dir is None:
        relax_candidates = [
            cwd / "RELAX",
            script_dir / "RELAX",
            script_dir.parent / "RELAX",
        ]
        for base in [cwd, *cwd.parents]:
            relax_candidates.append(base / "RELAX")
        relax_dir = first_existing_directory(relax_candidates)

    if rss_base is None:
        rss_candidates = []
        if relax_dir is not None:
            rss_candidates.append(relax_dir.parent / "RSS" / "Base")
        rss_candidates.extend(
            [
                cwd / "RSS" / "Base",
                script_dir / "RSS" / "Base",
                script_dir.parent / "RSS" / "Base",
            ]
        )
        if repo_root:
            rss_candidates.append(repo_root / "RSS" / "Base")
        for base in [cwd, *cwd.parents]:
            rss_candidates.append(base / "RSS" / "Base")
        rss_base = first_existing_directory(rss_candidates)

    return relax_dir, rss_base


def collect_targets(relax_dir, rss_base):
    relax_dir = relax_dir.resolve()
    rss_base = rss_base.resolve()
    seen = set()
    targets = []
    group_stats = []

    group_files = [path for path in sorted(relax_dir.rglob("group_*")) if path.is_file()]
    total_groups = len(group_files)

    for index, group_file in enumerate(group_files, 1):
        references = 0
        display_name = str(group_file.relative_to(relax_dir))
        try:
            lines = group_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        except OSError:
            group_stats.append({"path": group_file, "references": references})
            if total_groups:
                show_group_progress(index, total_groups, display_name)
            continue

        for line in lines:
            item = line.strip()
            if not item:
                continue
            path = Path(item)
            if not path.is_absolute():
                path = rss_base / item
            try:
                resolved = path.resolve(strict=True)
            except FileNotFoundError:
                continue
            try:
                resolved.relative_to(rss_base)
            except ValueError:
                continue
            references += 1
            if resolved not in seen:
                seen.add(resolved)
                targets.append(resolved)
        group_stats.append({"path": group_file, "references": references})
        if total_groups:
            show_group_progress(index, total_groups, display_name)
    return targets, group_stats


def filter_it_res_targets(targets, rss_base, keep_it_numbers):
    if not keep_it_numbers:
        return targets, []

    keep_set = set(keep_it_numbers)
    filtered = []
    skipped = []

    for path in targets:
        try:
            relative = path.relative_to(rss_base)
        except ValueError:
            filtered.append(path)
            continue

        parts = relative.parts
        if len(parts) >= 2 and parts[0].startswith("IT") and parts[1] == "res":
            suffix = parts[0][2:]
            if suffix.isdigit() and int(suffix) in keep_set:
                skipped.append(path)
                continue

        filtered.append(path)

    return filtered, skipped


def count_res_structures(rss_base):
    try:
        return sum(1 for _ in rss_base.rglob("*.res"))
    except OSError:
        return 0


def confirm_deletion():
    reply = input("Delete all listed files from RSS/Base? [y/N] ").strip().lower()
    return reply in {"y", "yes"}


def show_progress(current, total):
    filled = int(BAR_WIDTH * current / total)
    bar = "#" * filled + "-" * (BAR_WIDTH - filled)
    sys.stdout.write(f"\rDeleting [{bar}] {current}/{total}")
    if current == total:
        sys.stdout.write("\n")
    sys.stdout.flush()


def show_group_progress(current, total, group_name):
    filled = int(BAR_WIDTH * current / total)
    bar = "#" * filled + "-" * (BAR_WIDTH - filled)
    suffix = f" ({group_name})" if group_name else ""
    sys.stdout.write(f"\rScanning groups [{bar}] {current}/{total}{suffix}")
    if current == total:
        sys.stdout.write("\n")
    sys.stdout.flush()


def main(argv):
    args = parse_args(argv)
    try:
        relax_dir, rss_base = resolve_directories(args)
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1

    if relax_dir is None:
        print(
            "RELAX directory not found. Use --relax-dir or --project-root to specify it.",
            file=sys.stderr,
        )
        return 1

    if rss_base is None:
        print(
            "RSS/Base directory not found. Use --rss-base or --project-root to specify it.",
            file=sys.stderr,
        )
        return 1

    print(f"Using RELAX directory: {relax_dir}")
    print(f"Using RSS/Base directory: {rss_base}")

    targets, group_stats = collect_targets(relax_dir, rss_base)
    targets, skipped = filter_it_res_targets(targets, rss_base, args.keep_it)

    if group_stats:
        total_refs = sum(stat["references"] for stat in group_stats)
        print(f"Scanned {len(group_stats)} group file(s) with {total_refs} referenced entries.")
    else:
        print("No group_* files found under RELAX.")

    if skipped:
        print(f"Skipping {len(skipped)} file(s) in requested IT/res directories.")

    if not targets:
        print("No files to delete.")
        if args.apply or args.yes:
            remaining = count_res_structures(rss_base)
            print(f"RSS/Base structures remaining: {remaining}")
        return 0

    if args.save_list:
        Path(args.save_list).write_text("\n".join(str(path) for path in targets) + "\n")
        print("Target list saved.")

    print(f"Found {len(targets)} file(s) referenced.")

    apply_changes = args.apply or args.yes
    if not apply_changes:
        print("Dry-run: no files deleted.")
        return 0

    if not args.yes and not confirm_deletion():
        print("Aborted.")
        return 0

    deleted = 0
    failed = 0
    total = len(targets)

    for index, target in enumerate(targets, 1):
        try:
            target.unlink()
            deleted += 1
        except FileNotFoundError:
            failed += 1
        except OSError:
            failed += 1
        show_progress(index, total)

    print(f"Removed {deleted} file(s).", end="")
    if failed:
        print(f" Skipped {failed} file(s).")
    else:
        print()
    remaining = count_res_structures(rss_base)
    print(f"RSS/Base structures remaining: {remaining}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

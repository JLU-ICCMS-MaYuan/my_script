#!/usr/bin/env python3
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

DEFAULT_INPUT = "/home/mayuan/code/my_script/test/acnn/cam"
DEFAULT_OUTPUT_DIR = "/home/mayuan/code/my_script/test/acnn"

FORMULA_RE = re.compile(r"^([A-Z][a-z]?\d*)+$")
ELEMENT_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula_counts(formula: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for element, number_text in ELEMENT_RE.findall(formula):
        count = int(number_text) if number_text else 1
        counts[element] = counts.get(element, 0) + count
    return counts


def find_formula(tokens: List[str]) -> str | None:
    for token in tokens:
        if FORMULA_RE.match(token):
            return token
    return None


def normalize_enthalpy(value: float, raw: str) -> str:
    if abs(value) <= 1e-8:
        return "0"
    return raw


def parse_lines(lines: List[str]) -> Tuple[List[Tuple[int, str, Dict[str, int], float, str]], List[str]]:
    parsed: List[Tuple[int, str, Dict[str, int], float, str]] = []
    warnings: List[str] = []

    for index, raw_line in enumerate(lines, start=1):
        if not raw_line.strip():
            continue

        tokens = raw_line.strip().split()
        if len(tokens) < 7:
            warnings.append(f"第 {index} 行列数不足: {raw_line.rstrip()}")
            continue

        enthalpy_token = tokens[5]
        try:
            enthalpy_value = float(enthalpy_token)
        except ValueError:
            warnings.append(f"第 {index} 行第 6 列非数字: {enthalpy_token}")
            continue

        formula = find_formula(tokens[6:])
        if formula is None:
            warnings.append(f"第 {index} 行未找到化学式: {raw_line.rstrip()}")
            continue

        counts = parse_formula_counts(formula)
        enthalpy_text = normalize_enthalpy(enthalpy_value, enthalpy_token)
        parsed.append((index, formula, counts, enthalpy_value, enthalpy_text))

    return parsed, warnings


def main() -> int:
    input_path = Path(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT)
    output_dir = Path(sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUTPUT_DIR)

    if not input_path.exists():
        print(f"输入文件不存在: {input_path}", file=sys.stderr)
        return 1

    content = input_path.read_text(errors="ignore")
    lines = content.splitlines()

    parsed, warnings = parse_lines(lines)

    elements = sorted({element for _, _, counts, _, _ in parsed for element in counts.keys()})
    header = ["Number", "formula", *elements, "enthalpy"]

    stable_rows: List[str] = [",".join(header)]
    unstable_rows: List[str] = [",".join(header)]

    for number, formula, counts, enthalpy_value, enthalpy_text in parsed:
        row = [str(number), formula]
        row.extend(str(counts.get(element, 0)) for element in elements)
        row.append(enthalpy_text)
        line = ",".join(row)

        if abs(enthalpy_value) <= 1e-8:
            stable_rows.append(line)
        elif 0 < enthalpy_value <= 0.05:
            unstable_rows.append(line)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "stable.csv").write_text("\n".join(stable_rows) + "\n", encoding="utf-8")
    (output_dir / "unstable.csv").write_text("\n".join(unstable_rows) + "\n", encoding="utf-8")

    if warnings:
        print("\n".join(warnings), file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

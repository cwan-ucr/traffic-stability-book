#!/usr/bin/env python3
"""Prepare and merge a browser-based translation of a LaTeX manuscript.

This adapter exists for translators that accept only pasted text.  It keeps
mathematics, cross-references, citations, paths, and LaTeX command syntax out
of the text sent to the external service.  `extract` writes <= 4,300-character
jobs; `merge` reassembles their returned English prose into a new `.tex` file.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from translate_tex_deepl import HAN, TEXT_COMMAND, TOKEN, mask_opaque, unmask

UNIT = re.compile(r"@@UNIT(\d{6})@@")
SEP = re.compile(r"ZZZ_SEP_\d{6}_ZZZ")
SEP_ID = re.compile(r"ZZZ_SEP_(\d{6})_ZZZ")
WEB_KEEP = re.compile(r"KPRSV(\d{5})KPRSV")
FORMULA = re.compile(r"FORMULA_(\d{5})_TOKEN")


def split_unit(text: str, limit: int) -> list[str]:
    """Split an oversized prose unit only at readable sentence boundaries."""
    fragments: list[str] = []
    remaining = text
    while len(remaining) > limit:
        cut = max(remaining.rfind(mark, 0, limit) for mark in ("。", "！", "？", "；", "\n"))
        if cut < max(80, limit // 3):
            cut = limit
        else:
            cut += 1
        fragments.append(remaining[:cut])
        remaining = remaining[cut:]
    if remaining:
        fragments.append(remaining)
    return fragments


def collect_units(source: str, unit_limit: int) -> tuple[str, list[str], list[str]]:
    """Replace every Chinese prose span with a stable unit token."""
    masked, opaque = mask_opaque(source)
    # Google Translate treats leading @ characters as markup in some contexts.
    # Use a neutral, alphanumeric token while the text is outside LaTeX.
    masked = TOKEN.sub(lambda match: f"KPRSV{match.group(1)}KPRSV", masked)
    units: list[str] = []

    def add_unit(text: str) -> str:
        tokens: list[str] = []
        # Give the browser a readable, inert formula anchor.  It preserves
        # word order around mathematics; merge restores the original TeX.
        text = WEB_KEEP.sub(lambda match: f"FORMULA_{match.group(1)}_TOKEN", text)
        for fragment in split_unit(text, unit_limit):
            units.append(fragment)
            tokens.append(f"@@UNIT{len(units) - 1:06d}@@")
        return "".join(tokens)

    def protect_raw(text: str) -> str:
        # Retain every LaTeX control sequence (including `\\`, `\%`, and
        # `\&`), then turn only the intervening Chinese prose into units.
        # Table alignment marks and grouping braces are structural LaTeX too.
        # Keep them out of the browser request, otherwise web translators may
        # treat an entire table row as malformed input and omit individual cells.
        pieces = re.split(
            r"(\\(?:[A-Za-z@]+\*?(?:\[[^\]]*\])?|.)|[&{}])",
            text,
        )
        for index in range(0, len(pieces), 2):
            if HAN.search(pieces[index]):
                pieces[index] = add_unit(pieces[index])
        return "".join(pieces)

    # TEXT_COMMAND.sub would visit all headings before the intervening prose,
    # destroying reading order.  Walk the source instead, so each browser job
    # retains its original narrative context.
    output: list[str] = []
    position = 0
    for match in TEXT_COMMAND.finditer(masked):
        output.append(protect_raw(masked[position:match.start()]))
        body = match.group("body")
        output.append(match.group(0).replace(body, protect_raw(body) if HAN.search(body) else body, 1))
        position = match.end()
    output.append(protect_raw(masked[position:]))
    return "".join(output), units, opaque


def make_jobs(units: list[str], limit: int) -> list[dict[str, object]]:
    jobs: list[dict[str, object]] = []
    current: list[str] = []
    current_ids: list[int] = []
    current_length = 0
    for unit_id, text in enumerate(units):
        if len(text) > limit:
            raise ValueError(
                f"Translation unit {unit_id} has {len(text)} characters, above {limit}. "
                "Split that prose paragraph before preparing browser jobs."
            )
        separator = f"\n\nZZZ_SEP_{unit_id:06d}_ZZZ\n\n" if current else ""
        if current and current_length + len(separator) + len(text) > limit:
            jobs.append({"id": len(jobs), "unit_ids": current_ids, "text": "".join(current)})
            current, current_ids, current_length, separator = [], [], 0, ""
        current.append(separator + text)
        current_ids.append(unit_id)
        current_length += len(separator) + len(text)
    if current:
        jobs.append({"id": len(jobs), "unit_ids": current_ids, "text": "".join(current)})
    return jobs


def extract(source_path: Path, state_path: Path, jobs_path: Path, limit: int) -> None:
    source = source_path.read_text(encoding="utf-8")
    stage, units, opaque = collect_units(source, limit)
    jobs = make_jobs(units, limit)
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(json.dumps({"template": stage, "opaque": opaque, "units": units}, ensure_ascii=False), encoding="utf-8")
    jobs_path.write_text(json.dumps(jobs, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"units": len(units), "jobs": len(jobs), "han": len(HAN.findall(source))}))


def merge(state_path: Path, translations_path: Path, output_path: Path) -> None:
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if translations_path.is_dir():
        result_files = sorted(translations_path.glob("results-*.json"))
        if not result_files:
            raise ValueError("No results-*.json files found in translation directory.")
        results = [job for path in result_files for job in json.loads(path.read_text(encoding="utf-8"))]
    else:
        results = json.loads(translations_path.read_text(encoding="utf-8"))
    translated_units: dict[int, str] = {}
    for job in results:
        ids = [int(value) for value in job["unit_ids"]]
        parts = SEP.split(job["translation"])
        separators = SEP_ID.findall(job["translation"])
        if separators != [f"{value:06d}" for value in ids[1:]] or len(parts) != len(ids):
            raise ValueError(f"Job {job['id']} lost or altered a separator; do not merge it.")
        translated_units.update(dict(zip(ids, parts)))
    expected = set(range(len(state["units"])))
    if set(translated_units) != expected:
        raise ValueError("Incomplete translation results; every unit must appear exactly once.")
    assembled = UNIT.sub(lambda match: translated_units[int(match.group(1))], state["template"])
    assembled = FORMULA.sub(lambda match: f"KPRSV{match.group(1)}KPRSV", assembled)
    assembled = WEB_KEEP.sub(lambda match: f"@@KEEP{match.group(1)}@@", assembled)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(unmask(assembled, state["opaque"]), encoding="utf-8")
    print(json.dumps({"output": str(output_path), "han_after": len(HAN.findall(assembled))}))


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    extract_parser = subparsers.add_parser("extract")
    extract_parser.add_argument("source", type=Path)
    extract_parser.add_argument("state", type=Path)
    extract_parser.add_argument("jobs", type=Path)
    extract_parser.add_argument("--limit", type=int, default=4300)
    merge_parser = subparsers.add_parser("merge")
    merge_parser.add_argument("state", type=Path)
    merge_parser.add_argument("translations", type=Path)
    merge_parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "extract":
        extract(args.source, args.state, args.jobs, args.limit)
    else:
        merge(args.state, args.translations, args.output)


if __name__ == "__main__":
    main()

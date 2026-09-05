#!/usr/bin/env python3
"""Translate visible Chinese prose in a LaTeX manuscript through DeepL.

The script deliberately preserves mathematical environments, labels, references,
citations, paths, and verbatim blocks. It is a first-pass translator only: the
generated manuscript must be compiled and terminology-reviewed chapter by
chapter. Set DEEPL_AUTH_KEY in the environment; do not put an API key in git.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


MATH_ENVS = r"equation\*?|align\*?|gather\*?|multline\*?|cases|matrix|bmatrix|pmatrix|vmatrix|Vmatrix|verbatim|lstlisting"
OPAQUE = re.compile(
    rf"(\\begin\{{(?:{MATH_ENVS})\}}.*?\\end\{{(?:{MATH_ENVS})\}}|"
    r"(?<!\\)\\\[.*?(?<!\\)\\\]|\\\(.*?\\\)|\$\$.*?\$\$|(?<!\\)\$.*?(?<!\\)\$|"
    r"\\(?:label|ref|eqref|cite|citep|citet|includegraphics|input|url|path|verb)\*?(?:\[[^\]]*\])?\{[^{}]*\})",
    re.DOTALL,
)
TEXT_COMMAND = re.compile(
    r"\\(?P<name>chapter|section|subsection|subsubsection|part|caption|captionof|"
    r"chapterguide|sectionstudy|sectionwrapup|figurenote|tabletext|textbf|emph|"
    r"cnemph|texorpdfstring)\*?(?:\[[^\]]*\])?\{(?P<body>[^{}]*)\}"
)
HAN = re.compile(r"[\u3400-\u9fff]")
TOKEN = re.compile(r"@@KEEP(\d{5})@@")


def mask_opaque(text: str) -> tuple[str, list[str]]:
    saved: list[str] = []

    def replace(match: re.Match[str]) -> str:
        saved.append(match.group(0))
        return f"@@KEEP{len(saved) - 1:05d}@@"

    return OPAQUE.sub(replace, text), saved


def unmask(text: str, saved: list[str]) -> str:
    return TOKEN.sub(lambda match: saved[int(match.group(1))], text)


def split_for_api(text: str, limit: int = 3800) -> list[str]:
    if len(text) <= limit:
        return [text]
    chunks, remaining = [], text
    while len(remaining) > limit:
        cut = max(remaining.rfind("。", 0, limit), remaining.rfind("\n", 0, limit), remaining.rfind(" ", 0, limit))
        if cut < limit // 2:
            cut = limit
        chunks.append(remaining[:cut + 1])
        remaining = remaining[cut + 1:]
    if remaining:
        chunks.append(remaining)
    return chunks


def deepl_translate(text: str, key: str, endpoint: str) -> str:
    payload = urlencode(
        {
            "auth_key": key,
            "text": text,
            "source_lang": "ZH",
            "target_lang": "EN-US",
            "formality": "prefer_more",
            "preserve_formatting": "1",
        }
    ).encode()
    request = Request(endpoint, data=payload, headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urlopen(request, timeout=90) as response:  # nosec B310: explicit official endpoint supplied by user
        payload = json.loads(response.read().decode("utf-8"))
    return payload["translations"][0]["text"]


def translate_plain(text: str, key: str | None, endpoint: str, dry_run: bool) -> str:
    if not HAN.search(text):
        return text
    if dry_run:
        return text
    if not key:
        raise RuntimeError("DEEPL_AUTH_KEY is required unless --dry-run is used.")
    return "".join(deepl_translate(chunk, key, endpoint) for chunk in split_for_api(text))


def translate_document(text: str, key: str | None, endpoint: str, dry_run: bool) -> str:
    masked, saved = mask_opaque(text)

    def translate_command(match: re.Match[str]) -> str:
        body = translate_plain(match.group("body"), key, endpoint, dry_run)
        return match.group(0).replace(match.group("body"), body, 1)

    masked = TEXT_COMMAND.sub(translate_command, masked)
    # Remaining prose lies between TeX commands. Keep commands and their optional
    # arguments verbatim, and translate only the inter-command text spans.
    pieces = re.split(r"(\\[A-Za-z@]+\*?(?:\[[^\]]*\])?)", masked)
    for index in range(0, len(pieces), 2):
        pieces[index] = translate_plain(pieces[index], key, endpoint, dry_run)
    return unmask("".join(pieces), saved)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--dry-run", action="store_true", help="validate segmentation without calling DeepL")
    parser.add_argument("--endpoint", default=os.environ.get("DEEPL_ENDPOINT", "https://api.deepl.com/v2/translate"))
    args = parser.parse_args()
    source = args.source.read_text(encoding="utf-8")
    translated = translate_document(source, os.environ.get("DEEPL_AUTH_KEY"), args.endpoint, args.dry_run)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(translated, encoding="utf-8")
    print(json.dumps({"source": str(args.source), "output": str(args.output), "dry_run": args.dry_run, "han_before": len(HAN.findall(source)), "han_after": len(HAN.findall(translated))}))


if __name__ == "__main__":
    main()

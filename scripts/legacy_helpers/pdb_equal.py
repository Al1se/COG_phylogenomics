#!/usr/bin/env python3
"""Compare two files with comma/newline separated PDB_chain tokens."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


TOKEN_SPLIT_RE = re.compile(r"[,\n\r\t ]+")


def read_tokens(path: Path) -> set[str]:
    content = path.read_text(encoding="utf-8")
    return {token.strip() for token in TOKEN_SPLIT_RE.split(content) if token.strip()}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Diff PDB_chain token sets")
    parser.add_argument("file_a", help="First token file")
    parser.add_argument("file_b", help="Second token file")
    parser.add_argument(
        "--mode",
        choices=["b_minus_a", "a_minus_b", "symmetric"],
        default="b_minus_a",
        help="Diff mode (default: b_minus_a)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    set_a = read_tokens(Path(args.file_a))
    set_b = read_tokens(Path(args.file_b))

    if args.mode == "b_minus_a":
        result = set_b - set_a
    elif args.mode == "a_minus_b":
        result = set_a - set_b
    else:
        result = set_a ^ set_b

    for token in sorted(result):
        print(token)


if __name__ == "__main__":
    main()

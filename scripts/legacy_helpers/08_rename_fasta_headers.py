#!/usr/bin/env python3
"""Bulk rename FASTA headers to '<ID>@<suffix>' format."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator


Record = tuple[str, str]


def iter_fasta(path: Path) -> Iterator[Record]:
    header = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def write_fasta(path: Path, records: list[Record], width: int = 80) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def header_id(header: str) -> str:
    token = header.split()[0]
    if "|" in token:
        token = token.split("|", 1)[0]
    return token


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Rename FASTA headers to ID@suffix")
    parser.add_argument("input_fasta", help="Input FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA")
    parser.add_argument(
        "--suffix",
        default="phylogeny",
        help="Suffix after '@' (default: phylogeny)",
    )
    parser.add_argument(
        "--mapping",
        default=None,
        help="Optional mapping TSV: id<TAB>suffix; overrides --suffix per ID",
    )
    return parser


def load_mapping(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            mapping[parts[0]] = parts[1]
    return mapping


def main() -> None:
    args = build_parser().parse_args()

    records = list(iter_fasta(Path(args.input_fasta)))
    if not records:
        raise SystemExit("Input FASTA is empty.")

    mapping = load_mapping(Path(args.mapping)) if args.mapping else {}

    renamed: list[Record] = []
    for header, seq in records:
        rec_id = header_id(header)
        suffix = mapping.get(rec_id, args.suffix)
        renamed.append((f"{rec_id}@{suffix}", seq))

    write_fasta(Path(args.output), renamed)

    print(f"Input records: {len(records)}")
    print(f"Output records: {len(renamed)}")
    print(f"Output: {args.output}")
    if args.mapping:
        print(f"Mapping used: {args.mapping}")


if __name__ == "__main__":
    main()


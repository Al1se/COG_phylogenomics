#!/usr/bin/env python3
"""Reorder FASTA records by the leaf order in a Newick tree."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterator


Record = tuple[str, str]
LEAF_RE = re.compile(r"(?<=[(,])\s*('(?:[^']|'')*'|[^'():;,\s]+)\s*(?=[:),])")


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


def fasta_id(header: str) -> str:
    return header.split()[0]


def normalize_newick_label(label: str) -> str:
    label = label.strip()
    if label.startswith("'") and label.endswith("'"):
        label = label[1:-1].replace("''", "'")
    return label


def parse_newick_leaf_order(newick_text: str) -> list[str]:
    labels = [normalize_newick_label(m.group(1)) for m in LEAF_RE.finditer(newick_text)]
    deduped = list(dict.fromkeys(labels))
    return deduped


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Sort FASTA records by Newick leaf order")
    parser.add_argument("input_fasta", help="Input FASTA")
    parser.add_argument("tree_newick", help="Input Newick tree")
    parser.add_argument("-o", "--output", required=True, help="Output reordered FASTA")
    parser.add_argument(
        "--drop-unmatched",
        action="store_true",
        help="Drop FASTA records that are missing in the tree",
    )
    parser.add_argument(
        "--allow-partial-match",
        action="store_true",
        help="Allow tree/FASTA mismatches instead of failing (default: fail on mismatch)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    records = list(iter_fasta(Path(args.input_fasta)))
    if not records:
        raise SystemExit("Input FASTA is empty.")

    leaf_order = parse_newick_leaf_order(Path(args.tree_newick).read_text(encoding="utf-8"))
    if not leaf_order:
        raise SystemExit("No leaf labels found in Newick tree.")

    key_to_index: dict[str, int] = {}
    for idx, (header, _seq) in enumerate(records):
        keys = {header, fasta_id(header)}
        for key in keys:
            if key in key_to_index:
                continue
            key_to_index[key] = idx

    used_indices: set[int] = set()
    ordered: list[Record] = []
    missing_in_fasta: list[str] = []

    for leaf in leaf_order:
        idx = key_to_index.get(leaf)
        if idx is None:
            missing_in_fasta.append(leaf)
            continue
        if idx in used_indices:
            continue
        ordered.append(records[idx])
        used_indices.add(idx)

    not_in_tree = [record for i, record in enumerate(records) if i not in used_indices]

    if not args.allow_partial_match and (missing_in_fasta or not_in_tree):
        raise SystemExit(
            "Tree/FASTA mismatch detected. "
            f"Tree leaves missing in FASTA: {len(missing_in_fasta)}; "
            f"FASTA records missing in tree: {len(not_in_tree)}. "
            "Use --allow-partial-match only for exploratory runs."
        )

    if not args.drop_unmatched:
        ordered.extend(not_in_tree)

    write_fasta(Path(args.output), ordered)

    print(f"Input FASTA records: {len(records)}")
    print(f"Tree leaves parsed: {len(leaf_order)}")
    print(f"Matched records: {len(used_indices)}")
    print(f"Tree leaves missing in FASTA: {len(missing_in_fasta)}")
    print(f"FASTA records missing in tree: {len(not_in_tree)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Map residue contacts from structure numbering to alignment columns."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def iter_fasta(path: Path):
    header = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
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


def build_seq_to_col_map(aligned_seq: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    seq_pos = 0
    for col_idx, aa in enumerate(aligned_seq, start=1):
        if aa != "-":
            seq_pos += 1
            mapping[seq_pos] = col_idx
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(description="Map label_seq_id contacts onto MSA columns")
    parser.add_argument("contacts_tsv", help="TSV from analyze_metal_contacts.py")
    parser.add_argument("msa_fasta", help="Aligned FASTA with tokens present as first header token")
    parser.add_argument("-o", "--output", required=True, help="Output mapped TSV")
    args = parser.parse_args()

    msa = {header.split()[0]: seq for header, seq in iter_fasta(Path(args.msa_fasta))}
    maps = {token: build_seq_to_col_map(seq) for token, seq in msa.items()}

    with Path(args.contacts_tsv).open("r", encoding="utf-8") as src, Path(args.output).open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = list(reader.fieldnames or []) + ["alignment_column"]
        writer = csv.DictWriter(dst, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            token = row["token"]
            label_seq_id = row["label_seq_id"]
            alignment_column = ""
            if token in maps and label_seq_id not in {"", ".", "?"}:
                try:
                    alignment_column = str(maps[token].get(int(label_seq_id), ""))
                except ValueError:
                    alignment_column = ""
            row["alignment_column"] = alignment_column
            writer.writerow(row)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

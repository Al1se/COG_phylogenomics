#!/usr/bin/env python3
"""Map structure contacts onto alignment columns using a token->MSA sequence map."""

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
    parser = argparse.ArgumentParser(description="Map contact residues to MSA columns via token->MSA matches")
    parser.add_argument("contacts_tsv", help="TSV from analyze_metal_contacts.py")
    parser.add_argument("token_map_tsv", help="TSV from map_structure_sequences_to_alignment.py")
    parser.add_argument("msa_fasta", help="Aligned FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    args = parser.parse_args()

    msa = {header.split()[0]: seq for header, seq in iter_fasta(Path(args.msa_fasta))}
    msa_maps = {msa_id: build_seq_to_col_map(seq) for msa_id, seq in msa.items()}

    token_map: dict[str, dict[str, str]] = {}
    with Path(args.token_map_tsv).open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            token_map[row["token"]] = row

    with Path(args.contacts_tsv).open("r", encoding="utf-8") as src, Path(args.output).open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = list(reader.fieldnames or []) + ["matched_msa_id", "match_type", "alignment_column"]
        writer = csv.DictWriter(dst, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            token = row["token"]
            mapped = token_map.get(token, {})
            msa_id = mapped.get("representative_msa_id", "")
            match_type = mapped.get("match_type", "")
            alignment_column = ""
            label_seq_id = row.get("label_seq_id", "")

            if msa_id in msa_maps and label_seq_id not in {"", ".", "?"}:
                try:
                    alignment_column = str(msa_maps[msa_id].get(int(label_seq_id), ""))
                except ValueError:
                    alignment_column = ""

            row["matched_msa_id"] = msa_id
            row["match_type"] = match_type
            row["alignment_column"] = alignment_column
            writer.writerow(row)

    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Extract a compact subalignment by matching ungapped sequences to a source MSA."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterator


Record = tuple[str, str]


def iter_fasta(path: Path) -> Iterator[Record]:
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


def write_fasta(path: Path, records: list[Record], width: int = 80) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def fasta_id(header: str) -> str:
    return header.split()[0]


def normalize_sequence(seq: str) -> str:
    seq = seq.replace("-", "").replace(".", "").upper()
    return "".join(ch for ch in seq if "A" <= ch <= "Z")


def trim_all_gap_columns(records: list[Record]) -> list[Record]:
    if not records:
        return records
    seqs = [seq for _, seq in records]
    keep_cols = [
        idx
        for idx in range(len(seqs[0]))
        if any(seq[idx] != "-" for seq in seqs)
    ]
    trimmed: list[Record] = []
    for header, seq in records:
        trimmed_seq = "".join(seq[idx] for idx in keep_cols)
        trimmed.append((header, trimmed_seq))
    return trimmed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract query sequences as a compact subalignment")
    parser.add_argument("source_msa", help="Source aligned FASTA")
    parser.add_argument("query_fasta", help="Query FASTA with ungapped sequences to recover from the MSA")
    parser.add_argument("-o", "--output", required=True, help="Output aligned FASTA")
    parser.add_argument("--report", default=None, help="Optional TSV mapping report")
    parser.add_argument(
        "--trim-all-gap-columns",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Drop columns that are all-gap after extraction (default: on)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    source_records = list(iter_fasta(Path(args.source_msa)))
    if not source_records:
        raise SystemExit("Source MSA is empty.")
    source_by_seq: dict[str, list[tuple[str, str]]] = {}
    for header, aligned_seq in source_records:
        source_by_seq.setdefault(normalize_sequence(aligned_seq), []).append((header, aligned_seq))

    output_records: list[Record] = []
    report_rows: list[dict[str, str]] = []

    for header, seq in iter_fasta(Path(args.query_fasta)):
        rec_id = fasta_id(header)
        norm = normalize_sequence(seq)
        hits = source_by_seq.get(norm, [])
        row = {
            "query_id": rec_id,
            "query_length": str(len(norm)),
            "status": "missing",
            "matched_source_ids": "",
        }
        if hits:
            matched_ids = [fasta_id(src_header) for src_header, _ in hits]
            row["status"] = "exact"
            row["matched_source_ids"] = ";".join(matched_ids)
            output_records.append((rec_id, hits[0][1]))
        report_rows.append(row)

    if args.trim_all_gap_columns:
        output_records = trim_all_gap_columns(output_records)

    write_fasta(Path(args.output), output_records)

    if args.report:
        with Path(args.report).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["query_id", "query_length", "status", "matched_source_ids"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(report_rows)

    exact = sum(row["status"] == "exact" for row in report_rows)
    missing = sum(row["status"] == "missing" for row in report_rows)
    print(f"Queries: {len(report_rows)}")
    print(f"Exact matches: {exact}")
    print(f"Missing: {missing}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

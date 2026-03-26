#!/usr/bin/env python3
"""Filter FASTA records by length, motif presence, and optional sequence deduplication."""

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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filter FASTA records by simple QC rules")
    parser.add_argument("input_fasta", help="Input FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum length to keep")
    parser.add_argument("--max-length", type=int, default=0, help="Maximum length to keep (0 = no limit)")
    parser.add_argument(
        "--require-motif",
        default="",
        help="Require this motif to be present in the cleaned sequence",
    )
    parser.add_argument(
        "--dedup-by-sequence",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Drop exact duplicate cleaned sequences",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Optional TSV report",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    records = list(iter_fasta(Path(args.input_fasta)))
    if not records:
        raise SystemExit("Input FASTA is empty.")

    seen_seqs: set[str] = set()
    kept: list[Record] = []
    report_rows: list[dict[str, str]] = []

    for header, raw_seq in records:
        rec_id = fasta_id(header)
        seq = normalize_sequence(raw_seq)
        reason = "ok"
        keep = True
        motif_pos = ""

        if len(seq) < args.min_length:
            keep = False
            reason = "too_short"
        elif args.max_length and len(seq) > args.max_length:
            keep = False
            reason = "too_long"
        elif args.require_motif:
            idx = seq.find(args.require_motif)
            if idx < 0:
                keep = False
                reason = "missing_motif"
            else:
                motif_pos = str(idx + 1)
        if keep and args.dedup_by_sequence and seq in seen_seqs:
            keep = False
            reason = "duplicate_sequence"

        if keep:
            seen_seqs.add(seq)
            kept.append((rec_id, seq))

        report_rows.append(
            {
                "id": rec_id,
                "length": str(len(seq)),
                "motif_pos": motif_pos,
                "kept": str(int(keep)),
                "reason": reason,
            }
        )

    write_fasta(Path(args.output), kept)

    if args.report:
        with Path(args.report).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["id", "length", "motif_pos", "kept", "reason"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(report_rows)

    print(f"Input records: {len(records)}")
    print(f"Kept records: {len(kept)}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

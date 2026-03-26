#!/usr/bin/env python3
"""Prepare structure-derived FASTA before joint MSA with COG sequences."""

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


def fasta_id(header: str) -> str:
    return header.split()[0]


def normalize_sequence(seq: str) -> str:
    # Keep letters only, drop alignment gap chars and other noise.
    seq = seq.replace("-", "").replace(".", "").upper()
    return "".join(ch for ch in seq if "A" <= ch <= "Z")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Normalize and pre-filter structure FASTA (ungap + length + dedup)"
    )
    parser.add_argument("input_fasta", help="Input structure FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output cleaned FASTA")
    parser.add_argument(
        "--min-length",
        type=int,
        default=400,
        help="Minimum ungapped AA length to keep (default: 400)",
    )
    parser.add_argument(
        "--dedup-by-sequence",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Drop identical sequences (default: off)",
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

    kept: list[Record] = []
    report_lines: list[str] = ["id\torig_len\tclean_len\tstatus\treason"]

    seen_ids: set[str] = set()
    seen_seqs: set[str] = set()

    for header, seq in records:
        rec_id = fasta_id(header)
        clean = normalize_sequence(seq)
        orig_len = len(seq)
        clean_len = len(clean)

        status = "drop"
        reason = ""

        if rec_id in seen_ids:
            reason = "duplicate_id"
        elif clean_len < args.min_length:
            reason = "too_short"
        elif not clean:
            reason = "empty_after_cleanup"
        elif args.dedup_by_sequence and clean in seen_seqs:
            reason = "duplicate_sequence"
        else:
            status = "keep"
            reason = "ok"
            seen_ids.add(rec_id)
            seen_seqs.add(clean)
            kept.append((rec_id, clean))

        report_lines.append(f"{rec_id}\t{orig_len}\t{clean_len}\t{status}\t{reason}")

    write_fasta(Path(args.output), kept)

    if args.report:
        Path(args.report).write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Input records: {len(records)}")
    print(f"Kept records: {len(kept)}")
    print(f"Min length: {args.min_length}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

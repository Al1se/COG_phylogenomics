#!/usr/bin/env python3
"""Filter structure-only MSA records by per-sequence gap fraction."""

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


def ungap(seq: str) -> str:
    return seq.replace("-", "").replace(".", "")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Keep structure MSA sequences with acceptable per-sequence gap fraction"
    )
    parser.add_argument("input_msa", help="Aligned FASTA with only structure-derived records")
    parser.add_argument("-o", "--output", required=True, help="Output ungapped FASTA of kept records")
    parser.add_argument(
        "--max-gap-fraction",
        type=float,
        default=0.5,
        help="Maximum allowed per-sequence gap fraction in the alignment (default: 0.5)",
    )
    parser.add_argument(
        "--min-ungapped-length",
        type=int,
        default=400,
        help="Minimum ungapped length to keep after filtering (default: 400)",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Optional TSV report path",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    records = list(iter_fasta(Path(args.input_msa)))
    if not records:
        raise SystemExit("Input MSA is empty.")

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise SystemExit("Input FASTA is not an alignment: sequence lengths differ.")

    aln_len = next(iter(lengths))
    kept: list[Record] = []
    report_lines = [
        "id\taligned_len\tgap_count\tgap_fraction\tungapped_len\tkept\treason"
    ]

    for header, aligned_seq in records:
        rec_id = fasta_id(header)
        gap_count = aligned_seq.count("-") + aligned_seq.count(".")
        gap_fraction = gap_count / aln_len if aln_len else 0.0
        ungapped_seq = ungap(aligned_seq)
        ungapped_len = len(ungapped_seq)

        keep = True
        reason = "ok"
        if gap_fraction > args.max_gap_fraction:
            keep = False
            reason = "too_many_gaps"
        elif ungapped_len < args.min_ungapped_length:
            keep = False
            reason = "too_short_after_ungap"

        if keep:
            kept.append((rec_id, ungapped_seq))

        report_lines.append(
            f"{rec_id}\t{aln_len}\t{gap_count}\t{gap_fraction:.4f}\t{ungapped_len}\t{int(keep)}\t{reason}"
        )

    write_fasta(Path(args.output), kept)

    if args.report:
        Path(args.report).write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Input records: {len(records)}")
    print(f"Kept after gap filtering: {len(kept)}")
    print(f"Aligned length: {aln_len}")
    print(f"Max gap fraction: {args.max_gap_fraction}")
    print(f"Min ungapped length: {args.min_ungapped_length}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

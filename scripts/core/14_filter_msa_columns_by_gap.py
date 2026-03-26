#!/usr/bin/env python3
"""Trim alignment columns by minimum occupancy."""

from __future__ import annotations

import argparse
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Keep alignment columns with sufficient occupancy")
    parser.add_argument("input_msa", help="Aligned FASTA input")
    parser.add_argument("-o", "--output", required=True, help="Output aligned FASTA")
    parser.add_argument(
        "--min-occupancy",
        type=float,
        default=0.5,
        help="Minimum non-gap fraction required to keep a column (default: 0.5)",
    )
    parser.add_argument("--report", default=None, help="Optional text report path")
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
    n_records = len(records)
    keep_cols: list[int] = []
    occupancy_values: list[float] = []

    for idx in range(aln_len):
        non_gap = sum(seq[idx] not in "-." for _, seq in records)
        occupancy = non_gap / n_records
        occupancy_values.append(occupancy)
        if occupancy >= args.min_occupancy:
            keep_cols.append(idx)

    trimmed: list[Record] = []
    for header, seq in records:
        trimmed_seq = "".join(seq[idx] for idx in keep_cols)
        trimmed.append((header, trimmed_seq))

    write_fasta(Path(args.output), trimmed)

    if args.report:
        kept = len(keep_cols)
        dropped = aln_len - kept
        min_occ = min(occupancy_values) if occupancy_values else 0.0
        max_occ = max(occupancy_values) if occupancy_values else 0.0
        mean_occ = sum(occupancy_values) / len(occupancy_values) if occupancy_values else 0.0
        report = [
            f"input_records\t{n_records}",
            f"input_aligned_len\t{aln_len}",
            f"min_occupancy_threshold\t{args.min_occupancy}",
            f"kept_columns\t{kept}",
            f"dropped_columns\t{dropped}",
            f"output_aligned_len\t{kept}",
            f"min_column_occupancy\t{min_occ:.6f}",
            f"mean_column_occupancy\t{mean_occ:.6f}",
            f"max_column_occupancy\t{max_occ:.6f}",
        ]
        Path(args.report).write_text("\n".join(report) + "\n", encoding="utf-8")

    print(f"Input records: {n_records}")
    print(f"Input aligned length: {aln_len}")
    print(f"Min occupancy threshold: {args.min_occupancy}")
    print(f"Kept columns: {len(keep_cols)}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

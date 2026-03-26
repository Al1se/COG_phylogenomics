#!/usr/bin/env python3
"""Pick one or several representative query sequences from a COG FASTA file."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from statistics import median


def iter_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def normalize_sequence(sequence: str) -> str:
    sequence = sequence.replace("-", "").replace(".", "").upper()
    return "".join(ch for ch in sequence if "A" <= ch <= "Z")


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def dedup_by_sequence(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    seen: set[str] = set()
    deduped: list[tuple[str, str]] = []
    for header, sequence in records:
        if sequence in seen:
            continue
        seen.add(sequence)
        deduped.append((header, sequence))
    return deduped


def select_by_length_spread(
    records: list[tuple[str, str]],
    max_queries: int,
) -> list[tuple[int, str, str]]:
    if not records or max_queries <= 0:
        return []

    lengths = [len(seq) for _, seq in records]
    median_len = float(median(lengths))
    ranked = sorted(
        enumerate(records),
        key=lambda item: (
            abs(len(item[1][1]) - median_len),
            -len(item[1][1]),
            item[1][0],
        ),
    )
    return ranked[:max_queries]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Select representative query sequence(s) from an existing COG FASTA. "
            "The default strategy keeps unique sequences and prefers lengths closest "
            "to the median family length."
        )
    )
    parser.add_argument("input_fasta", help="Input FASTA for the current COG analysis set")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA with query sequence(s)")
    parser.add_argument("--report", required=True, help="Output TSV report")
    parser.add_argument(
        "--max-queries",
        type=int,
        default=3,
        help="Maximum number of representative query sequences to emit (default: 3)",
    )
    parser.add_argument(
        "--motif",
        default="",
        help="Optional preferred motif; if provided, motif-containing sequences are preferred",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    motif = args.motif.strip().upper()

    raw_records = [(header, normalize_sequence(seq)) for header, seq in iter_fasta_records(Path(args.input_fasta))]
    raw_records = [(header, seq) for header, seq in raw_records if seq]
    if not raw_records:
        raise SystemExit("No usable sequences found in input FASTA.")

    deduped = dedup_by_sequence(raw_records)
    motif_records = []
    for header, seq in deduped:
        if motif and motif in seq:
            motif_records.append((header, seq))

    pool = motif_records if motif_records else deduped
    selected = select_by_length_spread(pool, args.max_queries)
    if len(selected) < args.max_queries and pool is not deduped:
        taken_headers = {header for _, (header, _seq) in selected}
        leftovers = [(h, s) for h, s in deduped if h not in taken_headers]
        selected.extend(select_by_length_spread(leftovers, args.max_queries - len(selected)))

    lengths = [len(seq) for _, seq in deduped]
    median_len = float(median(lengths))

    out_lines: list[str] = []
    report_rows: list[dict[str, str]] = []
    for rank, (_original_idx, (header, seq)) in enumerate(selected, start=1):
        out_lines.append(f">{header}\n{wrap_sequence(seq)}")
        report_rows.append(
            {
                "rank": str(rank),
                "header": header,
                "length": str(len(seq)),
                "contains_motif": "yes" if motif and motif in seq else "no",
                "distance_from_median": f"{abs(len(seq) - median_len):.1f}",
            }
        )

    Path(args.output).write_text("\n".join(out_lines) + ("\n" if out_lines else ""), encoding="utf-8")
    with Path(args.report).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "header", "length", "contains_motif", "distance_from_median"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(report_rows)

    print(f"Input records: {len(raw_records)}")
    print(f"Unique sequences: {len(deduped)}")
    print(f"Selected queries: {len(report_rows)}")
    print(f"Representative FASTA: {args.output}")
    print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

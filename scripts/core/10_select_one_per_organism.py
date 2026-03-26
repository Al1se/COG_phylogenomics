#!/usr/bin/env python3
"""Select one representative sequence per organism from a FASTA file."""

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


def organism_key(header: str) -> str:
    token = fasta_id(header)
    return token.split("|", 1)[1] if "|" in token else token


def normalize_sequence(seq: str) -> str:
    seq = seq.replace("-", "").replace(".", "").upper()
    return "".join(ch for ch in seq if "A" <= ch <= "Z")


def motif_pos(seq: str, motif: str) -> int:
    idx = seq.find(motif)
    return idx + 1 if idx >= 0 else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Select one representative sequence per organism")
    parser.add_argument("input_fasta", help="Input FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA")
    parser.add_argument(
        "--motif",
        default="NADFDGD",
        help="Preferred motif for ranking (default: NADFDGD)",
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

    grouped: dict[str, list[tuple[str, str, str]]] = {}
    for header, seq in records:
        token = fasta_id(header)
        clean = normalize_sequence(seq)
        grouped.setdefault(organism_key(header), []).append((token, header, clean))

    chosen: list[Record] = []
    report_rows: list[dict[str, str]] = []

    for organism, candidates in sorted(grouped.items()):
        scored = []
        for token, header, seq in candidates:
            pos = motif_pos(seq, args.motif)
            score = (
                1 if pos else 0,
                len(seq),
                -pos if pos else 0,
                token,
            )
            scored.append((score, token, header, seq, pos))
        scored.sort(reverse=True)
        _, token, header, seq, pos = scored[0]
        chosen.append((token, seq))
        report_rows.append(
            {
                "organism": organism,
                "chosen_id": token,
                "length": str(len(seq)),
                "motif_pos": str(pos) if pos else "",
                "candidate_count": str(len(candidates)),
                "candidate_ids": ";".join(item[1] for item in scored),
            }
        )

    write_fasta(Path(args.output), chosen)

    if args.report:
        with Path(args.report).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "organism",
                    "chosen_id",
                    "length",
                    "motif_pos",
                    "candidate_count",
                    "candidate_ids",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(report_rows)

    print(f"Input records: {len(records)}")
    print(f"Unique organisms: {len(grouped)}")
    print(f"Chosen representatives: {len(chosen)}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

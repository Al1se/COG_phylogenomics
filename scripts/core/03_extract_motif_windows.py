#!/usr/bin/env python3
"""Extract motif-centered windows from FASTA sequences."""

from __future__ import annotations

import argparse
import re
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


def motif_to_regex(motif: str, wildcard: str) -> str:
    parts: list[str] = []
    wildcard = wildcard.upper()
    for ch in motif.upper():
        if ch == wildcard:
            parts.append(".")
        else:
            parts.append(re.escape(ch))
    return "".join(parts)


def extract_window(seq: str, start: int, end: int, left: int, right: int) -> tuple[str, int, int]:
    wanted_start = start - left
    wanted_end = end + right

    real_start = max(0, wanted_start)
    real_end = min(len(seq), wanted_end)

    left_pad = "-" * (real_start - wanted_start)
    right_pad = "-" * (wanted_end - real_end)

    window = left_pad + seq[real_start:real_end] + right_pad
    return window, start + 1, end


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract motif windows from FASTA")
    parser.add_argument("input_fasta", help="Input FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA with motif windows")
    parser.add_argument("--motif", required=True, help="Motif pattern (e.g. NADFDGD or XXAA)")
    parser.add_argument("--left", type=int, default=2, help="Left flank length (default: 2)")
    parser.add_argument("--right", type=int, default=2, help="Right flank length (default: 2)")
    parser.add_argument(
        "--wildcard",
        default="X",
        help="Wildcard character inside motif pattern (default: X)",
    )
    parser.add_argument(
        "--all-matches",
        action="store_true",
        help="Export all motif matches per sequence (default: first match only)",
    )
    parser.add_argument(
        "--strip-gaps",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Remove '-' before motif search (default: on)",
    )
    parser.add_argument(
        "--keep-unmatched",
        action="store_true",
        help="Keep sequences without motif as all-gap windows",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    records = list(iter_fasta(Path(args.input_fasta)))
    if not records:
        raise SystemExit("Input FASTA is empty.")

    motif_re = re.compile(motif_to_regex(args.motif, args.wildcard), flags=re.IGNORECASE)
    window_len = args.left + len(args.motif) + args.right

    out_records: list[Record] = []
    matched_sequences = 0
    total_hits = 0

    for header, seq in records:
        search_seq = seq.replace("-", "") if args.strip_gaps else seq
        hits = list(motif_re.finditer(search_seq))

        if not hits:
            if args.keep_unmatched:
                out_records.append((f"{header}|motif=NA", "-" * window_len))
            continue

        matched_sequences += 1

        selected_hits = hits if args.all_matches else hits[:1]
        for hit_index, hit in enumerate(selected_hits, start=1):
            window, motif_start, motif_end = extract_window(
                seq=search_seq,
                start=hit.start(),
                end=hit.end(),
                left=args.left,
                right=args.right,
            )
            total_hits += 1
            out_header = (
                f"{header}|motif={args.motif}|pos={motif_start}-{motif_end}|hit={hit_index}/{len(hits)}"
            )
            out_records.append((out_header, window))

    write_fasta(Path(args.output), out_records)

    print(f"Input records: {len(records)}")
    print(f"Sequences with motif: {matched_sequences}")
    print(f"Extracted windows: {total_hits}")
    print(f"Output records: {len(out_records)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

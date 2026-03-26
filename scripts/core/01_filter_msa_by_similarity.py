#!/usr/bin/env python3
"""Filter aligned FASTA records by similarity to reference profile columns."""

from __future__ import annotations

import argparse
from collections import Counter
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


def build_consensus(ref_seqs: list[str], ignore_chars: set[str]) -> str:
    length = len(ref_seqs[0])
    consensus_chars: list[str] = []

    for i in range(length):
        counts = Counter()
        for seq in ref_seqs:
            aa = seq[i]
            if aa == "-" or aa in ignore_chars:
                continue
            counts[aa] += 1

        if not counts:
            consensus_chars.append("-")
            continue

        # deterministic tie-break: by count desc, residue asc
        consensus_chars.append(sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0])

    return "".join(consensus_chars)


def similarity_to_profile(
    seq: str,
    profile: str,
    informative_cols: list[int],
    ignore_chars: set[str],
    penalize_gaps: bool,
) -> tuple[float, int, int]:
    matches = 0
    total = 0

    for i in informative_cols:
        ref_aa = profile[i]
        aa = seq[i]

        if aa in ignore_chars:
            continue

        if aa == "-" and not penalize_gaps:
            continue

        total += 1
        if aa == ref_aa:
            matches += 1

    if total == 0:
        return 0.0, matches, total

    return matches / total, matches, total


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filter MSA sequences by similarity to reference profile")
    parser.add_argument("input_msa", help="Aligned FASTA input")
    parser.add_argument("-o", "--output", required=True, help="Filtered FASTA output")
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.45,
        help="Minimum similarity to keep sequence (default: 0.45)",
    )
    parser.add_argument(
        "--reference-id",
        action="append",
        default=[],
        help="Reference sequence ID (can be repeated). If omitted, all sequences build the consensus profile",
    )
    parser.add_argument(
        "--ignore-chars",
        default="X?",
        help="Characters ignored in scoring (default: X?)",
    )
    parser.add_argument(
        "--penalize-gaps",
        action="store_true",
        help="Count gaps in target sequence as mismatches",
    )
    parser.add_argument(
        "--no-keep-reference",
        action="store_true",
        help="Do not force-keep reference IDs below threshold",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Optional TSV score report path",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    records = list(iter_fasta(Path(args.input_msa)))
    if not records:
        raise SystemExit("Input MSA is empty.")

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise SystemExit("Input FASTA is not an alignment: sequence lengths are different.")

    ignore_chars = set(args.ignore_chars)
    id_to_record = {fasta_id(header): (header, seq) for header, seq in records}

    ref_ids = list(dict.fromkeys(args.reference_id))
    if ref_ids:
        missing = [ref_id for ref_id in ref_ids if ref_id not in id_to_record]
        if missing:
            raise SystemExit(f"Reference IDs are missing in MSA: {', '.join(missing)}")
        ref_seqs = [id_to_record[ref_id][1] for ref_id in ref_ids]
    else:
        ref_seqs = [seq for _, seq in records]

    profile = build_consensus(ref_seqs, ignore_chars=ignore_chars)
    informative_cols = [i for i, aa in enumerate(profile) if aa != "-" and aa not in ignore_chars]

    if not informative_cols:
        raise SystemExit("No informative profile columns found. Check references/ignore chars.")

    kept: list[Record] = []
    report_lines: list[str] = ["id\tscore\tmatches\ttotal\tkept"]

    ref_id_set = set(ref_ids)
    for header, seq in records:
        seq_id = fasta_id(header)
        score, matches, total = similarity_to_profile(
            seq=seq,
            profile=profile,
            informative_cols=informative_cols,
            ignore_chars=ignore_chars,
            penalize_gaps=args.penalize_gaps,
        )

        keep = score >= args.threshold
        if not args.no_keep_reference and seq_id in ref_id_set:
            keep = True

        if keep:
            kept.append((header, seq))

        report_lines.append(f"{seq_id}\t{score:.4f}\t{matches}\t{total}\t{int(keep)}")

    write_fasta(Path(args.output), kept)

    if args.report:
        Path(args.report).write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Input records: {len(records)}")
    print(f"Kept records: {len(kept)}")
    print(f"Threshold: {args.threshold}")
    print(f"Informative columns: {len(informative_cols)}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

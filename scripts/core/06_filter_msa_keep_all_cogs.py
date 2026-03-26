#!/usr/bin/env python3
"""Filter only structure sequences in MSA while keeping all COG sequences."""

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


def build_consensus(seqs: list[str], ignore_chars: set[str]) -> str:
    length = len(seqs[0])
    out: list[str] = []
    for i in range(length):
        counts = Counter()
        for seq in seqs:
            aa = seq[i]
            if aa == "-" or aa in ignore_chars:
                continue
            counts[aa] += 1
        if not counts:
            out.append("-")
        else:
            out.append(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0])
    return "".join(out)


def score_vs_profile(
    seq: str,
    profile: str,
    informative_cols: list[int],
    ignore_chars: set[str],
    penalize_gaps: bool,
) -> tuple[float, int, int]:
    matches = 0
    total = 0
    for idx in informative_cols:
        aa = seq[idx]
        ref = profile[idx]
        if aa in ignore_chars:
            continue
        if aa == "-" and not penalize_gaps:
            continue
        total += 1
        if aa == ref:
            matches += 1
    if total == 0:
        return 0.0, matches, total
    return matches / total, matches, total


def load_id_set_from_fasta(path: Path) -> set[str]:
    ids: set[str] = set()
    for header, _seq in iter_fasta(path):
        ids.add(fasta_id(header))
    return ids


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Keep all COG sequences and filter only structure records in aligned FASTA"
    )
    parser.add_argument("input_msa", help="Aligned combined FASTA")
    parser.add_argument(
        "--cog-fasta",
        required=True,
        help="FASTA with COG records that must always be kept",
    )
    parser.add_argument("-o", "--output", required=True, help="Output filtered MSA")
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.22,
        help="Similarity threshold for structure records (default: 0.22)",
    )
    parser.add_argument(
        "--target-total",
        type=int,
        default=500,
        help="Target total number of sequences after filtering (default: 500)",
    )
    parser.add_argument(
        "--ignore-chars",
        default="X?",
        help="Characters ignored in scoring (default: X?)",
    )
    parser.add_argument(
        "--penalize-gaps",
        action="store_true",
        help="Count gaps as mismatches during structure scoring",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Optional TSV report",
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

    cog_ids = load_id_set_from_fasta(Path(args.cog_fasta))
    alignment_ids = {fasta_id(header) for header, _seq in records}
    missing_cog_ids = sorted(cog_ids - alignment_ids)
    if missing_cog_ids:
        preview = ", ".join(missing_cog_ids[:10])
        extra = "" if len(missing_cog_ids) <= 10 else f" ... (+{len(missing_cog_ids) - 10} more)"
        raise SystemExit(
            "Some COG records from --cog-fasta are missing in the alignment: "
            f"{preview}{extra}"
        )

    ignore_chars = set(args.ignore_chars)

    cog_records: list[Record] = []
    structure_records: list[Record] = []
    for header, seq in records:
        if fasta_id(header) in cog_ids:
            cog_records.append((header, seq))
        else:
            structure_records.append((header, seq))

    if not cog_records:
        raise SystemExit("No COG records from --cog-fasta were found in the alignment.")

    profile = build_consensus([seq for _, seq in cog_records], ignore_chars=ignore_chars)
    informative_cols = [i for i, aa in enumerate(profile) if aa != "-" and aa not in ignore_chars]
    if not informative_cols:
        raise SystemExit("No informative columns in COG profile.")

    scored_structures: list[tuple[float, int, int, Record]] = []
    report_lines: list[str] = ["id\tgroup\tscore\tmatches\ttotal\tkept\treason"]

    for header, seq in cog_records:
        rec_id = fasta_id(header)
        report_lines.append(f"{rec_id}\tcog\t1.0000\tNA\tNA\t1\talways_keep_cog")

    for rec in structure_records:
        header, seq = rec
        rec_id = fasta_id(header)
        score, matches, total = score_vs_profile(
            seq=seq,
            profile=profile,
            informative_cols=informative_cols,
            ignore_chars=ignore_chars,
            penalize_gaps=args.penalize_gaps,
        )
        scored_structures.append((score, matches, total, rec))
        report_lines.append(f"{rec_id}\tstructure\t{score:.4f}\t{matches}\t{total}\t0\tpending")

    scored_structures.sort(key=lambda x: x[0], reverse=True)

    max_structures = max(0, args.target_total - len(cog_records))
    kept_structures: list[Record] = []
    kept_structure_ids: set[str] = set()

    for i, (score, _matches, _total, rec) in enumerate(scored_structures):
        rec_id = fasta_id(rec[0])
        keep_by_rank = i < max_structures
        keep_by_threshold = score >= args.threshold
        if keep_by_rank and keep_by_threshold:
            kept_structures.append(rec)
            kept_structure_ids.add(rec_id)

    # Update report with final keep/drop decisions for structures.
    final_report: list[str] = [report_lines[0]]
    for line in report_lines[1 : 1 + len(cog_records)]:
        final_report.append(line)
    for score, matches, total, rec in scored_structures:
        rec_id = fasta_id(rec[0])
        kept = rec_id in kept_structure_ids
        reason = "kept" if kept else "filtered_structure"
        final_report.append(
            f"{rec_id}\tstructure\t{score:.4f}\t{matches}\t{total}\t{int(kept)}\t{reason}"
        )

    out_records = cog_records + kept_structures
    write_fasta(Path(args.output), out_records)

    if args.report:
        Path(args.report).write_text("\n".join(final_report) + "\n", encoding="utf-8")

    print(f"COG records kept (forced): {len(cog_records)}")
    print(f"Structure records input: {len(structure_records)}")
    print(f"Structure records kept: {len(kept_structures)}")
    print(f"Target total: {args.target_total}")
    print(f"Total output records: {len(out_records)}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

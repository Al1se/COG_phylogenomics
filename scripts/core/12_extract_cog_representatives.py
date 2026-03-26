#!/usr/bin/env python3
"""Extract one representative sequence per genome for a selected COG."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from cog_extractor_common import iter_fasta_records, parse_gene_and_organism_from_header, wrap_sequence


def load_short_list(path: Path) -> dict[str, str]:
    genome_to_name: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            genome_id, _, name = line.partition("@")
            genome_to_name[genome_id.strip()] = name.strip() if name else genome_id.strip()
    return genome_to_name


def collect_cog_genes_by_genome(
    cog_csv_path: Path,
    cog_id: str,
    allowed_genomes: set[str],
) -> tuple[dict[str, str], dict[str, list[str]]]:
    gene_to_genome: dict[str, str] = {}
    genome_to_genes: dict[str, list[str]] = {}

    with cog_csv_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if len(row) < 8 or row[7] != cog_id:
                continue
            genome_id = row[1].strip()
            if genome_id not in allowed_genomes:
                continue
            gene_id = row[0].strip()
            gene_to_genome[gene_id] = genome_id
            genome_to_genes.setdefault(genome_id, []).append(gene_id)

    return gene_to_genome, genome_to_genes


def normalize_sequence(seq: str) -> str:
    seq = seq.replace("-", "").replace(".", "").upper()
    return "".join(ch for ch in seq if "A" <= ch <= "Z")


def motif_pos(seq: str, motif: str) -> int:
    idx = seq.find(motif)
    return idx + 1 if idx >= 0 else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract one representative sequence per genome for one COG")
    parser.add_argument("cog_id", help="COG identifier, e.g. COG0086")
    parser.add_argument("--short-list", default="275list.txt", help="Genome short list")
    parser.add_argument("--cog-csv", default="cog-24.cog.csv", help="COG CSV table")
    parser.add_argument("--proteins-fasta", default="2296Genomes.prot.fasta", help="Protein FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA")
    parser.add_argument("--motif", default="NADFDGD", help="Preferred motif for representative ranking")
    parser.add_argument("--report", default=None, help="Optional TSV report")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    genome_to_name = load_short_list(Path(args.short_list))
    gene_to_genome, genome_to_genes = collect_cog_genes_by_genome(
        Path(args.cog_csv), args.cog_id, set(genome_to_name)
    )
    if not gene_to_genome:
        raise SystemExit("No matching genes found for the requested COG and genome list.")

    candidates: dict[str, list[tuple[str, str, str, str]]] = {}
    for header, raw_seq in iter_fasta_records(Path(args.proteins_fasta)):
        gene_id, organism = parse_gene_and_organism_from_header(header)
        genome_id = gene_to_genome.get(gene_id)
        if genome_id is None:
            continue
        seq = normalize_sequence(raw_seq)
        name = genome_to_name.get(genome_id, organism or genome_id)
        candidates.setdefault(genome_id, []).append((gene_id, genome_id, name, seq))

    out_lines: list[str] = []
    report_rows: list[dict[str, str]] = []

    for genome_id in sorted(genome_to_name):
        group = candidates.get(genome_id, [])
        if not group:
            report_rows.append(
                {
                    "genome_id": genome_id,
                    "organism": genome_to_name[genome_id],
                    "chosen_id": "",
                    "length": "",
                    "motif_pos": "",
                    "candidate_count": "0",
                    "status": "missing",
                    "candidate_ids": "",
                }
            )
            continue

        scored = []
        for gene_id, _genome_id, name, seq in group:
            pos = motif_pos(seq, args.motif)
            score = (
                1 if pos else 0,
                len(seq),
                -pos if pos else 0,
                gene_id,
            )
            scored.append((score, gene_id, name, seq, pos))
        scored.sort(reverse=True)
        _, gene_id, name, seq, pos = scored[0]
        header = f"{gene_id}|{genome_id}|{name}"
        out_lines.append(f">{header}\n{wrap_sequence(seq)}")
        report_rows.append(
            {
                "genome_id": genome_id,
                "organism": name,
                "chosen_id": gene_id,
                "length": str(len(seq)),
                "motif_pos": str(pos) if pos else "",
                "candidate_count": str(len(group)),
                "status": "chosen",
                "candidate_ids": ";".join(item[1] for item in scored),
            }
        )

    Path(args.output).write_text("\n".join(out_lines) + ("\n" if out_lines else ""), encoding="utf-8")

    if args.report:
        with Path(args.report).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "genome_id",
                    "organism",
                    "chosen_id",
                    "length",
                    "motif_pos",
                    "candidate_count",
                    "status",
                    "candidate_ids",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(report_rows)

    chosen = sum(row["status"] == "chosen" for row in report_rows)
    missing = sum(row["status"] == "missing" for row in report_rows)
    print(f"Requested genomes: {len(genome_to_name)}")
    print(f"Chosen representatives: {chosen}")
    print(f"Missing genomes: {missing}")
    print(f"Output: {args.output}")
    if args.report:
        print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

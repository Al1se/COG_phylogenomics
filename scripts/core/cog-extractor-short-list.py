#!/usr/bin/env python3
"""Extract proteins for one COG restricted to a short genome list."""

from __future__ import annotations

import argparse
from pathlib import Path

from cog_extractor_common import extract_cog_sequences, load_short_genome_ids


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract proteins for one COG using a restricted genome list (e.g. 275list.txt)"
    )
    parser.add_argument("cog_id", help="COG identifier, for example COG0085")
    parser.add_argument(
        "--short-list",
        default="275list.txt",
        help="Path to short genome list (<genome_id>@<name> per line)",
    )
    parser.add_argument("--cog-csv", default="cog-24.cog.csv", help="Path to cog CSV table")
    parser.add_argument(
        "--proteins-fasta",
        default="2296Genomes.prot.fasta",
        help="Path to full protein FASTA database",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output FASTA path (default: 275_<COG_ID>.fasta)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output file if it already exists",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    output_path = Path(args.output) if args.output else Path(f"275_{args.cog_id}.fasta")
    short_genomes = load_short_genome_ids(Path(args.short_list))

    stats = extract_cog_sequences(
        cog_id=args.cog_id,
        cog_csv_path=Path(args.cog_csv),
        proteins_fasta_path=Path(args.proteins_fasta),
        output_fasta_path=output_path,
        allowed_genome_ids=short_genomes,
        overwrite=args.overwrite,
    )

    print(f"COG: {stats.cog_id}")
    print(f"Short-list genomes: {len(short_genomes)}")
    print(f"Matched gene IDs in COG table: {stats.matched_genes}")
    print(f"Written FASTA records: {stats.written_records}")
    print(f"Output: {output_path}")


if __name__ == "__main__":
    main()

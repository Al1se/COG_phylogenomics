#!/usr/bin/env python3
"""Shared helpers for extracting COG protein sequences."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional


@dataclass
class ExtractStats:
    cog_id: str
    matched_genes: int
    written_records: int


def iter_fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    """Yield FASTA records as (header_without_gt, sequence)."""
    header: Optional[str] = None
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


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def load_short_genome_ids(path: Path) -> set[str]:
    """Read 275list-style file with '<genome_id>@<name>' per line."""
    genome_ids: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            genome_ids.add(line.split("@", 1)[0].strip())
    return genome_ids


def collect_gene_ids(
    cog_csv_path: Path,
    cog_id: str,
    allowed_genome_ids: Optional[set[str]] = None,
) -> set[str]:
    """Collect gene IDs from cog-24.cog.csv for the selected COG."""
    gene_ids: set[str] = set()

    with cog_csv_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if len(row) < 8:
                continue
            if row[7] != cog_id:
                continue
            if allowed_genome_ids is not None and row[1] not in allowed_genome_ids:
                continue
            gene_ids.add(row[0])

    return gene_ids


def parse_gene_and_organism_from_header(header: str) -> tuple[str, str]:
    """
    Parse source FASTA header from 2296Genomes.prot.fasta.

    Expected format example:
    gi|Unk|ref|AAR38856.1 NEQ001|NEQ001|...|Nanoarchaeum equitans Kin4-M|...
    """
    parts = header.split("|")

    gene_id = ""
    organism = ""

    if len(parts) > 3:
        gene_field = parts[3].strip()
        gene_tokens = gene_field.split()
        if len(gene_tokens) >= 2:
            gene_id = gene_tokens[1]
        elif gene_tokens:
            gene_id = gene_tokens[0]

    if len(parts) > 6:
        organism = parts[6].strip()

    if not gene_id:
        gene_id = header.split()[0]

    return gene_id, organism


def extract_cog_sequences(
    *,
    cog_id: str,
    cog_csv_path: Path,
    proteins_fasta_path: Path,
    output_fasta_path: Path,
    allowed_genome_ids: Optional[set[str]] = None,
    overwrite: bool = False,
) -> ExtractStats:
    if output_fasta_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output file already exists: {output_fasta_path}. Use --overwrite to replace it."
        )

    gene_ids = collect_gene_ids(
        cog_csv_path=cog_csv_path,
        cog_id=cog_id,
        allowed_genome_ids=allowed_genome_ids,
    )

    written_records = 0
    with output_fasta_path.open("w", encoding="utf-8") as out_handle:
        for header, sequence in iter_fasta_records(proteins_fasta_path):
            gene_id, organism = parse_gene_and_organism_from_header(header)
            if gene_id not in gene_ids:
                continue

            out_header = f"{gene_id}|{organism}" if organism else gene_id
            out_handle.write(f">{out_header}\n")
            out_handle.write(wrap_sequence(sequence) + "\n")
            written_records += 1

    return ExtractStats(cog_id=cog_id, matched_genes=len(gene_ids), written_records=written_records)

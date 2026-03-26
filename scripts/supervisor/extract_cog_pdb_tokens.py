#!/usr/bin/env python3
"""Extract PDB_chain tokens for a requested COG from supervisor mapping files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def normalize_token(token: str) -> str:
    token = token.strip()
    if not token or "_" not in token:
        return ""
    pdb_id, chain_id = token.split("_", 1)
    return f"{pdb_id.upper()}_{chain_id.upper()}"


def read_mapping_tokens(mapping_path: Path, cog_id: str) -> list[str]:
    tokens: list[str] = []
    with mapping_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            if parts[0] != cog_id:
                continue
            token = normalize_token(parts[1])
            if token:
                tokens.append(token)
    return tokens


def read_ligand_filtered_tokens(
    ligands_path: Path,
    cog_id: str,
    requested_ligands: list[str],
) -> list[str]:
    with ligands_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header: list[str] | None = None
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#COG"):
                header = [item.lstrip("#").strip() for item in row]
                continue
            if header is None or row[0] != cog_id:
                continue

            ligand_to_index = {name: idx for idx, name in enumerate(header)}
            seen: set[str] = set()
            tokens: list[str] = []
            for ligand in requested_ligands:
                idx = ligand_to_index.get(ligand)
                if idx is None or idx >= len(row):
                    continue
                cell = row[idx].strip()
                if not cell:
                    continue
                for raw_token in cell.split(","):
                    token = normalize_token(raw_token)
                    if token and token not in seen:
                        seen.add(token)
                        tokens.append(token)
            return tokens
    return []


def dedup_keep_order(tokens: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for token in tokens:
        if token in seen:
            continue
        seen.add(token)
        ordered.append(token)
    return ordered


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve a per-COG PDB_chain list from supervisor mapping files. "
            "If a ligand summary table is supplied, ligand-matching chains are "
            "preferred; otherwise all chains from COG_to_PDB_chains.txt are used."
        )
    )
    parser.add_argument("cog_id", help="COG ID, for example COG0086")
    parser.add_argument(
        "-m",
        "--mapping",
        required=True,
        help="Path to COG_to_PDB_chains.txt",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output text file with one PDB_chain token per line",
    )
    parser.add_argument(
        "--ligands-table",
        default=None,
        help="Optional supervisor *.ligands.txt file to prefer ligand-positive chains",
    )
    parser.add_argument(
        "--ligands",
        default="ZN,MG",
        help="Comma-separated ligand IDs used with --ligands-table (default: ZN,MG)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    cog_id = args.cog_id.strip().upper()
    mapping_path = Path(args.mapping)
    output_path = Path(args.output)

    requested_ligands = [item.strip().upper() for item in args.ligands.split(",") if item.strip()]

    mapping_tokens = dedup_keep_order(read_mapping_tokens(mapping_path, cog_id))
    if not mapping_tokens:
        output_path.write_text("", encoding="utf-8")
        print(f"No PDB_chain tokens found for {cog_id} in {mapping_path}")
        return

    chosen_tokens = mapping_tokens
    if args.ligands_table:
        ligand_tokens = dedup_keep_order(
            read_ligand_filtered_tokens(Path(args.ligands_table), cog_id, requested_ligands)
        )
        if ligand_tokens:
            chosen_tokens = ligand_tokens

    output_path.write_text("".join(f"{token}\n" for token in chosen_tokens), encoding="utf-8")
    print(
        f"Wrote {len(chosen_tokens)} tokens for {cog_id} "
        f"to {output_path}"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Map structure-chain tokens to an alignment by chain sequence."""

from __future__ import annotations

import argparse
import csv
import difflib
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def listify(value):
    if isinstance(value, list):
        return value
    return [value]


def iter_fasta(path: Path):
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


def read_tokens(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def extract_chain_sequence(mmcif_path: Path, token: str) -> str:
    entry_id, chain = token.split("_", 1)
    d = MMCIF2Dict(str(mmcif_path))
    entity_ids = listify(d.get("_entity_poly.entity_id", []))
    strand_ids = listify(d.get("_entity_poly.pdbx_strand_id", []))
    seqs = listify(d.get("_entity_poly.pdbx_seq_one_letter_code_can", []))

    for _entity_id, strands, raw_seq in zip(entity_ids, strand_ids, seqs):
        chains = [part.strip() for part in str(strands).split(",") if part.strip()]
        if chain in chains:
            return "".join(str(raw_seq).split())
    raise KeyError(f"Chain sequence not found for token {token} in {mmcif_path}")


def build_alignment_dict(msa_fasta: Path) -> dict[str, str]:
    return {header.split()[0]: seq.replace("-", "") for header, seq in iter_fasta(msa_fasta)}


def diff_count(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for aa, bb in zip(a, b) if aa != bb)


def main() -> None:
    parser = argparse.ArgumentParser(description="Map structure-chain tokens to MSA entries by sequence")
    parser.add_argument("token_list", help="Text file with one PDB_chain token per line")
    parser.add_argument("--mmcif-dir", required=True, help="Directory with mmCIF files")
    parser.add_argument("--msa-fasta", required=True, help="Target aligned FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument(
        "--near-threshold",
        type=float,
        default=0.99,
        help="Minimum similarity ratio for a near match when exact match is absent",
    )
    args = parser.parse_args()

    msa = build_alignment_dict(Path(args.msa_fasta))
    rows: list[dict[str, str]] = []

    for token in read_tokens(Path(args.token_list)):
        entry_id, _chain = token.split("_", 1)
        seq = extract_chain_sequence(Path(args.mmcif_dir) / f"{entry_id}.cif", token)
        exact_hits = [msa_id for msa_id, msa_seq in msa.items() if msa_seq == seq]

        row = {
            "token": token,
            "structure_length": str(len(seq)),
            "match_type": "none",
            "representative_msa_id": "",
            "all_msa_ids": "",
            "ratio": "",
            "diff_count": "",
        }

        if exact_hits:
            row["match_type"] = "exact"
            row["representative_msa_id"] = exact_hits[0]
            row["all_msa_ids"] = ";".join(exact_hits)
            row["ratio"] = "1.000000"
            row["diff_count"] = "0"
            rows.append(row)
            continue

        best_ratio = -1.0
        best_ids: list[str] = []
        best_seq = ""
        for msa_id, msa_seq in msa.items():
            if abs(len(msa_seq) - len(seq)) > 10:
                continue
            ratio = difflib.SequenceMatcher(None, seq, msa_seq).ratio()
            if ratio > best_ratio + 1e-12:
                best_ratio = ratio
                best_ids = [msa_id]
                best_seq = msa_seq
            elif abs(ratio - best_ratio) <= 1e-12:
                best_ids.append(msa_id)

        if best_ratio >= args.near_threshold and best_ids:
            row["match_type"] = "near"
            row["representative_msa_id"] = best_ids[0]
            row["all_msa_ids"] = ";".join(best_ids)
            row["ratio"] = f"{best_ratio:.6f}"
            row["diff_count"] = str(diff_count(seq, best_seq))

        rows.append(row)

    with Path(args.output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "token",
                "structure_length",
                "match_type",
                "representative_msa_id",
                "all_msa_ids",
                "ratio",
                "diff_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    exact_count = sum(row["match_type"] == "exact" for row in rows)
    near_count = sum(row["match_type"] == "near" for row in rows)
    print(f"Exact matches: {exact_count}")
    print(f"Near matches: {near_count}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

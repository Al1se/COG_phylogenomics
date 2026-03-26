#!/usr/bin/env python3
"""Summarize structure-mapped metal-binding sites in a reordered alignment."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


PREFERRED_ATOMS = {
    "CYS": {"SG"},
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "HIS": {"ND1", "NE2"},
    "MET": {"SD"},
    "ASN": {"OD1"},
    "GLN": {"OE1"},
    "SER": {"OG"},
    "THR": {"OG1"},
    "TYR": {"OH"},
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
}


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


def read_query_ids(path: Path | None) -> list[str]:
    if path is None or not path.exists():
        return []
    return [header.split()[0] for header, _seq in iter_fasta(path)]


def read_alignment(path: Path) -> dict[str, str]:
    return {header.split()[0]: seq for header, seq in iter_fasta(path)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize mapped structure sites in reordered MSA")
    parser.add_argument("mapped_contacts_tsv", help="TSV from map_contacts_via_sequence_matches.py")
    parser.add_argument("--sites-tsv", default=None, help="Optional TSV from analyze_metal_contacts.py")
    parser.add_argument("--msa-fasta", required=True, help="Reordered aligned FASTA")
    parser.add_argument(
        "--query-fasta",
        default=None,
        help="Optional FASTA with query/non-structure sequences to summarize site patterns",
    )
    parser.add_argument(
        "--title",
        default="STRUCTURE SITE COORDINATES IN REORDERED MSA",
        help="Title for the text report",
    )
    parser.add_argument("--output-txt", required=True, help="Human-readable output report")
    parser.add_argument("--output-tsv", default=None, help="Optional TSV with one row per summarized site")
    return parser.parse_args()


def choose_column(rows: list[dict[str, str]]) -> str:
    if not rows:
        return ""
    resname = rows[0]["resname"]
    preferred = PREFERRED_ATOMS.get(resname, set())
    filtered = [row for row in rows if row["atom_name"] in preferred and row["alignment_column"]]
    if preferred:
        working = filtered
    else:
        working = [row for row in rows if row["alignment_column"]]
    if not working:
        return ""
    counts = Counter(row["alignment_column"] for row in working)
    best_col, _count = sorted(counts.items(), key=lambda x: (-x[1], int(x[0])))[0]
    return best_col


def main() -> None:
    args = parse_args()

    msa = read_alignment(Path(args.msa_fasta))
    query_ids = [qid for qid in read_query_ids(Path(args.query_fasta)) if qid in msa]

    site_meta: dict[tuple[str, str], dict[str, str]] = {}
    if args.sites_tsv:
        with Path(args.sites_tsv).open("r", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                site_meta[(row["token"], row["metal_site_id"])] = row

    grouped_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    residue_rows: dict[tuple[str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)

    with Path(args.mapped_contacts_tsv).open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row.get("alignment_column"):
                continue
            site_key = (row["token"], row["metal_site_id"])
            grouped_rows[site_key].append(row)
            residue_key = (
                row["token"],
                row["metal_site_id"],
                row["resname"],
                row["auth_seq_id"],
                row["label_seq_id"],
            )
            residue_rows[residue_key].append(row)

    raw_site_summaries: list[dict[str, str]] = []
    for site_key, rows in sorted(grouped_rows.items()):
        token, metal_site_id = site_key
        exemplar = rows[0]
        residue_entries = []
        residue_columns: list[int] = []

        residue_keys = sorted(
            [
                key for key in residue_rows
                if key[0] == token and key[1] == metal_site_id
            ],
            key=lambda key: (
                int(key[4]) if key[4] not in {"", ".", "?"} else 10**9,
                int(key[3]) if key[3] not in {"", ".", "?"} else 10**9,
                key[2],
            ),
        )

        for res_key in residue_keys:
            res_rows = residue_rows[res_key]
            column = choose_column(res_rows)
            if not column:
                continue
            residue_entries.append(
                f"{res_key[2]}{res_key[3]}[seq={res_key[4]}]->col={column}"
            )
            residue_columns.append(int(column))

        if not residue_columns:
            continue

        raw_site_summaries.append(
            {
                "token": token,
                "pdb_id": exemplar["pdb_id"],
                "target_chain": exemplar["target_chain"],
                "metal": exemplar["metal"],
                "metal_site_id": metal_site_id,
                "reference_residues": ";".join(residue_entries),
                "columns": ",".join(str(col) for col in residue_columns),
                "min_column": str(min(residue_columns)),
                "residue_count": str(len(residue_columns)),
                "site_residue_summary": site_meta.get(site_key, {}).get("residues", ""),
            }
        )

    merged: dict[tuple[str, str], dict[str, str]] = {}
    merged_tokens: dict[tuple[str, str], list[str]] = defaultdict(list)
    merged_pdb_chains: dict[tuple[str, str], list[str]] = defaultdict(list)
    merged_raw_residues: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in raw_site_summaries:
        key = (row["metal"], row["columns"])
        if key not in merged:
            merged[key] = dict(row)
        merged_tokens[key].append(row["token"])
        merged_pdb_chains[key].append(f"{row['pdb_id']}_{row['target_chain']}")
        if row["site_residue_summary"]:
            merged_raw_residues[key].append(row["site_residue_summary"])

    site_summaries: list[dict[str, str]] = []
    for key, row in merged.items():
        row["supporting_tokens"] = ";".join(sorted(dict.fromkeys(merged_tokens[key])))
        row["supporting_pdb_chains"] = ";".join(sorted(dict.fromkeys(merged_pdb_chains[key])))
        row["site_residue_summary"] = ";".join(sorted(dict.fromkeys(merged_raw_residues[key])))
        site_summaries.append(row)

    site_summaries.sort(key=lambda row: (int(row["min_column"]), row["columns"], row["metal"]))

    output_lines = [args.title, ""]
    output_lines.append(f"Reordered alignment: {args.msa_fasta}")
    output_lines.append(f"Mapped structure sites: {len(site_summaries)}")
    output_lines.append("")

    for idx, row in enumerate(site_summaries, start=1):
        output_lines.append(f"Site {idx}")
        output_lines.append(f"Representative token: {row['token']}")
        output_lines.append(f"Supporting tokens: {row['supporting_tokens']}")
        output_lines.append(f"Supporting PDB/chains: {row['supporting_pdb_chains']}")
        output_lines.append(f"Metal: {row['metal']}")
        output_lines.append(f"Representative metal site id: {row['metal_site_id']}")
        if row["site_residue_summary"]:
            output_lines.append(f"Site residues (raw): {row['site_residue_summary']}")
        output_lines.append(f"Mapped residues: {row['reference_residues']}")
        output_lines.append(f"Reordered-MSA columns: {row['columns']}")

        if query_ids:
            cols = [int(x) for x in row["columns"].split(",") if x]
            pattern_counts = Counter(
                "".join(msa[qid][col - 1] for col in cols)
                for qid in query_ids
            )
            top_patterns = ", ".join(f"{pat}={count}" for pat, count in pattern_counts.most_common(8))
            output_lines.append(f"Top query patterns: {top_patterns}")
        output_lines.append("")

    Path(args.output_txt).write_text("\n".join(output_lines) + "\n", encoding="utf-8")

    if args.output_tsv:
        with Path(args.output_tsv).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "token",
                    "pdb_id",
                    "target_chain",
                    "metal",
                    "metal_site_id",
                    "supporting_tokens",
                    "supporting_pdb_chains",
                    "reference_residues",
                    "columns",
                    "min_column",
                    "residue_count",
                    "site_residue_summary",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(site_summaries)

    print(f"Sites summarized: {len(site_summaries)}")
    print(f"Output TXT: {args.output_txt}")
    if args.output_tsv:
        print(f"Output TSV: {args.output_tsv}")


if __name__ == "__main__":
    main()

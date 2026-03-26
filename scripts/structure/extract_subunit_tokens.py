#!/usr/bin/env python3
"""Extract beta and beta' chain tokens from downloaded mmCIF files by entity description."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def read_entries(path: Path) -> list[str]:
    entries: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        entries.append(line.upper())
    return entries


def as_list(value):
    if isinstance(value, list):
        return value
    return [value]


def desc_to_kind(desc: str) -> str | None:
    text = desc.lower().strip()
    if "subunit beta'" in text or "subunit beta prime" in text:
        return "betaprime"
    if "subunit beta" in text and "beta'" not in text and "beta prime" not in text:
        return "beta"
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Build token lists for RNAP beta and beta' chains")
    parser.add_argument("entries_txt", help="Text file with one 4-letter PDB entry per line")
    parser.add_argument("--mmcif-dir", required=True, help="Directory with downloaded .cif files")
    parser.add_argument("--beta-out", required=True, help="Output token list for beta chains")
    parser.add_argument("--betaprime-out", required=True, help="Output token list for beta' chains")
    parser.add_argument("--report", required=True, help="Output TSV summary")
    args = parser.parse_args()

    beta_tokens: list[str] = []
    betaprime_tokens: list[str] = []
    report_lines = ["pdb_id\tkind\tentity_id\tdescription\tchains"]

    for entry in read_entries(Path(args.entries_txt)):
        cif_path = Path(args.mmcif_dir) / f"{entry}.cif"
        if not cif_path.exists():
            continue
        d = MMCIF2Dict(str(cif_path))
        entity_ids = as_list(d.get("_entity_poly.entity_id", []))
        strand_ids = as_list(d.get("_entity_poly.pdbx_strand_id", []))
        descriptions = as_list(d.get("_entity.pdbx_description", []))
        desc_map = {str(i + 1): desc for i, desc in enumerate(descriptions)}

        for ent_id, strands in zip(entity_ids, strand_ids):
            desc = desc_map.get(ent_id, "")
            kind = desc_to_kind(desc)
            if kind is None:
                continue
            chains = [s.strip() for s in strands.split(",") if s.strip()]
            report_lines.append(f"{entry}\t{kind}\t{ent_id}\t{desc}\t{','.join(chains)}")
            for chain in chains:
                token = f"{entry}_{chain}"
                if kind == "beta":
                    beta_tokens.append(token)
                else:
                    betaprime_tokens.append(token)

    Path(args.beta_out).write_text("\n".join(beta_tokens) + ("\n" if beta_tokens else ""), encoding="utf-8")
    Path(args.betaprime_out).write_text(
        "\n".join(betaprime_tokens) + ("\n" if betaprime_tokens else ""), encoding="utf-8"
    )
    Path(args.report).write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"beta tokens: {len(beta_tokens)}")
    print(f"beta' tokens: {len(betaprime_tokens)}")
    print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

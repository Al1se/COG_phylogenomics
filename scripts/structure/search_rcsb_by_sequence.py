#!/usr/bin/env python3
"""Search RCSB PDB by protein sequence similarity and emit PDB_chain tokens."""

from __future__ import annotations

import argparse
import csv
import json
import urllib.request
from pathlib import Path
from typing import Any


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}"


def iter_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def post_json(url: str, payload: dict[str, Any], timeout: float) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def get_json(url: str, timeout: float) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def walk_match_context(node: Any) -> list[dict[str, Any]]:
    contexts: list[dict[str, Any]] = []
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "match_context" and isinstance(value, list):
                contexts.extend([item for item in value if isinstance(item, dict)])
            else:
                contexts.extend(walk_match_context(value))
    elif isinstance(node, list):
        for item in node:
            contexts.extend(walk_match_context(item))
    return contexts


def parse_float(value: Any) -> float | None:
    if value in (None, "", ".", "?"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def build_payload(sequence: str, evalue_cutoff: float, identity_cutoff: float, max_rows: int) -> dict[str, Any]:
    return {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "target": "pdb_protein_sequence",
                "value": sequence,
                "identity_cutoff": identity_cutoff,
                "evalue_cutoff": evalue_cutoff,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "scoring_strategy": "sequence",
            "results_content_type": ["experimental"],
            "results_verbosity": "verbose",
            "paginate": {"start": 0, "rows": max_rows},
        },
    }


def resolve_auth_chains(entry_id: str, entity_id: str, timeout: float) -> list[str]:
    data = get_json(ENTITY_URL.format(entry_id=entry_id, entity_id=entity_id), timeout=timeout)
    ids = data.get("rcsb_polymer_entity_container_identifiers", {})
    auth_asym_ids = ids.get("auth_asym_ids") or []
    if isinstance(auth_asym_ids, str):
        auth_asym_ids = [auth_asym_ids]
    cleaned = []
    for chain_id in auth_asym_ids:
        chain = str(chain_id).strip()
        if chain:
            cleaned.append(chain.upper())
    return cleaned


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Search RCSB PDB for experimental protein structures similar to one or more "
            "representative query sequences and emit PDB_chain tokens."
        )
    )
    parser.add_argument("query_fasta", help="Representative query FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output token list")
    parser.add_argument("--report", required=True, help="Output TSV report")
    parser.add_argument(
        "--evalue-cutoff",
        type=float,
        default=1e-10,
        help="Maximum E-value for the RCSB sequence search (default: 1e-10)",
    )
    parser.add_argument(
        "--identity-cutoff",
        type=float,
        default=0.30,
        help="Minimum sequence identity fraction for the RCSB search (default: 0.30)",
    )
    parser.add_argument(
        "--max-hits-per-query",
        type=int,
        default=100,
        help="Maximum polymer-entity hits retrieved per query (default: 100)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="HTTP timeout in seconds (default: 30)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    query_records = iter_fasta_records(Path(args.query_fasta))
    if not query_records:
        raise SystemExit("No query sequences found.")

    report_rows: list[dict[str, str]] = []
    token_order: list[str] = []
    token_seen: set[str] = set()

    for query_header, sequence in query_records:
        payload = build_payload(
            sequence=sequence,
            evalue_cutoff=args.evalue_cutoff,
            identity_cutoff=args.identity_cutoff,
            max_rows=args.max_hits_per_query,
        )
        response = post_json(SEARCH_URL, payload, timeout=args.timeout)
        results = response.get("result_set", [])

        for result in results:
            identifier = str(result.get("identifier", "")).strip()
            if not identifier or "_" not in identifier:
                continue
            entry_id, entity_id = identifier.split("_", 1)
            match_contexts = walk_match_context(result)
            sequence_identity = None
            evalue = None
            bitscore = None
            for context in match_contexts:
                if sequence_identity is None:
                    sequence_identity = parse_float(context.get("sequence_identity"))
                if evalue is None:
                    evalue = parse_float(context.get("evalue"))
                if bitscore is None:
                    bitscore = parse_float(context.get("bitscore"))

            if sequence_identity is not None and sequence_identity < args.identity_cutoff:
                continue
            if evalue is not None and evalue > args.evalue_cutoff:
                continue

            chains = resolve_auth_chains(entry_id.upper(), entity_id, timeout=args.timeout)
            if not chains:
                continue

            chain_tokens = []
            for chain in chains:
                token = f"{entry_id.upper()}_{chain}"
                chain_tokens.append(token)
                if token not in token_seen:
                    token_seen.add(token)
                    token_order.append(token)

            report_rows.append(
                {
                    "query_header": query_header,
                    "query_length": str(len(sequence)),
                    "polymer_entity": identifier,
                    "sequence_identity": "" if sequence_identity is None else f"{sequence_identity:.4f}",
                    "evalue": "" if evalue is None else f"{evalue:.3g}",
                    "bitscore": "" if bitscore is None else f"{bitscore:.3f}",
                    "chain_tokens": ",".join(chain_tokens),
                }
            )

    Path(args.output).write_text("".join(f"{token}\n" for token in token_order), encoding="utf-8")
    with Path(args.report).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "query_header",
                "query_length",
                "polymer_entity",
                "sequence_identity",
                "evalue",
                "bitscore",
                "chain_tokens",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(report_rows)

    print(f"Queries: {len(query_records)}")
    print(f"Matched polymer entities: {len(report_rows)}")
    print(f"Unique PDB_chain tokens: {len(token_order)}")
    print(f"Output token list: {args.output}")
    print(f"Report: {args.report}")


if __name__ == "__main__":
    main()

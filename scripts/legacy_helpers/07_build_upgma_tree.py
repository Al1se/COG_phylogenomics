#!/usr/bin/env python3
"""Build an approximate UPGMA tree from aligned FASTA and save Newick."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterator

import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist


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


def fasta_id(header: str) -> str:
    return header.split()[0]


def quote_newick_label(label: str) -> str:
    safe = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-|")
    if label and all(ch in safe for ch in label):
        return label
    return "'" + label.replace("'", "''") + "'"


def to_newick(node, parent_dist: float, leaf_names: list[str]) -> str:
    if node.is_leaf():
        branch = max(parent_dist - node.dist, 0.0)
        return f"{quote_newick_label(leaf_names[node.id])}:{branch:.8f}"
    left = to_newick(node.left, node.dist, leaf_names)
    right = to_newick(node.right, node.dist, leaf_names)
    if parent_dist < 0:
        return f"({left},{right})"
    branch = max(parent_dist - node.dist, 0.0)
    return f"({left},{right}):{branch:.8f}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build UPGMA tree from aligned FASTA")
    parser.add_argument("input_msa", help="Aligned FASTA input")
    parser.add_argument("-o", "--output", required=True, help="Output Newick path")
    parser.add_argument(
        "--max-seqs",
        type=int,
        default=500,
        help="Maximum number of sequences to use (default: 500, keep input order)",
    )
    parser.add_argument(
        "--allow-truncate",
        action="store_true",
        help="Allow truncation when input exceeds --max-seqs (default: fail instead of truncating)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    records = list(iter_fasta(Path(args.input_msa)))
    if not records:
        raise SystemExit("Input MSA is empty.")

    if args.max_seqs > 0 and len(records) > args.max_seqs:
        if not args.allow_truncate:
            raise SystemExit(
                f"Input has {len(records)} sequences, which exceeds --max-seqs={args.max_seqs}. "
                "Increase --max-seqs, reduce input size earlier in the pipeline, or pass "
                "--allow-truncate for an exploratory run."
            )
        records = records[: args.max_seqs]

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise SystemExit("Input is not aligned: sequence lengths are different.")

    ids = [fasta_id(h) for h, _ in records]
    seq_len = len(records[0][1])
    n = len(records)

    if n == 1:
        newick = f"{quote_newick_label(ids[0])}:0.00000000;\n"
        Path(args.output).write_text(newick, encoding="utf-8")
        print("Input records used: 1")
        print(f"Aligned length: {seq_len}")
        print("Method: trivial single-leaf tree")
        print(f"Output: {args.output}")
        return

    # Encode sequences as uint8 matrix for fast Hamming distance via scipy.pdist.
    mat = np.zeros((n, seq_len), dtype=np.uint8)
    for i, (_h, seq) in enumerate(records):
        mat[i, :] = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)

    dist = pdist(mat, metric="hamming")
    link = linkage(dist, method="average")
    root = to_tree(link, rd=False)
    newick = to_newick(root, -1.0, ids) + ";\n"
    Path(args.output).write_text(newick, encoding="utf-8")

    print(f"Input records used: {n}")
    print(f"Aligned length: {seq_len}")
    print("Method: UPGMA on Hamming distances")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()

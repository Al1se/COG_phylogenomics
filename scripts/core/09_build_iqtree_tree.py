#!/usr/bin/env python3
"""Build a phylogenetic tree from aligned FASTA using IQ-TREE."""

from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
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


def prepare_iqtree_input(input_msa: Path, out_fasta: Path, out_idmap_tsv: Path) -> tuple[int, int]:
    records = list(iter_fasta(input_msa))
    if not records:
        raise SystemExit("Input MSA is empty.")

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise SystemExit("Input FASTA is not an alignment: sequence lengths differ.")
    aln_len = next(iter(lengths))

    normalized: list[Record] = []
    seen_ids: set[str] = set()
    idmap_lines = ["id\toriginal_header"]

    for header, seq in records:
        rec_id = fasta_id(header)
        if not rec_id:
            raise SystemExit(f"Failed to derive ID from header: {header!r}")
        if rec_id in seen_ids:
            raise SystemExit(
                "Duplicate FASTA ID after normalization (first token). "
                f"Please fix input headers: {rec_id}"
            )
        seen_ids.add(rec_id)
        normalized.append((rec_id, seq))
        idmap_lines.append(f"{rec_id}\t{header}")

    write_fasta(out_fasta, normalized)
    out_idmap_tsv.write_text("\n".join(idmap_lines) + "\n", encoding="utf-8")
    return len(normalized), aln_len


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build IQ-TREE from aligned FASTA")
    parser.add_argument("input_msa", help="Aligned FASTA input")
    parser.add_argument("-o", "--output", required=True, help="Output Newick path (.nwk)")
    parser.add_argument(
        "--iqtree-bin",
        default="iqtree2",
        help="IQ-TREE executable (default: iqtree2)",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="IQ-TREE output prefix (default: output path without suffix)",
    )
    parser.add_argument(
        "--model",
        default="MFP",
        help="Substitution model, e.g. LG+G4 or MFP (default: MFP)",
    )
    parser.add_argument(
        "--threads",
        default="AUTO",
        help="Number of threads for IQ-TREE (default: AUTO)",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="Ultrafast bootstrap replicates, 0 to disable (default: 1000)",
    )
    parser.add_argument(
        "--alrt",
        type=int,
        default=1000,
        help="SH-aLRT replicates, 0 to disable (default: 1000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for IQ-TREE (default: 42)",
    )
    parser.add_argument(
        "--redo",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Pass -redo to IQ-TREE (default: on)",
    )
    parser.add_argument(
        "--bnni",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use --bnni with UFBoot (default: on)",
    )
    parser.add_argument(
        "--keep-ident",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Pass --keep-ident to keep identical sequences in final tree (default: on)",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Optional log path; default: <prefix>.cmd.log",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Pass -fast to IQ-TREE for a quicker approximate search",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    input_msa = Path(args.input_msa)
    output_nwk = Path(args.output)
    output_nwk.parent.mkdir(parents=True, exist_ok=True)

    prefix = Path(args.prefix) if args.prefix else output_nwk.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)

    iqtree_input = Path(str(prefix) + ".ids.fasta")
    idmap_tsv = Path(str(prefix) + ".idmap.tsv")
    n_records, aln_len = prepare_iqtree_input(input_msa, iqtree_input, idmap_tsv)

    cmd = [
        args.iqtree_bin,
        "-s",
        str(iqtree_input),
        "-st",
        "AA",
        "-m",
        args.model,
        "-nt",
        str(args.threads),
        "-pre",
        str(prefix),
    ]

    if args.bootstrap > 0:
        cmd += ["-B", str(args.bootstrap)]
    if args.alrt > 0:
        cmd += ["--alrt", str(args.alrt)]
    if args.bnni and args.bootstrap > 0:
        cmd += ["--bnni"]
    if args.keep_ident:
        cmd += ["--keep-ident"]
    if args.fast:
        cmd += ["-fast"]
    if args.seed is not None:
        cmd += ["-seed", str(args.seed)]
    if args.redo:
        cmd += ["-redo"]

    log_path = Path(args.log) if args.log else Path(str(prefix) + ".cmd.log")
    quoted_cmd = " ".join(shlex.quote(part) for part in cmd)
    print("Running:", quoted_cmd)

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write("COMMAND:\n" + quoted_cmd + "\n\n")
        log_handle.write("STDOUT+STDERR:\n")
        log_handle.flush()
        run = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, text=True)

    if run.returncode != 0:
        raise SystemExit(
            "IQ-TREE failed.\n"
            f"Exit code: {run.returncode}\n"
            f"Check log: {log_path}"
        )

    treefile = Path(str(prefix) + ".treefile")
    if not treefile.exists():
        raise SystemExit(f"IQ-TREE finished but tree file is missing: {treefile}")

    shutil.copyfile(treefile, output_nwk)

    print(f"Input records used: {n_records}")
    print(f"Aligned length: {aln_len}")
    print(f"IQ-TREE input (ID-normalized): {iqtree_input}")
    print(f"ID map TSV: {idmap_tsv}")
    print(f"Treefile: {treefile}")
    print(f"Output Newick: {output_nwk}")
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()

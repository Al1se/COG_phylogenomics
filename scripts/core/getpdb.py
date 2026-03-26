#!/usr/bin/env python3
"""Fetch FASTA sequences for PDB chains listed in a text file."""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from typing import Iterator
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


TOKEN_RE = re.compile(r"^[A-Za-z0-9]{4}_[A-Za-z0-9]+$")


def iter_fasta_records(fasta_text: str) -> Iterator[tuple[str, str]]:
    header = None
    seq_chunks: list[str] = []

    for raw_line in fasta_text.splitlines():
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


def parse_pdb_chain_list(path: Path) -> list[tuple[str, str]]:
    content = path.read_text(encoding="utf-8")
    raw_tokens = [token.strip() for token in re.split(r"[,\n\r\t ]+", content) if token.strip()]

    pairs: list[tuple[str, str]] = []
    invalid_tokens: list[str] = []

    for token in raw_tokens:
        if not TOKEN_RE.match(token):
            invalid_tokens.append(token)
            continue
        pdb_id, chain_id = token.split("_", 1)
        pairs.append((pdb_id.upper(), chain_id.upper()))

    if invalid_tokens:
        for token in invalid_tokens:
            print(f"Warning: skipped malformed token '{token}'", file=sys.stderr)

    return pairs


def fetch_url_text(url: str, timeout: float) -> str:
    with urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_pdb_fasta(pdb_id: str, timeout: float, retries: int) -> str:
    urls = [
        f"https://www.rcsb.org/fasta/entry/{pdb_id}/download",
        f"https://www.rcsb.org/fasta/entry/{pdb_id}",
    ]

    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        for url in urls:
            try:
                return fetch_url_text(url=url, timeout=timeout)
            except (HTTPError, URLError, TimeoutError) as exc:
                last_error = exc
                continue
        time.sleep(min(2.0, 0.4 * attempt))

    if last_error is not None:
        raise RuntimeError(f"Failed to fetch FASTA for {pdb_id}: {last_error}") from last_error
    raise RuntimeError(f"Failed to fetch FASTA for {pdb_id}")


def extract_chains_from_header(header: str) -> set[str]:
    chains: set[str] = set()
    for match in re.finditer(r"\bChains?\s+([^|]+)", header, flags=re.IGNORECASE):
        block = match.group(1)
        block = re.sub(r"\band\b", ",", block, flags=re.IGNORECASE)
        for token in re.split(r"[,/ ]+", block):
            token = token.strip().upper()
            if token:
                chains.add(token)
    return chains


def header_matches_chain(header: str, pdb_id: str, chain_id: str) -> bool:
    upper_header = header.upper()

    if f"{pdb_id}_{chain_id}" in upper_header:
        return True

    if re.search(rf"\b{re.escape(pdb_id)}[ _-]{re.escape(chain_id)}\b", upper_header):
        return True

    chains = extract_chains_from_header(upper_header)
    if chain_id in chains:
        return True

    return False


def fetch_chain_sequence(pdb_id: str, chain_id: str, timeout: float, retries: int) -> str:
    fasta_text = fetch_pdb_fasta(pdb_id=pdb_id, timeout=timeout, retries=retries)
    for header, sequence in iter_fasta_records(fasta_text):
        if header_matches_chain(header=header, pdb_id=pdb_id, chain_id=chain_id):
            return sequence
    return ""


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Download FASTA sequences for PDB_chain IDs")
    parser.add_argument("input_file", help="Text file with tokens like 5TJG_D, 6C6T_J")
    parser.add_argument(
        "-o",
        "--output",
        default="sequences.fasta",
        help="Output FASTA path (default: sequences.fasta)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.3,
        help="Delay between processed records in seconds (default: 0.3)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=15.0,
        help="HTTP timeout in seconds (default: 15)",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of fetch retry rounds (default: 3)",
    )
    parser.add_argument(
        "--failed",
        default=None,
        help="Optional file to save failed PDB_chain tokens",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    pairs = parse_pdb_chain_list(Path(args.input_file))
    if not pairs:
        raise SystemExit("No valid PDB_chain tokens found in input file.")

    successes = 0
    failures: list[str] = []

    with Path(args.output).open("w", encoding="utf-8") as out_handle:
        for idx, (pdb_id, chain_id) in enumerate(pairs, start=1):
            token = f"{pdb_id}_{chain_id}"
            print(f"[{idx}/{len(pairs)}] Fetching {token}...")
            try:
                sequence = fetch_chain_sequence(
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    timeout=args.timeout,
                    retries=args.retries,
                )
            except Exception as exc:  # noqa: BLE001
                print(f"  Error: {exc}", file=sys.stderr)
                failures.append(token)
                time.sleep(args.delay)
                continue

            if not sequence:
                print(f"  Not found: chain {chain_id} in PDB {pdb_id}", file=sys.stderr)
                failures.append(token)
                time.sleep(args.delay)
                continue

            out_handle.write(f">{token}\n{wrap_sequence(sequence)}\n")
            successes += 1
            time.sleep(args.delay)

    if args.failed:
        Path(args.failed).write_text("\n".join(failures) + ("\n" if failures else ""), encoding="utf-8")

    print(f"Done. Success: {successes}, Failed: {len(failures)}")
    print(f"Output: {args.output}")
    if args.failed:
        print(f"Failed list: {args.failed}")


if __name__ == "__main__":
    main()

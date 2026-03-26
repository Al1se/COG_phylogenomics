#!/usr/bin/env python3
"""Download PDB-format files from RCSB for tokens like 5TJG_D or 6UU8_DDD."""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


TOKEN_RE = re.compile(r"^([A-Za-z0-9]{4})_(.+)$")


def read_tokens(path: Path) -> list[str]:
    tokens: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tokens.append(line)
    return tokens


def parse_entry_id(token: str) -> str:
    match = TOKEN_RE.match(token)
    if not match:
        raise ValueError(f"Malformed token: {token}")
    return match.group(1).upper()


def fetch(url: str, timeout: float) -> bytes:
    with urlopen(url, timeout=timeout) as response:
        return response.read()


def main() -> None:
    parser = argparse.ArgumentParser(description="Download RCSB PDB files for token lists")
    parser.add_argument("token_list", help="Text file with one PDB_chain token per line")
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory for downloaded .pdb files",
    )
    parser.add_argument("--timeout", type=float, default=30.0, help="HTTP timeout in seconds")
    parser.add_argument("--retries", type=int, default=3, help="Retry rounds per entry")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tokens = read_tokens(Path(args.token_list))
    if not tokens:
        raise SystemExit("Token list is empty.")

    seen_entries: set[str] = set()
    failures: list[str] = []

    for token in tokens:
        try:
            entry_id = parse_entry_id(token)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            failures.append(token)
            continue

        if entry_id in seen_entries:
            continue
        seen_entries.add(entry_id)

        out_path = output_dir / f"{entry_id}.pdb"
        if out_path.exists():
            print(f"Skip existing: {out_path.name}")
            continue

        url = f"https://files.rcsb.org/download/{entry_id}.pdb"
        last_error: Exception | None = None
        for attempt in range(1, args.retries + 1):
            try:
                payload = fetch(url, timeout=args.timeout)
                out_path.write_bytes(payload)
                print(f"Downloaded: {entry_id}.pdb")
                break
            except (HTTPError, URLError, TimeoutError) as exc:
                last_error = exc
                time.sleep(min(2.0, 0.5 * attempt))
        else:
            failures.append(entry_id)
            print(f"Failed: {entry_id} ({last_error})", file=sys.stderr)

    print(f"Unique entries requested: {len(seen_entries)}")
    print(f"Failures: {len(failures)}")
    if failures:
        print("Failed entries:", ", ".join(failures), file=sys.stderr)


if __name__ == "__main__":
    main()

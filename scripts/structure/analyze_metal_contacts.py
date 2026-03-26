#!/usr/bin/env python3
"""Analyze Zn/Mg contacts around a target chain in RCSB mmCIF files."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


TOKEN_RE = re.compile(r"^([A-Za-z0-9]{4})_(.+)$")
DEFAULT_METALS = ("ZN", "MG")
AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


@dataclass(frozen=True)
class AtomRow:
    group: str
    atom_id: str
    comp_id: str
    auth_asym_id: str
    label_asym_id: str
    auth_seq_id: str
    label_seq_id: str
    alt_id: str
    x: float
    y: float
    z: float
    occupancy: str
    element: str


def listify(value: object) -> list[str]:
    if isinstance(value, list):
        return [str(v) for v in value]
    return [str(value)]


def parse_token(token: str) -> tuple[str, str]:
    match = TOKEN_RE.match(token)
    if not match:
        raise ValueError(f"Malformed token: {token}")
    return match.group(1).upper(), match.group(2)


def read_tokens(path: Path) -> list[str]:
    tokens: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tokens.append(line)
    return tokens


def iter_atom_rows(mmcif_path: Path) -> list[AtomRow]:
    mmcif = MMCIF2Dict(str(mmcif_path))
    group = listify(mmcif["_atom_site.group_PDB"])
    atom_id = listify(mmcif["_atom_site.label_atom_id"])
    comp_id = listify(mmcif["_atom_site.label_comp_id"])
    auth_asym = listify(mmcif["_atom_site.auth_asym_id"])
    label_asym = listify(mmcif["_atom_site.label_asym_id"])
    auth_seq = listify(mmcif["_atom_site.auth_seq_id"])
    label_seq = listify(mmcif["_atom_site.label_seq_id"])
    alt_id = listify(mmcif["_atom_site.label_alt_id"])
    xs = listify(mmcif["_atom_site.Cartn_x"])
    ys = listify(mmcif["_atom_site.Cartn_y"])
    zs = listify(mmcif["_atom_site.Cartn_z"])
    occ = listify(mmcif["_atom_site.occupancy"])
    element = listify(mmcif["_atom_site.type_symbol"])

    rows: list[AtomRow] = []
    for idx in range(len(group)):
        alt = alt_id[idx]
        if alt not in {".", "?", "A", "1"}:
            continue
        rows.append(
            AtomRow(
                group=group[idx],
                atom_id=atom_id[idx],
                comp_id=comp_id[idx],
                auth_asym_id=auth_asym[idx],
                label_asym_id=label_asym[idx],
                auth_seq_id=auth_seq[idx],
                label_seq_id=label_seq[idx],
                alt_id=alt,
                x=float(xs[idx]),
                y=float(ys[idx]),
                z=float(zs[idx]),
                occupancy=occ[idx],
                element=element[idx],
            )
        )
    return rows


def dist(a: AtomRow, b: AtomRow) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def choose_target_atoms(rows: list[AtomRow], target_chain: str) -> list[AtomRow]:
    direct = [row for row in rows if row.auth_asym_id == target_chain]
    if direct:
        return direct
    direct = [row for row in rows if row.label_asym_id == target_chain]
    if direct:
        return direct
    upper = target_chain.upper()
    direct = [row for row in rows if row.auth_asym_id.upper() == upper]
    if direct:
        return direct
    return [row for row in rows if row.label_asym_id.upper() == upper]


def make_site_id(metal: AtomRow) -> str:
    return f"{metal.comp_id}:{metal.auth_asym_id}:{metal.auth_seq_id}:{metal.atom_id}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract Zn/Mg contacts for target PDB chains")
    parser.add_argument("token_list", help="Text file with one PDB_chain token per line")
    parser.add_argument(
        "--mmcif-dir",
        required=True,
        help="Directory with downloaded .cif files",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=3.2,
        help="Distance cutoff in angstroms for contact detection (default: 3.2)",
    )
    parser.add_argument(
        "--metals",
        default=",".join(DEFAULT_METALS),
        help="Comma-separated metal residue names to inspect (default: ZN,MG)",
    )
    parser.add_argument(
        "--contacts-tsv",
        required=True,
        help="Output TSV with per-atom contacts",
    )
    parser.add_argument(
        "--sites-tsv",
        required=True,
        help="Output TSV with per-site residue summaries",
    )
    parser.add_argument(
        "--pymol-dir",
        required=True,
        help="Directory for generated .pml helper scripts",
    )
    args = parser.parse_args()

    metals = {m.strip().upper() for m in args.metals.split(",") if m.strip()}
    mmcif_dir = Path(args.mmcif_dir)
    pymol_dir = Path(args.pymol_dir)
    pymol_dir.mkdir(parents=True, exist_ok=True)

    contact_rows: list[dict[str, str]] = []
    site_rows: list[dict[str, str]] = []

    for token in read_tokens(Path(args.token_list)):
        entry_id, target_chain = parse_token(token)
        mmcif_path = mmcif_dir / f"{entry_id}.cif"
        if not mmcif_path.exists():
            continue

        rows = iter_atom_rows(mmcif_path)
        target_atoms = [
            row for row in choose_target_atoms(rows, target_chain)
            if row.group == "ATOM" and row.comp_id in AA3
        ]
        if not target_atoms:
            continue

        metals_in_structure = [
            row for row in rows
            if row.comp_id.upper() in metals and row.group == "HETATM"
        ]

        site_contacts: dict[str, list[dict[str, str]]] = defaultdict(list)

        for metal in metals_in_structure:
            for atom in target_atoms:
                distance = dist(metal, atom)
                if distance > args.radius:
                    continue
                row = {
                    "token": token,
                    "pdb_id": entry_id,
                    "target_chain": target_chain,
                    "metal": metal.comp_id.upper(),
                    "metal_site_id": make_site_id(metal),
                    "metal_chain": metal.auth_asym_id,
                    "metal_auth_seq_id": metal.auth_seq_id,
                    "metal_atom": metal.atom_id,
                    "resname": atom.comp_id,
                    "auth_seq_id": atom.auth_seq_id,
                    "label_seq_id": atom.label_seq_id,
                    "atom_name": atom.atom_id,
                    "distance": f"{distance:.3f}",
                }
                contact_rows.append(row)
                site_contacts[make_site_id(metal)].append(row)

        for site_id, rows_for_site in sorted(site_contacts.items()):
            residues = []
            seen = set()
            for row in sorted(
                rows_for_site,
                key=lambda r: (
                    r["label_seq_id"] if r["label_seq_id"] not in {".", "?"} else "999999",
                    r["auth_seq_id"],
                    r["atom_name"],
                ),
            ):
                key = (row["resname"], row["auth_seq_id"], row["label_seq_id"])
                if key in seen:
                    continue
                seen.add(key)
                residues.append(f"{row['resname']}{row['auth_seq_id']}[seq={row['label_seq_id']}]")
            exemplar = rows_for_site[0]
            site_rows.append(
                {
                    "token": token,
                    "pdb_id": entry_id,
                    "target_chain": target_chain,
                    "metal": exemplar["metal"],
                    "metal_site_id": site_id,
                    "residue_count": str(len(seen)),
                    "residues": ";".join(residues),
                }
            )

        pml_path = pymol_dir / f"{token}.pml"
        pml_path.write_text(
            "\n".join(
                [
                    f"load {mmcif_path}",
                    f"hide everything, {entry_id}",
                    f"show cartoon, chain {target_chain}",
                    "select metals, resn ZN+MG",
                    "show spheres, metals",
                    f"select near_target, (chain {target_chain}) within {args.radius} of metals",
                    "show sticks, near_target",
                    f"zoom (chain {target_chain}) or metals, 12",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    with Path(args.contacts_tsv).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "token", "pdb_id", "target_chain", "metal", "metal_site_id",
                "metal_chain", "metal_auth_seq_id", "metal_atom",
                "resname", "auth_seq_id", "label_seq_id", "atom_name", "distance",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(contact_rows)

    with Path(args.sites_tsv).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "token", "pdb_id", "target_chain", "metal", "metal_site_id",
                "residue_count", "residues",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(site_rows)

    print(f"Contact rows: {len(contact_rows)}")
    print(f"Site rows: {len(site_rows)}")
    print(f"Contacts TSV: {args.contacts_tsv}")
    print(f"Sites TSV: {args.sites_tsv}")
    print(f"PyMOL scripts dir: {pymol_dir}")


if __name__ == "__main__":
    main()

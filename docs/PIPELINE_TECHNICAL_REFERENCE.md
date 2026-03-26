# Pipeline Technical Reference

This document describes the active analysis path used by `sh/run_full_analysis.sh`, the identifiers passed between stages, the files produced at each step, and the exact role of the active scripts.

## 1. Input Data And Identifiers

Primary local input files:
- `input/cog-24.cog.csv`
  Main COG assignment table.
- `input/cog-24.org.csv`
  Organism metadata keyed by assembly/genome identifiers.
- `input/cog-24.tax.csv`
  Taxonomy metadata.
- `input/2296Genomes.prot.fasta`
  Protein sequences for the COG source database.
- `resources/shortlists/275list.txt`
  Optional representative genome shortlist.

Identifier types used by the pipeline:
- `COG_ID`
  Example: `COG0086`.
- `genome_id`
  Example: `GCF_...` or `GCA_...`; taken from `cog-24.cog.csv` and shortlist files.
- `gene_id`
  Example: locus-like gene identifiers from `cog-24.cog.csv` and `2296Genomes.prot.fasta`.
- FASTA record header
  Usually `gene_id|organism` for extracted COG sequences.
- `PDB_chain` token
  Example: `6GOV_Y`, `1I50_I`; this is the standard structure identifier used inside the pipeline.
- `polymer_entity`
  RCSB entity identifier like `6OVT_1`; used internally during similarity search before conversion to `PDB_chain`.
- `metal_site_id`
  Internal structure-site key like `ZN:D:123:ZN`.
- MSA column coordinate
  1-based column index in the final reordered alignment.

## 2. Master Pipeline Order

The active end-to-end order is:
1. extract COG family sequences,
2. choose the COG FASTA that corresponds to the selected mode,
3. discover or resolve relevant `PDB_chain` tokens,
4. download and clean structure sequences,
5. build the combined sequence + structure MSA,
6. trim alignment columns for tree building,
7. build an IQ-TREE tree,
8. reorder the full combined MSA by that tree,
9. download structure files,
10. detect metal-contact residues in those structures,
11. map structure residues back into the reordered alignment,
12. summarize site coordinates in final MSA space,
13. optionally run the untouched supervisor ligand scan.

The master wrapper is:
- `sh/run_full_analysis.sh`

## 3. Stage 1: COG Sequence Extraction

Wrapper:
- `sh/01_extract_inputs.sh`

Modes:
- `full`
  Script: `scripts/core/cog-extractor-full.py`
  Output: `<COG_ID>.full.raw.fasta`
- `short`
  Script: `scripts/core/cog-extractor-short-list.py`
  Output: `<COG_ID>.short.raw.fasta`
- `representative`
  Script: `scripts/core/12_extract_cog_representatives.py`
  Output: `<COG_ID>.representatives.fasta`
  Report: `<COG_ID>.representatives.tsv`

Shared helper:
- `scripts/core/cog_extractor_common.py`

Data flow:
- `COG_ID` selects rows in `cog-24.cog.csv`.
- rows yield `gene_id` and `genome_id`.
- `gene_id` is matched against `2296Genomes.prot.fasta`.
- output FASTA headers are converted to `gene_id|organism` or `gene_id|genome_id|name`.

## 4. Stage 2: Choosing Structure Tokens

Primary route:
- `sh/00_discover_pdb_tokens.sh`

Scripts:
- `scripts/structure/select_reference_queries.py`
- `scripts/structure/search_rcsb_by_sequence.py`

Files produced:
- `reference_queries.fasta`
- `reference_queries.tsv`
- `pdb_sequence_hits.tsv`
- `pdb_tokens.similarity.txt`

Logic:
- the current COG FASTA is reduced to 1-3 representative query sequences,
- those queries are searched against RCSB by sequence similarity,
- returned `polymer_entity` hits are converted into `PDB_chain` tokens using RCSB entity metadata,
- an E-value and identity cutoff are applied before tokens are accepted.

Fallback route:
- `resources/pdb_lists/<COG_ID>`
- `resources/pdb_lists/ZN<suffix>`
- `resources/pdb_lists/MG<suffix>`

Last fallback:
- `sh/00_resolve_pdb_tokens.sh`
- `scripts/supervisor/extract_cog_pdb_tokens.py`
- `resources/supervisor/COG_to_PDB_chains.txt`

Final provenance file:
- `pdb_token_source.txt`

## 5. Stage 3: Downloading And Cleaning Structure Sequences

Wrapper:
- `sh/02_fetch_and_prepare_structures.sh`

Scripts:
- `scripts/core/getpdb.py`
- `scripts/core/04_prepare_structure_fasta.py`

Inputs:
- `PDB_chain` token list

Outputs:
- `structures.raw.fasta`
- `structures.failed.txt`
- `structures.cleaned.fasta`
- `structures.cleaned.tsv`

Data flow:
- `getpdb.py` fetches FASTA sequences for each `PDB_chain`,
- `04_prepare_structure_fasta.py` normalizes sequences, optionally filters by minimum length, and optionally deduplicates by exact sequence.

## 6. Stage 4: Building The MSA

Wrapper:
- `sh/03_build_msa.sh`

Inputs:
- chosen COG FASTA,
- `structures.cleaned.fasta`

Outputs:
- `combined.input.fasta`
- `structures.only.msa.fasta`
- `combined.msa.fasta`

Logic:
- COG and structure sequences are concatenated into `combined.input.fasta`,
- if the structure FASTA has at least 2 records, a structure-only MSA is built,
- if the combined FASTA has at least 2 records, the final combined MSA is built with `MUSCLE`,
- 0- or 1-record cases are copied through without calling `MUSCLE`.

## 7. Stage 5: Core Alignment, Tree, And Reordered Full MSA

Wrapper:
- `sh/04_tree_and_reorder.sh`

Scripts:
- `scripts/core/14_filter_msa_columns_by_gap.py`
- `scripts/core/09_build_iqtree_tree.py`
- `scripts/core/02_sort_fasta_by_newick.py`

Outputs:
- `combined.core.msa.fasta`
- `combined.core.report.txt`
- `tree.nwk`
- `reordered.fasta`

Data flow:
- columns in `combined.msa.fasta` are filtered by minimum occupancy to create the tree-building alignment,
- IQ-TREE is run on the trimmed alignment,
- the resulting Newick order is applied back to the full combined alignment,
- the final MSA used for interpretation is `reordered.fasta`.

## 8. Stage 6: Structure File Download

Triggered only when structure tokens exist and structure mapping is enabled.

Scripts:
- `scripts/structure/download_rcsb_mmcif.py`
- `scripts/structure/download_rcsb_pdb.py`

Outputs:
- `mmcif/*.cif`
- `pdb/*.pdb`

`mmCIF` is the required parsing source for residue-contact detection. `PDB` files are supplementary and may fail for some entries without stopping the pipeline.

## 9. Stage 7: Detecting Metal-Binding Residues In Structures

Script:
- `scripts/structure/analyze_metal_contacts.py`

Inputs:
- `PDB_chain` token list,
- downloaded `mmCIF` directory,
- metal set from `METALS`, default `ZN,MG`

Outputs:
- `structure_contacts.tsv`
- `structure_sites.tsv`
- `pymol/*.pml`

Data flow:
- each `PDB_chain` is located in its `mmCIF`,
- metal atoms matching the requested ligand names are collected,
- amino-acid atoms within the distance cutoff of each metal are recorded,
- per-contact and per-site summaries are written.

## 10. Stage 8: Mapping Structure Residues Into The Final Alignment

Wrapper:
- `sh/05_map_structure_sites.sh`

Scripts:
- `scripts/structure/map_structure_sequences_to_alignment.py`
- `scripts/structure/map_contacts_via_sequence_matches.py`
- `scripts/structure/summarize_structure_sites_in_alignment.py`

Outputs:
- `structure_sequence_mapping.tsv`
- `contacts_in_reordered.tsv`
- `site_coordinates_in_reordered.txt`
- `site_coordinates_in_reordered.tsv`

Data flow:
- each structure chain sequence is matched to its best sequence in `reordered.fasta`,
- per-residue contact records are transferred from structure residue numbering to MSA column coordinates,
- equivalent sites from multiple structures are merged,
- the final site report is written in human-readable and TSV form.

The final coordinate files are the main bridge between 3D structure analysis and Jalview-style inspection of the family alignment.

## 11. Stage 9: Optional Supervisor Ligand Scan

Wrapper:
- `sh/06_supervisor_ligand_scan.sh`

Scripts:
- `scripts/supervisor/find_ligands_in_PDB.py`
- `scripts/supervisor/udav_pdb.py`
- `scripts/supervisor/udav_soft.py`

Resources:
- `resources/supervisor/req_ligands.txt`
- `resources/supervisor/COG_to_PDB_chains.txt`

Outputs:
- `<prefix>.ligands.txt`
- `<prefix>.ligands.num.txt`
- `<prefix>.log`

This branch is optional, uses Python 2, and is treated as a legacy validation/reporting step rather than the main structure-mapping engine.

## 12. Auxiliary Active Scripts

These scripts are part of the active toolkit but are not always used by the default one-command run:
- `scripts/core/01_filter_msa_by_similarity.py`
  Similarity-based filtering inside an alignment.
- `scripts/core/03_extract_motif_windows.py`
  Extract motif-centered sequence windows.
- `scripts/core/05_filter_structure_msa_by_gaps.py`
  Filter structure-derived MSA fragments by gap content.
- `scripts/core/10_select_one_per_organism.py`
  Select one sequence per organism from a FASTA.
- `scripts/core/11_filter_fasta_records.py`
  Subset FASTA records by ID list or predicates.
- `scripts/core/13_extract_subalignment_by_sequence.py`
  Extract exact sequence subsets from an existing alignment.
- `sh/07_filter_structure_similarity.sh`
  Shell wrapper for similarity filtering.
- `sh/08_extract_motif_windows.sh`
  Shell wrapper for motif-window extraction.

## 13. Legacy Strict Workflow

Wrapper:
- `sh/run_strict_keep_all_cogs.sh`

Purpose:
- older `keep-all-COGs` workflow retained for compatibility with previous `COG0085/COG0086` analyses.

Key scripts:
- `scripts/core/04_prepare_structure_fasta.py`
- `scripts/core/06_filter_msa_keep_all_cogs.py`
- `scripts/core/14_filter_msa_columns_by_gap.py`
- `scripts/core/09_build_iqtree_tree.py`
- `scripts/core/02_sort_fasta_by_newick.py`

Legacy helper scripts not used by the main path:
- `scripts/legacy_helpers/07_build_upgma_tree.py`
- `scripts/legacy_helpers/08_rename_fasta_headers.py`
- `scripts/legacy_helpers/map_contacts_to_alignment.py`
- `scripts/legacy_helpers/pdb_equal.py`

## 14. Minimal Reproducible Execution Order

For the active standard pipeline, the strict file/script order is:
1. `sh/run_full_analysis.sh`
2. `sh/01_extract_inputs.sh`
3. choose `<COG_ID>.*.raw.fasta` or `<COG_ID>.representatives.fasta`
4. `sh/00_discover_pdb_tokens.sh`
5. optionally `sh/00_resolve_pdb_tokens.sh`
6. `sh/02_fetch_and_prepare_structures.sh`
7. `sh/03_build_msa.sh`
8. `sh/04_tree_and_reorder.sh`
9. `scripts/structure/download_rcsb_mmcif.py`
10. `scripts/structure/download_rcsb_pdb.py`
11. `scripts/structure/analyze_metal_contacts.py`
12. `sh/05_map_structure_sites.sh`
13. optionally `sh/06_supervisor_ligand_scan.sh`

The main interpretation files at the end are:
- `reordered.fasta`
- `tree.nwk`
- `site_coordinates_in_reordered.txt`
- `site_coordinates_in_reordered.tsv`

# Project Structure

## Active code

- `scripts/core`
  Main COG extraction, filtering, alignment, tree, and motif utilities.
- `scripts/structure`
  PDB/mmCIF downloaders, metal-contact detection, structure-to-alignment mapping,
  automatic site-coordinate summarization in the final reordered MSA, and
  representative-sequence-driven PDB discovery by similarity search.
- `scripts/supervisor`
  Supervisor-provided helper scripts kept as-is, plus a small Python 3 resolver
  that turns `COG_to_PDB_chains.txt` into per-COG token lists without modifying
  the original supervisor code.

## Non-core helpers

- `scripts/legacy_helpers`
  Older or narrower helpers retained for compatibility:
  - `07_build_upgma_tree.py`
  - `08_rename_fasta_headers.py`
  - `pdb_equal.py`
  - `map_contacts_to_alignment.py`

## Small versioned resources

- `resources/pdb_lists`
  Example PDB-chain token lists.
- `resources/shortlists`
  Genome short-list files.
- `resources/supervisor`
  Supervisor mapping and ligand files used by the automatic PDB-token resolver
  and optional ligand scan.

## Local heavy data

- `input`
  Large local datasets expected by the pipeline and gitignored in normal use.
- `results`
  Output folders from reproducible runs.

## Core documentation

- `README.md`
  GitHub-facing project overview, installation, inputs, outputs, and standard usage.
- `docs/PIPELINE_TECHNICAL_REFERENCE.md`
  Stage-by-stage technical description of identifiers, scripts, files, and data flow.
- `docs/PROJECT_STRUCTURE.md`
  High-level map of the repository.

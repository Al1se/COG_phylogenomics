# COG_phylogenomics

Reproducible phylogenomic pipeline for RNA polymerase-related COG families, with structure-aware site mapping onto final multiple sequence alignments.

The project combines:
- extraction of COG-defined protein families from local NCBI COG tables and proteomes,
- optional selection of one representative per genome,
- automatic discovery of relevant PDB chains by sequence similarity,
- addition of structure-derived sequences to the family alignment,
- tree-based reordering of the final MSA,
- mapping of structure-derived metal-binding residues onto the reordered alignment.

It is designed for workflows such as `COG0085`, `COG0086`, `COG1594`, `COG2093`, but is not hardcoded to those families.

## What The Pipeline Produces

For a standard run, the main wrapper creates:
- extracted COG FASTA,
- structure FASTA downloaded from `PDB_chain` identifiers,
- combined sequence + structure MSA,
- trimmed core alignment for tree building,
- IQ-TREE tree,
- full reordered MSA,
- structure-contact tables,
- site-coordinate reports in the final reordered alignment.

## Repository Layout

- `scripts/core`
  Main family-extraction, alignment, filtering, tree, and FASTA utilities.
- `scripts/structure`
  Structure discovery, download, residue-contact detection, and structure-to-alignment mapping.
- `scripts/supervisor`
  Supervisor-provided legacy tools kept unchanged.
- `scripts/legacy_helpers`
  Older helper scripts retained for compatibility.
- `sh`
  Shell wrappers for each stage and one master runner.
- `resources`
  Versioned small inputs such as shortlist files, curated PDB lists, and supervisor mapping tables.
- `input`
  Large local datasets expected by the pipeline.
- `results`
  Default destination for reproducible analysis runs.
- `docs`
  Project structure notes, legacy material, and technical documentation.

## Installation

Create a Python environment and install the packaged Python dependencies:

```bash
bash sh/install_python_env.sh
source .venv/bin/activate
```

Or install them manually:

```bash
pip install -r requirements.txt
```

Python dependencies in [`requirements.txt`](requirements.txt):
- `biopython` is required by the active structure-aware pipeline.
- `numpy` and `scipy` are only required by the legacy UPGMA helper.

External tools:
- `muscle`
  Required for MSA construction in the active pipeline.
- `iqtree3` or `iqtree2`
  Required for tree building.
- `pymol`
  Optional. Used only for manual inspection of generated `.pml` views.
- `python2`
  Optional. Required only for the untouched supervisor tools.

Example `conda` installation:

```bash
conda install -c conda-forge muscle iqtree biopython numpy scipy pymol-open-source
```

## Windows / WSL Setup

The project is designed around `bash` wrappers in `sh/`, so the recommended Windows path is:
- run the pipeline inside `WSL` (preferably Ubuntu),
- keep the repository and large input files inside the Linux filesystem,
- optionally open results from Windows editors or viewers.

Native `cmd.exe` / PowerShell execution is not the supported path for the current repository layout.

Short version for Windows users:

1. In PowerShell:

```powershell
wsl --install -d Ubuntu
```

2. In Ubuntu after cloning the repository:

```bash
bash sh/install_wsl_ubuntu.sh
source .venv/bin/activate
```

3. Put the required input files into `input/`, then run:

```bash
bash sh/run_full_analysis.sh COG0086 results/COG0086 short
```

Full Windows instructions are in [`docs/WINDOWS_WSL_SETUP.md`](docs/WINDOWS_WSL_SETUP.md).

Recommended first-time setup on Windows:

1. Install WSL and Ubuntu:

```powershell
wsl --install -d Ubuntu
```

2. Open Ubuntu and install the base system tools:

```bash
sudo apt update
sudo apt install -y python3 python3-venv python3-pip muscle iqtree git curl
```

3. Clone the repository inside WSL:

```bash
cd ~
git clone <YOUR_REPOSITORY_URL> COG_phylogenomics
cd COG_phylogenomics
```

4. Create the Python environment:

```bash
bash sh/install_wsl_ubuntu.sh
source .venv/bin/activate
```

5. Put the required large input files into `input/`:

```bash
ls input
```

Expected files:
- `cog-24.cog.csv`
- `cog-24.org.csv`
- `cog-24.tax.csv`
- `2296Genomes.prot.fasta`

6. Run the pipeline:

```bash
bash sh/run_full_analysis.sh COG0086 results/COG0086 short
```

Optional notes for Windows users:
- For best performance, keep the repository under the WSL filesystem, for example in `~/COG_phylogenomics`, not under `/mnt/c/...`.
- `PyMOL` is optional. The pipeline itself does not require it. If you want GUI visualization, the easiest route is often to run the analysis in WSL and inspect generated `.pml` / structure files later in native Windows PyMOL.
- On Windows 11 with `WSLg`, Linux GUI applications can also work, but this is not required for the main pipeline.

## Required Local Input Data

The active pipeline expects these files in `input/`:
- `cog-24.cog.csv`
- `cog-24.org.csv`
- `cog-24.tax.csv`
- `2296Genomes.prot.fasta`

These large files are gitignored even if local copies are present in the working tree.

## One-Command Analysis

Standard run:

```bash
bash sh/run_full_analysis.sh COG0086 results/COG0086 short
```

Explicit structure token list:

```bash
bash sh/run_full_analysis.sh COG0086 results/COG0086 short resources/pdb_lists/ZN0086
```

Available extraction modes:
- `full`
  All records for the requested COG from the local COG table.
- `short`
  All records for the requested COG restricted to the shortlist genomes.
- `representative`
  One selected record per shortlist genome, optionally motif-prioritized.

## How Structure Chains Are Chosen

If an explicit `pdb_token_list` is not provided, `run_full_analysis.sh` resolves structures in this order:
1. sequence-similarity search against RCSB PDB using representative query sequence(s) selected from the current COG FASTA,
2. curated token lists in `resources/pdb_lists/`,
3. supervisor mapping from `resources/supervisor/COG_to_PDB_chains.txt`.

Files written during this process:
- `reference_queries.fasta`
- `reference_queries.tsv`
- `pdb_sequence_hits.tsv`
- `pdb_tokens.similarity.txt`
- `pdb_tokens.auto.txt`
- `pdb_token_source.txt`

## Main Outputs Of A Standard Run

Typical important files in `results/<RUN_NAME>/`:
- `<COG_ID>.*.raw.fasta`
  Extracted family sequences.
- `structures.raw.fasta`
  Downloaded PDB-chain sequences.
- `structures.cleaned.fasta`
  Cleaned and optionally deduplicated structure sequences.
- `combined.input.fasta`
  Union of family and structure sequences.
- `combined.msa.fasta`
  Full MSA used as the final sequence space.
- `combined.core.msa.fasta`
  Gap-filtered alignment used for tree building.
- `tree.nwk`
  IQ-TREE result.
- `reordered.fasta`
  Final full MSA reordered by the tree.
- `structure_contacts.tsv`
  Per-atom structure-contact table.
- `structure_sites.tsv`
  Per-site structure summary.
- `contacts_in_reordered.tsv`
  Structure contacts mapped to the final reordered alignment.
- `site_coordinates_in_reordered.txt`
  Human-readable site summary in final MSA coordinates.
- `site_coordinates_in_reordered.tsv`
  Tabular version of the same mapping.

## Main Shell Entry Points

- `sh/install_python_env.sh`
  Create a virtual environment and install Python dependencies.
- `sh/01_extract_inputs.sh`
  Extract COG sequences in `full`, `short`, or `representative` mode.
- `sh/00_discover_pdb_tokens.sh`
  Build representative query sequence(s) and search RCSB PDB by sequence similarity.
- `sh/00_resolve_pdb_tokens.sh`
  Resolve a per-COG PDB token list from supervisor mapping files.
- `sh/02_fetch_and_prepare_structures.sh`
  Download PDB-chain FASTA and clean structure records.
- `sh/03_build_msa.sh`
  Build structure-only and combined MSA files.
- `sh/04_tree_and_reorder.sh`
  Trim alignment columns, build an IQ-TREE tree, and reorder the full MSA.
- `sh/05_map_structure_sites.sh`
  Map structure-derived residues onto the reordered MSA and summarize site coordinates.
- `sh/06_supervisor_ligand_scan.sh`
  Optional launch of the untouched supervisor ligand scanner.
- `sh/run_full_analysis.sh`
  Main end-to-end pipeline.
- `sh/run_strict_keep_all_cogs.sh`
  Legacy strict workflow retained for older analyses.

## Environment Variables

Useful overrides:

```bash
export MUSCLE_BIN=muscle
export IQTREE_BIN=iqtree3
export CORE_OCCUPANCY=0.7
export STRUCT_MIN_LENGTH=0
export STRUCT_DEDUP=yes
export RUN_STRUCTURE_MAPPING=1
export RUN_SUPERVISOR_SCAN=auto
export METALS=ZN,MG
export REFERENCE_QUERY_COUNT=3
export REFERENCE_FASTA=/path/to/reference_queries.fasta
export PDB_SEARCH_EVALUE=1e-10
export PDB_SEARCH_IDENTITY=0.30
export PDB_SEARCH_MAX_HITS=100
export SUPERVISOR_PYTHON_BIN=python2.7
```

If `IQTREE_BIN` is not set, the wrappers try:
- executable from `PATH`,
- `/tmp/iqtree-env/bin/iqtree3`,
- `iqtree2` from `PATH`.

## Technical Reference

For the exact order of stages, identifier flow, input/output file semantics, and the role of every active script, see [`docs/PIPELINE_TECHNICAL_REFERENCE.md`](docs/PIPELINE_TECHNICAL_REFERENCE.md).

## Supervisor And Legacy Components

These files were copied without code changes and are treated as external legacy tools:
- `scripts/supervisor/udav_pdb.py`
- `scripts/supervisor/udav_soft.py`
- `scripts/supervisor/find_ligands_in_PDB.py`
- `resources/supervisor/COG_to_PDB_chains.txt`

Older notes and runbooks remain under `docs/legacy/`.

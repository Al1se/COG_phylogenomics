# COG0085/COG0086 Pipeline Usage

Minimal local workflow for the coursework project.

For strict requirements and a full serious run scenario, see:
- `TZ_STRICT_COG_PIPELINE.txt`
- `RUNBOOK_COG0086_SERIOUS.txt`

## 1) Extract COG sequences

Full dataset:

```bash
python3 cog-extractor-full.py COG0085 -o COG0085.raw.fasta
python3 cog-extractor-full.py COG0086 -o COG0086.raw.fasta
```

Short-list dataset (275 genomes):

```bash
python3 cog-extractor-short-list.py COG0085 -o 275_COG0085.raw.fasta
python3 cog-extractor-short-list.py COG0086 -o 275_COG0086.raw.fasta
```

## 2) Fetch PDB reference sequences

```bash
python3 getpdb.py ZN0085 -o ZN0085.fasta --failed ZN0085.failed.txt
python3 getpdb.py ZN0086 -o ZN0086.fasta --failed ZN0086.failed.txt
```

## 3) Merge local FASTA files

Use shell concatenation (simple and explicit):

```bash
cat COG0085.raw.fasta ZN0085.fasta > COG0085.with_pdb.fasta
cat COG0086.raw.fasta ZN0086.fasta > COG0086.with_pdb.fasta
```

## 4) Run MSA (MUSCLE)

Confirmed launch format:

```bash
muscle -align ./ZN0085.fasta -output ./test_pdb.fasta
muscle -align ./ZN0086.fasta -output ./test_pdb.fasta
```

```bash
muscle -align ./COG0085.with_pdb.fasta -output ./COG0085.msa.fasta
muscle -align ./COG0086.with_pdb.fasta -output ./COG0086.msa.fasta
```

## 5) Filter MSA by similarity profile

Example with references and report:

```bash
python3 01_filter_msa_by_similarity.py \
  COG0085.msa.fasta \
  -o COG0085.msa.filtered.fasta \
  --reference-id 5TJG_D --reference-id 6C6T_J \
  --threshold 0.45 \
  --report COG0085.filter.tsv
```

```bash
python3 01_filter_msa_by_similarity.py \
  COG0086.msa.fasta \
  -o COG0086.msa.filtered.fasta \
  --reference-id 5TJG_D --reference-id 6C6T_J \
  --threshold 0.45 \
  --report COG0086.filter.tsv
```

## 6) Build tree externally and sort FASTA by Newick leaf order

After downloading/exporting Newick from ngphylogeny.fr:

```bash
python3 02_sort_fasta_by_newick.py \
  COG0085.msa.filtered.fasta COG0085.tree.nwk \
  -o COG0085.msa.filtered.reordered.fasta
```

```bash
python3 02_sort_fasta_by_newick.py \
  COG0086.msa.filtered.fasta COG0086.tree.nwk \
  -o COG0086.msa.filtered.reordered.fasta
```

By default, `02_sort_fasta_by_newick.py` now fails if the tree and FASTA do not contain the same IDs.
Use `--allow-partial-match` only for exploratory runs.

## 7) Extract motif windows

Example: keep windows around motif with 2 aa flanks (`XX[motif]XX` logic):

```bash
python3 03_extract_motif_windows.py \
  COG0086.msa.filtered.reordered.fasta \
  -o COG0086.NADFDGD.windows.fasta \
  --motif NADFDGD --left 2 --right 2
```

Wildcard motif example (`X` means any residue):

```bash
python3 03_extract_motif_windows.py \
  COG0085.msa.filtered.reordered.fasta \
  -o COG0085.XXAAXX.windows.fasta \
  --motif XXAAXX --left 2 --right 2 --wildcard X
```

Note: `XXAAXX` is only a technical wildcard example, not a scientific motif to interpret in conclusions.

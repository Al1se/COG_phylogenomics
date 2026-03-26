# Windows / WSL Setup

This project is built around `bash` wrappers and Linux command-line tools. The recommended Windows workflow is:
- install `WSL` with Ubuntu,
- clone the repository inside the Linux filesystem,
- run all analysis commands inside Ubuntu,
- optionally inspect results later from Windows tools.

Running the project directly from native `cmd.exe` or pure PowerShell is not the supported path.

## Why WSL Is Recommended

The active pipeline depends on:
- `bash`,
- `python3` with `venv`,
- `muscle`,
- `iqtree`,
- standard Linux-style file paths and shell wrappers.

`WSL` provides all of that with minimal adaptation.

## One-Time Setup

### 1. Install WSL and Ubuntu

Run in Windows PowerShell:

```powershell
wsl --install -d Ubuntu
```

Reboot if Windows asks for it, then open the new Ubuntu terminal and finish the initial user setup.

### 2. Clone The Repository Inside WSL

In Ubuntu:

```bash
cd ~
git clone <YOUR_REPOSITORY_URL> COG_phylogenomics
cd COG_phylogenomics
```

Recommended:
- keep the repository under `~/...` inside WSL,
- do not keep the working copy under `/mnt/c/...` if you want stable performance on large FASTA files.

### 3. Install All Required Runtime Dependencies

In Ubuntu:

```bash
bash sh/install_wsl_ubuntu.sh
source .venv/bin/activate
```

This installs:
- `python3`
- `python3-venv`
- `python3-pip`
- `muscle`
- `iqtree`
- Python packages from `requirements.txt`

Optional components not installed by default:
- `python2`
  only needed for the untouched supervisor scripts
- `pymol`
  only needed for manual structure inspection

If needed, they can be installed later.

## Required Input Files

Place these files into `input/`:
- `cog-24.cog.csv`
- `cog-24.org.csv`
- `cog-24.tax.csv`
- `2296Genomes.prot.fasta`

Quick check:

```bash
ls input
```

## First Test Run

After activation of the Python environment:

```bash
bash sh/run_full_analysis.sh COG0086 results/COG0086 short
```

This should create a new run directory under `results/COG0086/`.

## Typical Daily Usage

After the one-time installation, start Ubuntu and run:

```bash
cd ~/COG_phylogenomics
source .venv/bin/activate
bash sh/run_full_analysis.sh COG0086 results/COG0086 short
```

## Optional GUI Tools

`PyMOL` is not required for the pipeline itself.

Two practical options:
- run the pipeline in WSL and open the generated structure files later in native Windows PyMOL,
- or install a Linux `PyMOL` under WSL/WSLg if your Windows setup supports Linux GUI applications.

## Troubleshooting

### `iqtree3` not found

Check:

```bash
which iqtree
which iqtree2
which iqtree3
```

The wrappers can work with either `iqtree3` or `iqtree2`.

### `muscle` not found

Check:

```bash
which muscle
```

### Slow filesystem performance

If the repository lives under `/mnt/c/...`, move it into the Linux home directory:

```bash
cd ~
mv /mnt/c/path/to/COG_phylogenomics .
```

### Python environment not active

Activate it:

```bash
source .venv/bin/activate
```

## Minimal “Out Of The Box” Command Set

After cloning into WSL, the intended short path is:

```bash
bash sh/install_wsl_ubuntu.sh
source .venv/bin/activate
```

Then, after copying the input files into `input/`:

```bash
bash sh/run_full_analysis.sh COG0086 results/COG0086 short
```

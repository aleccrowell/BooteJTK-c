# BooteJTK

[![PyPI version](https://img.shields.io/pypi/v/bootjtk.svg)](https://pypi.org/project/bootjtk/)
[![Python Versions](https://img.shields.io/pypi/pyversions/bootjtk.svg)](https://pypi.org/project/bootjtk/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

BooteJTK is an implementation of empirical JTK (eJTK) on parametrically bootstrapped resamplings of time series, used for detecting circadian rhythms in genomic data.

Based on [BooteJTK](https://github.com/alanlhutchison/BooteJTK) by Alan Hutchison _et al._; this fork improves Python 3 compatibility and integration with [LIMBR](https://github.com/aleccrowell/LIMBR).

**References**

- Hutchison AL _et al._ (2016), "BooteJTK: Improved Rhythm Detection via Bootstrapping", bioRxiv.
- Hutchison AL, Maienschein-Cline M, Chiang AH _et al._ "Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data." _PLoS Computational Biology_ 2015 11(3): e1004094. [doi:10.1371/journal.pcbi.1004094](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004094)


## Installation

```
pip install bootjtk
```

Requires Python 3.8 or later. All dependencies (numpy, scipy, pandas, matplotlib, statsmodels) are installed automatically. No R installation is required.


## Quick start

Run 10 bootstrap resamplings on data with 2 replicates per timepoint:

```
bootejtk-calcp -f example/TestInput4.txt -x MYPREFIX -r 2 -z 10
```

The `-p` (period), `-s` (phases), and `-a` (asymmetries) ref-file arguments default to the standard 24 h files bundled with the package, so they can be omitted for typical circadian analyses.


## Usage

### `bootejtk-calcp` — full pipeline

This is the main entry point. It runs BooteJTK bootstrapping followed by p-value calculation.

```
bootejtk-calcp -f <input_file> -x <prefix> -r <replicates> -z <bootstraps> [options]
```

| Option | Description | Default |
|---|---|---|
| `-f` / `--filename` | Input data file (tab-delimited, header row starting with `#` or `ID`) | required |
| `-x` / `--prefix` | Output file prefix | required |
| `-r` / `--reps` | Number of replicates per timepoint | `2` |
| `-z` / `--size` | Number of bootstrap resamplings | `50` |
| `-j` / `--workers` | Worker processes (`0` = all CPUs) | `1` |
| `-w` / `--waveform` | Reference waveform shape (see below) | `cosine` |
| `-p` / `--period` | Period reference file | bundled 24 h file |
| `-s` / `--phase` | Phase reference file | bundled 0–22 h by 2 file |
| `-a` / `--width` | Asymmetry reference file | bundled 2–22 h by 2 file |
| `-B` / `--basic` | Skip variance shrinkage preprocessing | off |
| `-L` / `--limma` | Use limma vooma variance estimation only (without vash imputation) | off |
| `--vash` | Use vash NA imputation before variance estimation | off |
| `-U` / `--noreps` | No replicates mode: estimate variance from arrhythmic genes | off |
| `-R` / `--rnaseq` | RNA-seq mode (passed to preprocessing) | off |
| `-W` / `--write` | Write pickle output files (`.pkl`) from BooteJTK | off |

Run `bootejtk-calcp --help` to see all options and current defaults.

### Preprocessing / variance shrinkage

By default, `bootejtk-calcp` applies vooma-style variance estimation with NA imputation and empirical Bayes shrinkage (the equivalent of `--limma --vash`). This preprocessing is implemented entirely in Python — no R installation is required.

| Flag | Behaviour |
|---|---|
| *(default, no flag)* | Run limma vooma + vash NA imputation + eBayes shrinkage |
| `-L` / `--limma` | Run limma vooma + eBayes shrinkage (no NA imputation) |
| `--vash` | Enable vash NA imputation (implied by default) |
| `-B` / `--basic` | Skip all variance shrinkage; use raw values directly |

### `bootejtk` — core analysis only

Runs the BooteJTK analysis step without the CalcP p-value fitting step. Useful if you want to run CalcP separately or with custom settings.

```
bootejtk -f <input_file> -x <prefix> -r <replicates> -z <bootstraps> [options]
```

### Waveform shapes

| Value | Shape |
|---|---|
| `cosine` (default) | Smooth sinusoidal peak |
| `trough` | Triangular trough |
| `impulse` | Narrow spike |
| `step` | Rectangular step |

### Parallel processing

Use `-j` to speed up large datasets by distributing genes across CPUs:

```
bootejtk-calcp -f example/TestInput4.txt -r 2 -z 50 -j 8
```

| `-j` value | Behaviour |
|---|---|
| `1` (default) | Sequential, single process |
| `N > 1` | Use N worker processes |
| `0` | Use all available CPUs |


## Input format

Tab-delimited text file. The header row must start with `#` or `ID`; subsequent columns are zeitgeber time labels. Each data row begins with a gene/feature identifier.

The preferred header format encodes both the timepoint and replicate number:

```
#	ZT00_1	ZT00_2	ZT02_1	ZT02_2	ZT04_1	ZT04_2	...
gene1	1.23	1.31	2.45	2.38	3.10	3.05	...
gene2	5.01	4.95	4.87	4.91	3.92	4.01	...
```

Column labels follow the pattern `ZT{HH}_{rep}` where `HH` is the zero-padded hour (e.g. `00`, `02`, `04`) and `rep` is the replicate number. This format is shared with [LIMBR](https://github.com/aleccrowell/LIMBR) and [PIRS](https://github.com/aleccrowell/PIRS), so output from those tools can be passed directly to BooteJTK without reformatting.

The legacy formats `ZT0`, `ZT2`, `CT0`, `CT2` (no replicate suffix) remain fully supported for backwards compatibility, as do decimal values (e.g. `ZT14.7`). Time labels do not need to be evenly spaced.


## Output files

Running `bootejtk-calcp` produces output files prefixed with the value passed to `-x`:

| File | Contents |
|---|---|
| `*_GammaP.txt` | BooteJTK output with Gamma-fitted p-values |
| `*.txt` | Main BooteJTK output (best-matching waveform per gene, feeds into CalcP) |
| `*_NULL1000.txt` | Randomly generated null time series used to fit the null tau distribution |
| `*_order_probs.pkl` | *(requires `-W`)* Pickle: per-gene `[means, stds, ns]` and rank-order bootstrap frequencies |
| `*_order_probs_vars.pkl` | *(requires `-W`)* Pickle: per-gene tau and phase probability distributions |

> Running the example command on an already-existing output directory appends `_1` to output filenames.


## FAQ

**Can I use non-integer or uneven time intervals (e.g. ZT14.7)?**
Yes. The label just needs to start with `ZT` or `CT`; decimal values are read correctly.

**Does BooteJTK handle uneven sampling intervals?**
Yes. All timepoints in the header are used as given.

**Why does BooteJTK report phases like 14.4 that don't match my sampling intervals?**
BooteJTK runs bootstrap resamplings and reports the *mean* phase across those resamplings. For example, if 8 of 10 resamplings give phase 14 and 2 give phase 16, the reported mean phase is 14.4.

**Do the phase/asymmetry search intervals need to match the sampling intervals?**
No. You can sample every hour but only search for phases every two hours, for example.


## Development

```
git clone https://github.com/aleccrowell/BooteJTK-c
cd BooteJTK-c
pip install poetry
poetry install
poetry run pytest tests/ -v
```


## License

Released under the MIT License. See `LICENSE` for details.

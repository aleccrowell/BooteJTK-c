# Note

This is a fork of [BooteJTK](https://github.com/alanlhutchison/BooteJTK) by Alan Hutchison, intended to improve compatibility with [LIMBR](https://github.com/aleccrowell/LIMBR).

# BooteJTK

BooteJTK is an implementation of empirical JTK (eJTK) on parametrically bootstrapped resamplings of time series, used for detecting circadian rhythms in genomic data.

Information on BooteJTK can be found in [Hutchison AL _et al._ (2016), "BooteJTK: Improved Rhythm Detection via Bootstrapping"](), available at bioRxiv.

Information on eJTK can be found in [Hutchison AL, Maienscein-Cline M, Chiang AH, _et al._ "Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data." PLoS Computational Biology 2015 Mar. Vol. 11, No. 3, pp. e1004094. doi:10.1371/journal.pcbi.1004094](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004094)


## Requirements

* Python >= 3.8
* numpy >= 1.11.0
* scipy >= 0.15.1
* cython >= 0.24
* statsmodels >= 0.6.1
* matplotlib >= 2.0.0

Install dependencies with [Poetry](https://python-poetry.org/):

```
poetry install
```

Or with pip:

```
pip install numpy scipy cython statsmodels matplotlib
```


## Build

The Kendall's tau and waveform-matching core is a compiled Cython extension. Build it before first use:

```
cd bin
python setup.py build_ext --inplace
```

This produces `bin/get_stat_probs.cpython-*.so`. You must rebuild after any changes to `bin/get_stat_probs.pyx`.


## Running

### Basic example

This command runs 10 bootstrap resamplings on data with 2 replicates per timepoint. It searches phases 0–22h in 2h steps, asymmetries 2–22h in 2h steps, and a 24h period.

```
./BooteJTK-CalcP.py -f example/TestInput4.txt \
    -p ref_files/period24.txt \
    -s ref_files/phases_00-22_by2.txt \
    -a ref_files/asymmetries_02-22_by2.txt \
    -x OTHERTEXT -r 2 -z 10
```

### Parallel processing

Use `-j` to run genes in parallel across multiple CPUs. This is the largest available speedup for datasets with many genes:

```
./BooteJTK-CalcP.py -f example/TestInput4.txt \
    -p ref_files/period24.txt \
    -s ref_files/phases_00-22_by2.txt \
    -a ref_files/asymmetries_02-22_by2.txt \
    -x OTHERTEXT -r 2 -z 50 -j 8
```

| `-j` value | Behaviour |
|---|---|
| `1` (default) | Sequential, single process |
| `N > 1` | Use N worker processes |
| `0` | Use all available CPUs |


## Output files

1. **`*_GammaP.txt`** — BooteJTK output with a Gamma fit of the null dataset, used to assign p-values to Tau values.

2. **`*.txt`** — Main BooteJTK output. Contains the best-matching reference waveform for each time series (highest absolute Tau). This feeds into `CalcP.py`, which runs automatically.

3. **`*_order_probs.pkl`** — Pickle of two dictionaries:
   - Keys are gene IDs; values are lists of `[means, standard_deviations, replicate_counts]` per timepoint.
   - Keys are gene IDs; values are dicts mapping rank-ordered bootstrap samples to their frequency.

4. **`*_order_probs_vars.pkl`** — Pickle of two dictionaries:
   - Keys are gene IDs; values are dicts mapping Tau values to their probability.
   - Keys are gene IDs; values are dicts mapping phase values to their probability.

5. **`*_NULL1000.txt`** — Randomly generated time series used to fit the null Tau distribution. See `example/TestInput4_NULL1000.txt` for a reference.

> Running the example command as-is will append `_1` to output filenames since the example outputs already exist.


## Tests

```
python -m pytest tests/ -v
```


## License

Released under the MIT License. See `LICENSE.txt` for details.

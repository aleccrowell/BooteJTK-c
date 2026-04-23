# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.2.0] — 2026-04-22

### Changed
- Header parsers (`get_data`, `get_data2`, `get_data_multi`, `parse_timepoint_label`) now
  strip a `_{rep}` replicate suffix before converting the numeric portion, so columns in
  the shared `ZT{HH}_{rep}` format (e.g. `ZT02_1`, `ZT02_2`) are parsed correctly. The
  legacy `ZTn` / `CTn` formats remain fully supported.

## [1.1.0] — 2026-04-22

### Added
- `bootjtk/limma_voom.py`: pure-Python replacement for the `Limma_voom_core.R`
  and `Limma_voom_vash_core.R` R scripts, implementing vooma-style variance
  estimation and Smyth (2004) empirical Bayes shrinkage. The package now has
  **no R dependency** at any point in the pipeline.
- End-to-end integration tests with synthetically generated circadian
  time-series, covering rhythm detection accuracy (phase, period, tau) and the
  full BooteJTK → CalcP p-value pipeline (249 tests total).
- `bootjtk.__version__` attribute sourced from package metadata.

### Changed
- `pipeline.py` (`bootejtk-calcp`) no longer shells out to `Rscript`; the
  `--limma` and `--vash` preprocessing paths now call the new Python functions
  directly.
- Removed all eight legacy R scripts (`Limma_*.R`) from the package.
- Removed `write_preprocessed` helper from `limma_preprocess.py` (no longer
  needed now that data is not written to temp files for R).

### Fixed
- Python 3 compatibility bug in `BooteJTK.eBayes`: `np.array(dict.values())`
  produced a 0-dimensional object array in Python 3; wrapped with `list()`.

## [1.0.0] — 2026-04-21

First release on PyPI.

### Added
- `bootjtk` installable Python package with `bootejtk` and `bootejtk-calcp` console-script entry points.
- `bootjtk/get_stat_probs.py`: pure-Python replacement for the former Cython extension, using `scipy.stats.kendalltau` for Kendall's tau computation. No compiled extension required.
- `bootjtk/pipeline.py`: installable entry point for the full BooteJTK-CalcP pipeline (formerly the standalone `BooteJTK-CalcP.py` script).
- Parallel gene processing via `-j` / `--workers` argument (multiprocessing).
- Test suite covering core statistical functions, waveform generation, bootstrapping, and CLI utilities (144 tests).

### Changed
- Package source moved from `bin/` to `bootjtk/`; `ref_files/` moved inside `bootjtk/` so they are included in the installed package.
- All intra-package imports converted to relative imports (`from .module import ...`).
- `_REF_DIR` and R-script paths now resolved relative to `__file__` so they work correctly after installation.
- Dropped Cython as a dependency; `scipy` (already a core dependency) provides an equivalent `kendalltau` implementation.
- `pyproject.toml` updated with PyPI classifiers, keywords, homepage/repository URLs, and `[tool.poetry.scripts]` entry points.
- Renamed `LICENSE.txt` → `LICENSE`.

### Fixed
- Python 3 compatibility: removed `xrange` calls, fixed bare `print` statements, and ensured `map()` results are consumed as lists where needed.

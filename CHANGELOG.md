# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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

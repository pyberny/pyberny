# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.7.0] - 2026-07-17

### Added

- `symmetry` argument to `Berny`: point-group detection and opt-in symmetry breaking of start geometries ([#162](https://github.com/pyberny/pyberny/pull/162))
- `molsym` as a required dependency, backing the `symmetry` feature ([#162](https://github.com/pyberny/pyberny/pull/162), [#184](https://github.com/pyberny/pyberny/pull/184))
- NumPy 2 support ([#43](https://github.com/pyberny/pyberny/pull/43))
- Covalent radii for Ce–Yb, Po, At, and Fr–U ([#45](https://github.com/pyberny/pyberny/pull/45))
- `berny.solvers.XTBSolver`, a GFN-xTB backend via `tblite` (default GFN2-xTB) ([#139](https://github.com/pyberny/pyberny/pull/139))
- Linear-bend internal coordinates via dummy atoms ([#53](https://github.com/pyberny/pyberny/pull/53))
- `Ghost`, `X`, and `Bq` basis-function-only centres ([#53](https://github.com/pyberny/pyberny/pull/53))
- `berny.tests` subpackage of reusable, optimizer-agnostic end-to-end tests ([#101](https://github.com/pyberny/pyberny/pull/101))
- Benchmark suite: `scripts/benchmark.py` runner and a `berny.benchmarks` subpackage bundling the Birkholz–Schlegel ([#55](https://github.com/pyberny/pyberny/pull/55)), Baker ([#84](https://github.com/pyberny/pyberny/pull/84)), and oligomer ([#127](https://github.com/pyberny/pyberny/pull/127)) sets with a discovery API ([#116](https://github.com/pyberny/pyberny/pull/116))
- `energy_noise` optimizer parameter for noise-aware convergence ([#108](https://github.com/pyberny/pyberny/pull/108))
- `trace` argument to `Berny` writing a structured per-step JSON trace ([#112](https://github.com/pyberny/pyberny/pull/112))
- PEP 561 `py.typed` marker; fully typed under `mypy --strict` ([#119](https://github.com/pyberny/pyberny/pull/119))

### Changed

- `berny.Math.FindrootException` renamed to `berny.Math.FindrootError` ([#42](https://github.com/pyberny/pyberny/pull/42))
- Unknown keyword arguments to `Berny()` now raise `TypeError` ([#51](https://github.com/pyberny/pyberny/pull/51))
- Mid-run coordinate rebuilds preserve accumulated Hessian curvature ([#122](https://github.com/pyberny/pyberny/pull/122))
- Linear-bend rebuild also fires at higher-coordination centres ([#104](https://github.com/pyberny/pyberny/pull/104))
- Repository and documentation moved to the `pyberny` GitHub organization ([#183](https://github.com/pyberny/pyberny/pull/183))

### Removed

- Support for Python older than 3.10 ([#119](https://github.com/pyberny/pyberny/pull/119))
- The runtime `setuptools` (`pkg_resources`) dependency ([#42](https://github.com/pyberny/pyberny/pull/42))
- `berny.solvers.MopacSolver` and the `"mopac"` geometry output format (use `XTBSolver`) ([#142](https://github.com/pyberny/pyberny/pull/142))
- The module-level `berny.berny.defaults` dict; pass overrides as `Berny()` keyword arguments instead ([#51](https://github.com/pyberny/pyberny/pull/51))

### Fixed

- `get_property` raises a clear `KeyError` for missing species data ([#45](https://github.com/pyberny/pyberny/pull/45))
- Species lookup by atomic number ([#45](https://github.com/pyberny/pyberny/pull/45))

## [0.6.3] - 2021-02-22

### Fixed

- CLI

[0.7.0]: https://github.com/pyberny/pyberny/compare/0.6.3...0.7.0
[0.6.3]: https://github.com/pyberny/pyberny/releases/tag/0.6.3

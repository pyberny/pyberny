# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.7.0] - 2026-07-17

### Added

- `symmetry` argument to `Berny`: point-group detection and opt-in symmetry breaking of start geometries [#162]
- NumPy 2 support [#43]
- Covalent radii for Ce–Yb, Po, At, and Fr–U [#45]
- `berny.BernyParams` dataclass exposing all optimizer parameters [#51]
- `berny.solvers.XTBSolver`, a GFN-xTB backend via `tblite` (default GFN2-xTB) [#139]
- Linear-bend internal coordinates via dummy atoms [#53]
- `Ghost`, `X`, and `Bq` basis-function-only centres [#53]
- `berny.tests` subpackage of reusable, optimizer-agnostic end-to-end tests [#101]
- `berny.benchmarks` subpackage bundling the Birkholz–Schlegel [#55], Baker [#84], and oligomer [#127] sets with a discovery API [#116]
- `scripts/benchmark.py` benchmark runner [#55]
- Interactive 3D viewer of the benchmark molecules in the docs [#106], [#180]
- PEP 561 `py.typed` marker; fully typed under `mypy --strict` [#119]

### Changed

- Minimum supported Python raised to 3.10 [#119]
- `berny.Math.FindrootException` renamed to `berny.Math.FindrootError` [#42]
- Dropped the runtime `setuptools` (`pkg_resources`) dependency [#42]
- Unknown keyword arguments to `Berny()` now raise `TypeError` [#51]
- Mid-run coordinate rebuilds preserve accumulated Hessian curvature [#122]
- Linear-bend rebuild also fires at higher-coordination centres [#104]

### Removed

- `berny.solvers.MopacSolver` (use `XTBSolver`) [#142]
- The `"mopac"` geometry output format [#142]
- The module-level `berny.berny.defaults` dict (use `BernyParams`) [#51]

### Fixed

- `get_property` raises a clear `KeyError` for missing species data [#45]
- Species lookup by atomic number [#45]

## [0.6.3] - 2021-02-22

### Fixed

- CLI

[0.7.0]: https://github.com/pyberny/pyberny/compare/0.6.3...0.7.0
[0.6.3]: https://github.com/pyberny/pyberny/releases/tag/0.6.3
[#42]: https://github.com/pyberny/pyberny/pull/42
[#43]: https://github.com/pyberny/pyberny/pull/43
[#45]: https://github.com/pyberny/pyberny/pull/45
[#51]: https://github.com/pyberny/pyberny/pull/51
[#53]: https://github.com/pyberny/pyberny/pull/53
[#55]: https://github.com/pyberny/pyberny/pull/55
[#84]: https://github.com/pyberny/pyberny/pull/84
[#101]: https://github.com/pyberny/pyberny/pull/101
[#104]: https://github.com/pyberny/pyberny/pull/104
[#106]: https://github.com/pyberny/pyberny/pull/106
[#116]: https://github.com/pyberny/pyberny/pull/116
[#119]: https://github.com/pyberny/pyberny/pull/119
[#122]: https://github.com/pyberny/pyberny/pull/122
[#127]: https://github.com/pyberny/pyberny/pull/127
[#139]: https://github.com/pyberny/pyberny/pull/139
[#142]: https://github.com/pyberny/pyberny/pull/142
[#162]: https://github.com/pyberny/pyberny/pull/162
[#180]: https://github.com/pyberny/pyberny/pull/180

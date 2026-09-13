# APPLE 1.0.0

- Package the supplied September 2026 source with explicit dependencies and
  generated function documentation. Loading the package no longer installs
  or attaches dependencies.
- Export the current interfaces: `Remove.IP()`, `Map.Tail()`, `Tail.PCA()`,
  and `Tail.DiffPair()`.
- Import `makeTxDbFromGFF()` from `txdbmaker`, and declare data.table and
  tidy-evaluation namespace requirements.
- Update README examples for the new tail-analysis interfaces while preserving
  the title artwork, version badge, workflow diagram, and table formatting.
- Add package-level regression tests for BED loading, tail mapping, PD/RPP,
  PCA, and pairwise tail-length comparisons.

## Migration

`mapTail()`, `tail_pca()`, and `polyAlength()` are no longer exported by the
supplied source. Use `Map.Tail()`, `Tail.PCA()`, and `Tail.DiffPair()` instead.
PCA arguments and return values have changed. Pairwise analysis now selects
one treatment and one test method per call. See the README and function help.

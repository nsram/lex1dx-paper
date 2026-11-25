# LEX Scala reference

This is a Scala 3 implementation of the same algorithm as the C++ core: atoms, Collapse, Refine with
the R-aligned tie rule, weights, and the same three methods: `Lex.lexHybrid` (the default, coarse
exact start), `Lex.lexPlain` (Collapse start; the deck preset is `Lex.lexPlain(xs, k, Options.deck)`)
and `Lex.ckmeansExact` (the DP).

```sh
python ../fixtures/gen_fixtures.py    # optional: checks the committed shared fixture data
sbt test
```

`CrossCheckTest` reads `../fixtures/expected_cpp.txt`, which the C++ `lex_cli` writes. It requires
every case, weighted ones included, to be bit-identical to C++: cuts, SSE and centre bits, start, atom, sweep and
move counts. Set `LEX_FIXTURES` to use another fixtures directory, or `LEX_SIM` to run the optional
2024-deck reproduction test. Examples: `sbt "runMain examples.Quickstart"` and `sbt "runMain examples.Twelve"`.

Intended differences from C++:

- No presorted flag.
- Fields are named `counts` / `means` (C++ uses `sizes` / `centers`).
- Error messages print doubles the Java way.

Presorted input needs no special handling here: the JDK sorts (`Arrays.sort` on doubles for unweighted
data, TimSort on the index array for weighted data) detect sorted runs and finish in linear time. The
C++ core skips the sort explicitly and reports `presorted`; Scala does not report it.

# Building and testing lex1d

All commands run from the repository root unless a step says `cd`. No absolute paths are needed.
Every test uses generated data, and nothing is downloaded apart from build dependencies (pip
packages, sbt and Scala libraries).

## Minimum versions

| Tool | Minimum | Tested on Linux (the dev box) |
|---|---|---|
| CMake | 3.21 | 4.5.0 |
| C++ compiler | C++17 with floating-point `std::to_chars`: GCC 11, Clang 14 with libc++ / Xcode 14.3, MSVC 2019 16.4 (MSVC 2022 recommended, see Windows) | GCC (conda-forge, x86-64) |
| Python | 3.10, with development headers (`Python.h`) for building the extension | 3.12.15 |
| Python packages | numpy ≥ 1.21, pytest ≥ 7 (scikit-build-core ≥ 0.10 and nanobind ≥ 2.0 are fetched by pip) | |
| JDK | 17 | Temurin 17.0.21 |
| sbt | any recent launcher; the build pins sbt 1.10.7 (`scala/project/build.properties`) | launcher 2.0.12, sbt 1.10.7 |
| Scala | 3.3.4 (fetched by sbt) | 3.3.4 |
| R | a C++17 compiler (Rtools on Windows, Xcode command line tools on macOS); packages Rcpp, testthat ≥ 3 (Ckmeans.1d.dp optional) | R 4.5.3 (conda-forge), Rcpp 1.1.2, testthat 3.3.2, Ckmeans.1d.dp 4.3.6 |

## Linux (verified)

These steps were run literally, in this order, in a fresh clone on the dev box (x86-64). The only
difference was that the clone came from the local repository instead of GitHub.

Prerequisites (Debian/Ubuntu names; use your distribution's equivalents):

```sh
sudo apt install build-essential cmake python3-dev python3-venv openjdk-17-jdk
# sbt: see https://www.scala-sbt.org/download (apt repository, SDKMAN: `sdk install sbt`, or coursier: `cs install sbt`)
# If your distribution's CMake is older than 3.21: python3 -m pip install --user cmake
```

Build and test:

```sh
git clone https://github.com/nsram/lex1d
cd lex1d

# 1. Fixtures: checks the committed fixtures/data against fixtures/SHA256SUMS
python3 fixtures/gen_fixtures.py

# 2. C++: configure, build, test
cmake -S cpp -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure

# 3. Python: install from the checkout into a virtual environment, then test
python3 -m venv .venv
. .venv/bin/activate
pip install ./python pytest
pytest python/tests
#   (the tests build lex_cli with CMake for the bit-identity check; set LEX_CLI=build/lex_cli to reuse step 2's)

# 4. Scala
cd scala
sbt test
cd ..

# 5. R: install from the checkout (configure copies cpp/ into r/src), then test
Rscript -e 'install.packages(c("Rcpp", "testthat", "Ckmeans.1d.dp"))'
R CMD INSTALL r
(cd r/tests/testthat && Rscript -e 'testthat::test_dir(".", package = "lex1d", load_package = "installed")')

# 6. Examples (optional)
python examples/quickstart.py
python examples/embeddings.py
cmake -S examples/cpp -B build-examples -DCMAKE_BUILD_TYPE=Release
cmake --build build-examples -j 4
build-examples/quickstart
(cd scala && sbt "runMain examples.Quickstart")
deactivate
```

Expected results: `python3 fixtures/gen_fixtures.py` prints `ok: 20 data files match SHA256SUMS; 89
cases`; ctest reports 6/6 passed; pytest reports 12 passed; sbt reports Total 32, Passed 31, Skipped 1 (the
optional 2024-deck reproduction, which runs only when `LEX_SIM` points at that data). The Scala
cross-check prints `CROSSCHECK|scala vs c++|356 cases identical|64 weighted|16 with k clamped`.
The R tests report `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 2033 ]` (2033 expectations, including the
356 fixture runs checked bit for bit against `expected_cpp.txt`).

## macOS (not yet run on a Mac)

Prerequisites (Homebrew):

```sh
xcode-select --install          # Apple clang (C++17, floating-point to_chars needs Xcode 14.3+)
brew install cmake python@3.12 openjdk@17 sbt
export JAVA_HOME="$(brew --prefix openjdk@17)/libexec/openjdk.jdk/Contents/Home"
```

Then run the same steps as Linux, starting from `git clone`. Use `python3.12` in place of `python3`
for the venv if your default `python3` is older than 3.10. `-ffp-contract=off` is required on
Apple Silicon (Apple clang otherwise fuses multiply-adds), and CMake adds it automatically.

## Windows (untested until CI runs on GitHub)

Prerequisites:
- **Visual Studio 2022 Build Tools** with the "Desktop development with C++" workload (MSVC 19.30+).
  CMake compiles with `/fp:precise`, which does not contract multiply-adds from VS 2022 on. With
  older MSVC it falls back to `/fp:strict`.
- **CMake 3.21+** (bundled with Visual Studio, or from cmake.org).
- **Python 3.10+** from python.org (includes the headers).
- **JDK 17** (for example Temurin) and **sbt** (MSI from scala-sbt.org).
- Git configured to keep LF line endings in this repository. The repository's `.gitattributes` does
  this, and the fixtures are byte-exact.

In "x64 Native Tools Command Prompt for VS 2022" or PowerShell, from the repository root:

```bat
git clone https://github.com/nsram/lex1d
cd lex1d
python fixtures\gen_fixtures.py
cmake -S cpp -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
python -m venv .venv
.venv\Scripts\activate
pip install .\python pytest
pytest python\tests
cd scala
sbt test
cd ..
```

Visual Studio generators are multi-config, so binaries land in `build\Release\` (for example
`build\Release\lex_cli.exe`), and `ctest` needs `-C Release`.

## R source tarball and R CMD check

`R CMD build` does not run `configure`, so copy the C++ core into `r/src` first:

```sh
(cd r && sh tools/sync_core.sh)
R CMD build r
LEX1D_FIXTURES="$PWD/fixtures" R CMD check --as-cran lex1d_0.1.0.tar.gz
```

On the dev box this gives one package-specific WARNING, `Non-portable flags in variable
'PKG_CXXFLAGS': -ffp-contract=off`. The flag keeps results bit-identical to the C++ library on
platforms with fused multiply-add (Apple Silicon, aarch64). The other findings there came from the
box (no pdflatex, `checkbashisms` or HTML Tidy; conda's `-march=nocona`) and from the package not
being on CRAN yet (new submission, repository URL not public). CI runs the same check on Linux,
macOS and Windows with `error-on: "error"`.

On Windows, run the same commands from a shell that has `sh` (Git Bash, or Rtools' shell).

## Refreshing the C++ reference output

After a change to the C++ core that is meant to change results:

```sh
python3 fixtures/make_expected.py build/lex_cli      # Windows: build\Release\lex_cli.exe
```

The Python test `test_expected_file_is_current` fails until this has been run. The Scala
cross-check reads the refreshed file.

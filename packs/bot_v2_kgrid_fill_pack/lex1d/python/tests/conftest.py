import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
FIXTURES = REPO / "fixtures"


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    """The committed shared fixture set, checked against SHA256SUMS."""
    import hashlib
    for line in (FIXTURES / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split()
        got = hashlib.sha256((FIXTURES / name).read_bytes()).hexdigest()
        assert got == digest, f"fixture {name} differs from SHA256SUMS (fixtures/data is committed; restore it with git)"
    return FIXTURES


@pytest.fixture(scope="session")
def lex_cli(tmp_path_factory) -> Path:
    """Path to the C++ lex_cli: $LEX_CLI if set, else built from ../cpp with CMake."""
    if os.environ.get("LEX_CLI"):
        return Path(os.environ["LEX_CLI"])
    if not shutil.which("cmake"):
        pytest.skip("cmake not found and LEX_CLI not set")
    build = tmp_path_factory.mktemp("cli_build")
    subprocess.run(["cmake", "-S", str(REPO / "cpp"), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release",
                    "-DLEX_BUILD_TESTS=OFF", "-DLEX_BUILD_CLI=ON"], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["cmake", "--build", str(build), "--target", "lex_cli", "--config", "Release", "-j", "4"],
                   check=True, stdout=subprocess.DEVNULL)
    # build/lex_cli (single-config generators) or build/Release/lex_cli.exe (Visual Studio)
    return next(p for p in build.rglob("*") if p.name in ("lex_cli", "lex_cli.exe") and p.is_file())

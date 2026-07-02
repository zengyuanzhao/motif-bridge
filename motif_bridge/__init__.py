from pathlib import Path

from .core import Motif
from .io import (
    read_homer,
    read_json,
    read_meme,
    write_homer,
    write_json,
    write_meme,
)

__all__ = [
    "Motif",
    "read_homer",
    "read_json",
    "read_meme",
    "write_homer",
    "write_json",
    "write_meme",
    "__version__",
]

try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:  # pragma: no cover - importlib.metadata is available on Python 3.8+
    PackageNotFoundError = Exception  # type: ignore[assignment]
    version = None  # type: ignore[assignment]


def _source_tree_version() -> str:
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    try:
        for line in pyproject.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.startswith("version = "):
                return stripped.split("=", 1)[1].strip().strip('"')
    except OSError:
        pass
    return ""


__version__ = _source_tree_version()
if not __version__:
    try:
        if version is None:
            raise PackageNotFoundError
        __version__ = version("motif-bridge")
    except PackageNotFoundError:
        __version__ = "0.3.1"

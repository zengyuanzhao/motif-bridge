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

try:
    if version is None:
        raise PackageNotFoundError
    __version__ = version("motif-bridge")
except PackageNotFoundError:
    __version__ = "0.2.0"

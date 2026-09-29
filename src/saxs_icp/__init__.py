"""ICP-SAXS: rigid-body fitting of atomic structures into SAXS envelopes."""
__version__ = "1.0.0"

from .core import align, proper_transform  # noqa: E402,F401
from .io import read_coords  # noqa: E402,F401

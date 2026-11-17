# ====================================================================
# This file is part of FlexibleSUSY.
#
# FlexibleSUSY is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# FlexibleSUSY is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FlexibleSUSY.  If not, see
# <http://www.gnu.org/licenses/>.
# ====================================================================

"""Helpers for arrays: numpy arrays if numpy is installed, nested lists otherwise."""

from __future__ import annotations

import itertools
import numbers
from typing import Any, Dict, List, Sequence, Tuple

try:  # numpy is optional
    import numpy as _np
except ImportError:  # pragma: no cover - depends on the environment
    _np = None

HAVE_NUMPY = _np is not None


def _nested(shape: Sequence[int], fill: Any) -> Any:
    if not shape:
        return fill
    return [_nested(shape[1:], fill) for _ in range(shape[0])]


def from_entries(shape: Sequence[int], entries: Dict[Tuple[int, ...], Any],
                 fill: Any = 0.0) -> Any:
    """Build an array of the given shape from 0-based index -> value entries."""
    if HAVE_NUMPY:
        is_complex = any(isinstance(v, complex) for v in entries.values())
        result = _np.full(tuple(shape), fill, dtype=complex if is_complex else float)
        for idx, value in entries.items():
            result[idx] = value
        return result
    result = _nested(shape, fill)
    for idx, value in entries.items():
        target = result
        for i in idx[:-1]:
            target = target[i]
        target[idx[-1]] = value
    return result


def to_entries(value: Any, shape: Sequence[int], name: str) -> Dict[Tuple[int, ...], Any]:
    """Check that value has the given shape; return 0-based index -> number.

    Accepts numbers (for shape ()), nested sequences and numpy arrays.
    Raises ValueError with a readable message otherwise.
    """
    shape = tuple(shape)
    if HAVE_NUMPY and isinstance(value, _np.ndarray):
        if tuple(value.shape) != shape:
            raise ValueError(f"{name}: expected shape {shape}, got {tuple(value.shape)}")
        return {idx: value[idx].item() for idx in itertools.product(*(range(n) for n in shape))}
    if not shape:
        if isinstance(value, (list, tuple)) or not isinstance(value, numbers.Number):
            raise ValueError(f"{name}: expected a number, got {value!r}")
        return {(): value}
    entries = {}

    def collect(sub: Any, depth: int, prefix: Tuple[int, ...]) -> None:
        if depth == len(shape):
            if isinstance(sub, (list, tuple)) or not isinstance(sub, numbers.Number):
                raise ValueError(f"{name}: expected shape {shape} with numbers, got {value!r}")
            entries[prefix] = sub
            return
        if not isinstance(sub, (list, tuple)) or len(sub) != shape[depth]:
            raise ValueError(f"{name}: expected shape {shape}, got {value!r}")
        for i, s in enumerate(sub):
            collect(s, depth + 1, prefix + (i,))

    collect(value, 0, ())
    return entries


def to_list(value: Any) -> Any:
    """Convert arrays to nested lists for JSON; complex numbers become [re, im]."""
    if HAVE_NUMPY and isinstance(value, (_np.ndarray, _np.generic)):
        value = value.tolist()
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, list):
        return [to_list(v) for v in value]
    return value

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

"""Minimal reader and writer for the SLHA format as written by FlexibleSUSY.

Blocks map integer keys (tuples for matrix entries) to values.  Values are
floats where possible; entries that are not numbers (e.g. the problem and
warning strings in SPINFO) are kept as strings.  Entries with the same key
are kept in order of appearance, so SPINFO can hold several warnings.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple, Union

Key = Union[int, Tuple[int, ...]]
Value = Union[float, str]


@dataclass
class Block:
    """An SLHA block: name, optional scale Q and entries in input order."""

    name: str
    scale: Optional[float] = None
    entries: List[Tuple[Key, Value]] = field(default_factory=list)
    comments: List[str] = field(default_factory=list)

    def __getitem__(self, key: Key) -> Value:
        for k, v in self.entries:
            if k == key:
                return v
        raise KeyError(f"no entry {key!r} in block {self.name}")

    def __contains__(self, key: Key) -> bool:
        return any(k == key for k, _ in self.entries)

    def get(self, key: Key, default: Optional[Value] = None) -> Optional[Value]:
        try:
            return self[key]
        except KeyError:
            return default

    def get_all(self, key: Key) -> List[Value]:
        """All values with the given key (e.g. all warnings in SPINFO 3)."""
        return [v for k, v in self.entries if k == key]

    def keys(self) -> List[Key]:
        return [k for k, _ in self.entries]

    def __iter__(self) -> Iterator[Tuple[Key, Value]]:
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)


@dataclass
class Channel:
    """A decay channel: branching ratio and PDG codes of the final state."""

    br: float
    pdg: Tuple[int, ...]
    comment: str = ""


@dataclass
class Decay:
    """An SLHA DECAY block."""

    pdg: int
    width: float
    channels: List[Channel] = field(default_factory=list)
    comment: str = ""


@dataclass
class SLHA:
    """Contents of an SLHA file: blocks by upper-case name, decays by PDG code."""

    blocks: Dict[str, Block] = field(default_factory=dict)
    decays: Dict[int, Decay] = field(default_factory=dict)

    def block(self, name: str) -> Optional[Block]:
        return self.blocks.get(name.upper())

    def __contains__(self, name: str) -> bool:
        return name.upper() in self.blocks


def _split_comment(line: str) -> Tuple[str, str]:
    data, _, comment = line.partition("#")
    return data.strip(), comment.strip()


def _to_value(token: str) -> Value:
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return token


def _is_int(token: str) -> bool:
    try:
        int(token)
        return True
    except ValueError:
        return False


def _parse_block_header(data: str) -> Block:
    tokens = data.split()
    name = tokens[1].upper()
    scale = None
    rest = " ".join(tokens[2:])
    if rest:
        rest = rest.replace(" ", "")
        if rest.upper().startswith("Q="):
            scale = float(rest[2:])
    return Block(name=name, scale=scale)


# blocks whose entries are a single integer key followed by free text
_TEXT_BLOCKS = {"SPINFO", "DCINFO"}


def _parse_entry(data: str, comment: str, block: Block) -> None:
    tokens = data.split()
    if block.name in _TEXT_BLOCKS and len(tokens) > 1 and _is_int(tokens[0]):
        block.entries.append((int(tokens[0]), data.split(None, 1)[1]))
        block.comments.append(comment)
        return
    # leading integer tokens form the key; the rest is the value
    n_keys = 0
    while n_keys < len(tokens) - 1 and _is_int(tokens[n_keys]):
        n_keys += 1
    if n_keys == 0:
        # entry without key (e.g. Block ALPHA)
        block.entries.append(((), _to_value(" ".join(tokens))))
    else:
        ints = tuple(int(t) for t in tokens[:n_keys])
        key: Key = ints[0] if n_keys == 1 else ints
        value_str = " ".join(tokens[n_keys:])
        block.entries.append((key, _to_value(value_str)))
    block.comments.append(comment)


def read(text: str) -> SLHA:
    """Parse SLHA text."""
    result = SLHA()
    block: Optional[Block] = None
    decay: Optional[Decay] = None
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        data, comment = _split_comment(raw)
        head = data.split()[0].upper()
        if head == "BLOCK":
            block = _parse_block_header(data)
            decay = None
            result.blocks[block.name] = block
        elif head == "DECAY":
            tokens = data.split()
            decay = Decay(pdg=int(tokens[1]), width=float(tokens[2]), comment=comment)
            block = None
            result.decays[decay.pdg] = decay
        elif decay is not None:
            tokens = data.split()
            n = int(tokens[1])
            decay.channels.append(
                Channel(br=float(tokens[0]),
                        pdg=tuple(int(t) for t in tokens[2:2 + n]),
                        comment=comment))
        elif block is not None:
            _parse_entry(data, comment, block)
    return result


def _format_value(value: Value) -> str:
    if isinstance(value, str):
        return value
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"cannot write non-finite value {v} to SLHA")
    return f"{v: .17e}"


def format_block(name: str, entries: Dict[Key, Value]) -> str:
    """Format a block with the given entries as SLHA text."""
    lines = [f"Block {name}"]
    for key, value in entries.items():
        keys = key if isinstance(key, tuple) else (key,)
        lines.append("   " + "   ".join(str(k) for k in keys) + "   " + _format_value(value))
    return "\n".join(lines) + "\n"


def write(blocks: Dict[str, Dict[Key, Value]]) -> str:
    """Format several blocks as SLHA text, in the given order."""
    return "".join(format_block(name, entries) for name, entries in blocks.items())

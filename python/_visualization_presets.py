"""Display presets: material colours, volume layers, saved views and sections.

A preset is a mapping or a JSON file supplied by the caller. It changes only how
a model is shown. Lengths are metres in the geometry frame.
"""
import fnmatch
import json
import os
from pathlib import Path
import re

import numpy as np

BACKGROUND = (.914, .929, .941)
PANEL = (.969, .976, .980)
INK = (.110, .165, .200)
SOFT = (.333, .400, .439)
RULE = (.788, .827, .851)
ACCENT = (.059, .435, .549)


def clean_name(name):
    """Drop the pointer suffix Geant4 writes into exported names."""
    return re.sub(r"0x[0-9a-fA-F]+$", "", str(name))


def rgb(colour):
    """(r, g, b) in 0..1 from "#rrggbb" or three numbers."""
    if isinstance(colour, str):
        if not re.fullmatch(r"#[0-9a-fA-F]{6}", colour):
            raise ValueError("colour must be '#rrggbb' or three numbers: %r" % colour)
        return tuple(int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    out = tuple(float(c) for c in colour)
    if len(out) != 3 or not all(0 <= c <= 1 for c in out):
        raise ValueError("colour must have three components between 0 and 1")
    return out


def material_row(name, rows):
    """The caller's first row matching a material, or None."""
    plain = clean_name(name)
    return next((row for row in rows
                 if any(fnmatch.fnmatchcase(plain, glob) for glob in row["match"])), None)


def _vector(value, name, size=3):
    v = np.asarray(value, dtype=float)
    if v.shape != (size,) or not np.isfinite(v).all():
        raise ValueError("%s must be %d finite numbers" % (name, size))
    return [float(c) for c in v]


def _globs(value):
    return [value] if isinstance(value, str) else [str(g) for g in value]


def _known(row, allowed, where):
    extra = set(row) - set(allowed)
    if extra:
        raise ValueError("unknown %s entries: %s" % (where, ", ".join(sorted(extra))))


def normalize_section(row):
    """Validate one section; *row* needs an origin and a normal."""
    _known(row, ("label", "origin", "normal", "up", "centre", "width", "height"), "section")
    if "origin" not in row or "normal" not in row:
        raise ValueError("a section needs 'origin' and 'normal'")
    out = dict(label=str(row.get("label", "Section")), origin=_vector(row["origin"], "origin"),
               normal=_vector(row["normal"], "normal"))
    if not np.any(out["normal"]):
        raise ValueError("section normal must be a nonzero vector")
    for key in ("up", "centre"):
        out[key] = _vector(row[key], key) if row.get(key) is not None else None
    for key in ("width", "height"):
        out[key] = float(row[key]) if row.get(key) is not None else None
    return out


def load_presets(presets):
    """Validated presets from a mapping, a JSON path or None."""
    if presets is None:
        presets = {}
    elif isinstance(presets, (str, os.PathLike)):
        presets = json.loads(Path(presets).read_text())
    _known(presets, ("title", "note", "up", "home", "opacity", "materials", "layers", "sections",
                     "hidden_volumes"),
           "preset")
    out = dict(title=str(presets.get("title", "")), note=str(presets.get("note", "")),
               up=_vector(presets.get("up", (0., 0., 1.)), "up"), home=None,
               hidden_volumes=_globs(presets.get("hidden_volumes", ())))
    if not np.any(out["up"]):
        raise ValueError("preset 'up' must be a nonzero vector")
    out["opacity"] = float(presets.get("opacity", .9))
    if not 0 < out["opacity"] <= 1:
        raise ValueError("preset opacity must be above 0 and at most 1")
    home = presets.get("home")
    if home is not None:
        _known(home, ("centre", "direction", "width", "height"), "home")
        out["home"] = dict(centre=_vector(home["centre"], "home centre"),
                           direction=_vector(home["direction"], "home direction"),
                           width=float(home["width"]), height=float(home["height"]))
    out["materials"] = []
    for row in presets.get("materials", ()):
        _known(row, ("label", "match", "colour"), "material")
        out["materials"].append(dict(label=str(row["label"]), match=_globs(row["match"]),
                                     colour=rgb(row["colour"])))
    out["layers"] = []
    for row in presets.get("layers", ()):
        _known(row, ("label", "match", "shown", "edges"), "layer")
        out["layers"].append(dict(label=str(row["label"]), match=_globs(row["match"]),
                                  shown=bool(row.get("shown", True)),
                                  edges=bool(row.get("edges", True))))
    out["sections"] = [normalize_section(row) for row in presets.get("sections", ())]
    for key in ("layers", "sections"):
        labels = [row["label"] for row in out[key]]
        if len(set(labels)) != len(labels):
            raise ValueError("preset %s need different labels" % key)
    return out


OTHER = "Everything else"


def assign_layers(scene, layers):
    """Give every instance a layer label: its own match, else its mother's."""
    cache = {}
    def match(name):
        if name not in cache:
            cache[name] = next((row["label"] for row in layers for glob in row["match"]
                                if fnmatch.fnmatchcase(name, glob)), None)
        return cache[name]
    assigned = []
    for row in scene["instances"]:
        parent = row.get("parent", -1)
        label = (match(clean_name(row["name"].rsplit("/", 1)[-1])) or
                 match(clean_name(scene["prototypes"][row["prototype"]]["name"])))
        # A layer named for the world itself would claim everything else.
        assigned.append(label or (assigned[parent] if 0 < parent < len(assigned) else OTHER))
    return assigned

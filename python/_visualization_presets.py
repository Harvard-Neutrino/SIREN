"""Display presets: material colors, volume layers, saved views and sections.

A preset is a mapping or a JSON file supplied by the caller. It changes only how
a model is shown. Lengths are meters in the geometry frame.
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


def rgb(color):
    """(r, g, b) in 0..1 from "#rrggbb" or three numbers."""
    if isinstance(color, str):
        if not re.fullmatch(r"#[0-9a-fA-F]{6}", color):
            raise ValueError("color must be '#rrggbb' or three numbers: %r" % color)
        return tuple(int(color[i:i + 2], 16) / 255 for i in (1, 3, 5))
    out = tuple(float(c) for c in color)
    if len(out) != 3 or not all(0 <= c <= 1 for c in out):
        raise ValueError("color must have three components between 0 and 1")
    return out


def material_row(name, rows):
    """The caller's first row matching a material, or None."""
    plain = clean_name(name)
    return next((row for row in rows
                 if any(fnmatch.fnmatchcase(plain, glob) for glob in row["match"])), None)


SURFACES = ("solid", "see-through")


def _choice(value, allowed, name):
    if value not in allowed:
        raise ValueError("%s must be one of %s, not %r" % (name, ", ".join(allowed), value))
    return value


def _look(row):
    """Color, solid or see-through, and shown: each optional in a row."""
    return dict(color=rgb(row["color"]) if row.get("color") is not None else None,
                opacity=(_choice(row["opacity"], SURFACES, "opacity")
                         if row.get("opacity") is not None else None),
                shown=bool(row.get("shown", True)))


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
    _known(row, ("label", "origin", "normal", "up", "center", "width", "height"), "section")
    if "origin" not in row or "normal" not in row:
        raise ValueError("a section needs 'origin' and 'normal'")
    out = dict(label=str(row.get("label", "Section")), origin=_vector(row["origin"], "origin"),
               normal=_vector(row["normal"], "normal"))
    if not np.any(out["normal"]):
        raise ValueError("section normal must be a nonzero vector")
    for key in ("up", "center"):
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
    _known(presets, ("title", "note", "up", "home", "opacity", "see_through", "color_mode",
                     "surfaces", "materials", "layers", "volumes", "sections", "hidden_volumes"),
           "preset")
    out = dict(title=str(presets.get("title", "")), note=str(presets.get("note", "")),
               up=_vector(presets.get("up", (0., 0., 1.)), "up"), home=None,
               hidden_volumes=_globs(presets.get("hidden_volumes", ())))
    if not np.any(out["up"]):
        raise ValueError("preset 'up' must be a nonzero vector")
    for key, default in (("opacity", .9), ("see_through", .25)):
        out[key] = float(presets.get(key, default))
        if not 0 < out[key] <= 1:
            raise ValueError("preset %s must be above 0 and at most 1" % key)
    out["color_mode"] = _choice(presets.get("color_mode", "custom"), ("custom", "density"), "color_mode")
    out["surfaces"] = _choice(presets.get("surfaces", "solid"), SURFACES, "surfaces")
    home = presets.get("home")
    if home is not None:
        _known(home, ("center", "direction", "width", "height"), "home")
        out["home"] = dict(center=_vector(home["center"], "home center"),
                           direction=_vector(home["direction"], "home direction"),
                           width=float(home["width"]), height=float(home["height"]))
    out["materials"] = []
    for row in presets.get("materials", ()):
        _known(row, ("label", "match", "color", "opacity", "shown"), "material")
        out["materials"].append(dict(label=str(row["label"]), match=_globs(row["match"]), **_look(row)))
    out["layers"] = []
    for row in presets.get("layers", ()):
        _known(row, ("label", "match", "shown", "edges", "color", "opacity"), "layer")
        out["layers"].append(dict(label=str(row["label"]), match=_globs(row["match"]),
                                  edges=bool(row.get("edges", True)), **_look(row)))
    out["volumes"] = []
    for row in presets.get("volumes", ()):
        _known(row, ("label", "match", "color", "opacity", "shown"), "volume")
        globs = _globs(row["match"])
        out["volumes"].append(dict(label=str(row.get("label", globs[0])), match=globs, **_look(row)))
    out["sections"] = [normalize_section(row) for row in presets.get("sections", ())]
    for key in ("materials", "layers", "volumes", "sections"):
        labels = [row["label"] for row in out[key]]
        if len(set(labels)) != len(labels):
            raise ValueError("preset %s need different labels" % key)
    return out


def to_json(presets):
    """A normalized preset as JSON-ready data, without the defaults."""
    def color(c):
        return "#%02x%02x%02x" % tuple(int(round(255 * x)) for x in c)

    def looks(row, *keys):
        out = {k: row[k] for k in ("label", "match") + keys}
        if row.get("color") is not None:
            out["color"] = color(row["color"])
        if row.get("opacity") is not None:
            out["opacity"] = row["opacity"]
        if not row.get("shown", True):
            out["shown"] = False
        return out
    out = {k: presets[k] for k in ("title", "note", "up", "opacity", "see_through", "color_mode",
                                   "surfaces", "hidden_volumes")}
    if presets["home"] is not None:
        out["home"] = presets["home"]
    out["materials"] = [looks(row) for row in presets["materials"]]
    out["layers"] = [dict(looks(row), shown=row["shown"], **({} if row["edges"] else {"edges": False}))
                     for row in presets["layers"]]
    out["volumes"] = [looks(row) for row in presets["volumes"]]
    out["sections"] = [{k: v for k, v in row.items() if v is not None} for row in presets["sections"]]
    return {k: v for k, v in out.items() if v not in ([], "")}


def assign_volumes(scene, volumes):
    """The first volume rule naming each instance (placed or logical name), or None."""
    def match(name):
        return next((row["label"] for row in volumes for glob in row["match"]
                     if fnmatch.fnmatchcase(name, glob)), None)
    cache = {}
    out = []
    for row in scene["instances"]:
        key = (row["name"].rsplit("/", 1)[-1], row["prototype"])
        if key not in cache:
            cache[key] = (match(clean_name(key[0])) or
                          match(clean_name(scene["prototypes"][row["prototype"]]["name"])))
        out.append(cache[key])
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

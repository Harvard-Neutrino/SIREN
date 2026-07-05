"""Command-line converter: SIREN ``.siren_events`` -> HepMC3 ascii.

Run as::

    python -m siren.hepmc3_convert <input>.siren_events [-o <output>.hepmc3]

The compiled extension module is ``siren.hepmc3``; its name cannot be reused by a
pure-Python module, so the converter CLI lives here as ``siren.hepmc3_convert``.
"""


def main():
    import argparse
    from . import _util
    p = argparse.ArgumentParser(
        prog="python -m siren.hepmc3_convert",
        description="Convert a SIREN .siren_events file to a HepMC3 ascii file.")
    p.add_argument("input", help="path to a .siren_events file")
    p.add_argument("-o", "--output", default=None,
                   help="output .hepmc3 path (default: <input>.hepmc3)")
    args = p.parse_args()
    out = _util.convert_siren_events_to_hepmc3(args.input, args.output)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Top-level ``siren`` namespace size guard.

This lives with the full assembled package because it
checks the complete top-level export set, not the Simulation-API subset.

The count is measured in a fresh subprocess so it reflects the *intended*
public API and is not inflated by submodules that other tests in the same
session have imported into the ``siren`` namespace (which makes ``dir(siren)``
grow and the check order-dependent).
"""

import subprocess
import sys


def test_top_level_count_reasonable():
    """The top-level siren namespace should not be bloated."""
    code = (
        "import siren; "
        "print(len([x for x in dir(siren) if not x.startswith('_')]))"
    )
    out = subprocess.check_output([sys.executable, "-c", code], text=True)
    count = int(out.strip())
    # The spec vocabulary (Vertex, Directed, channels, expand, errors,
    # generate, Propagated, Fixed) is part of the intended public surface.
    assert count < 45, f"Too many top-level names: {count}"

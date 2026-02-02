#!/usr/bin/env python3
from pathlib import Path

from snakebids import bidsapp, plugins

# For Sphinx doc generation, __file__ may not be defined
try:
    app_path = Path(__file__).resolve().parent
except NameError:
    # Fallback for Sphinx - use current working directory or find snakedwi module
    import snakedwi
    app_path = Path(snakedwi.__file__).parent

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(app_path),
        plugins.Version(distribution="snakedwi"),
    ]
)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()

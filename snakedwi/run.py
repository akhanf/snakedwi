#!/usr/bin/env python3
import os

from snakebids.app import SnakeBidsApp


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    app = SnakeBidsApp("../snakedwi", skip_parse_args=True)
    return app.parser


def main():
    app = SnakeBidsApp(
        os.path.abspath(os.path.dirname(__file__)),
        configfile_path="config/snakebids.yml",
    )
    app.run_snakemake()


if __name__ == "__main__":
    main()

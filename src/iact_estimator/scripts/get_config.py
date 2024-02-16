"""Command line script to get the example configuration file."""

import argparse
from pathlib import Path
import shutil

from .. import RESOURCES_PATH

parser = argparse.ArgumentParser()

parser.add_argument(
    "--output-path",
    default=None,
    type=str,
    help="Where to save the configuration file (default: current working directory).",
)


def main():
    args = parser.parse_args()

    output_path = Path.cwd() if args.output_path is None else Path(args.output_path)

    shutil.copy(RESOURCES_PATH / "config.yml", output_path)


if __name__ == "__main__":
    main()

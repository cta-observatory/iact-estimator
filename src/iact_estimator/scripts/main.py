"""Script to estimate telescopes performance."""

import argparse
import logging
from pathlib import Path
import shutil

from astroplan import FixedTarget, Observer
from astropy.time import Time
import astropy.units as u
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt

from .. import __version__
from ..io import read_yaml
from ..core import (
    setup_logging,
    initialize_model,
    check_input_configuration,
    prepare_data,
    source_detection,
    calculate,
)
from ..plots import (
    plot_spectrum,
    plot_sed,
    plot_transit,
    plot_altitude_airmass,
    plot_observability_constraints_grid,
)
from ..observability import define_constraints, check_observability
from .. import RESOURCES_PATH

parser = argparse.ArgumentParser()

parser.add_argument("--version", action="store_true")

subparsers = parser.add_subparsers(dest="command")

parser_config = subparsers.add_parser(
    "config",
    help="Copy the default configuration file to your new project directory.",
)

parser_config.add_argument(
    "--to",
    default=None,
    type=str,
    help="Where to save the configuration file (default: current working directory).",
)

parser_run = subparsers.add_parser(
    "run",
    help="Launch the estimation of observability and physical properties of the target source.",
)

parser_run.add_argument(
    "--config", required=True, type=str, help="Path to configuration file."
)
parser_run.add_argument(
    "--source-name",
    default="test_source",
    type=str,
    help="Name of the source to estimate.",
)
parser_run.add_argument(
    "--output-path",
    default=None,
    type=str,
    help="Path where output will be saved (defaults to current working directory)",
)
parser_run.add_argument(
    "--log-level",
    default="INFO",
    type=str.upper,
    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    help="Logging level.",
)
parser_run.add_argument(
    "--overwrite", action="store_true", help="Overwrite any generated output."
)


def main():
    global_args = parser.parse_args()

    if global_args.version:
        parser.exit(0, message=f"{__version__}\n")

    namespace = parser.parse_args()

    if namespace.command == "config":
        args_config = parser.parse_args()
        config_output_path = (
            Path.cwd() if args_config.to is None else Path(args_config.to)
        )
        shutil.copy(RESOURCES_PATH / "config.yml", config_output_path)

        parser.exit(0, message="You're halfway there! Have fun!\n")

    if namespace.command == "run":
        args = parser.parse_args()

        output_path = (
            Path(args.output_path) if args.output_path is not None else Path.cwd()
        )
        source_name = args.source_name

        previous_output = len(
            [file for file in output_path.rglob(f"**/{source_name}*")]
        )
        if (previous_output > 0) and (not args.overwrite):
            raise ValueError(
                "Previous results are present and --overwrite option was not used."
            )

        logger = setup_logging(args.log_level, source_name)

        logger.info("Loading configuration file")
        config = read_yaml(args.config)

        logging.info("Validating input configuration")
        if not check_input_configuration(config):
            logging.critical(
                "One or more invalid configuration settings have been found."
            )
            parser.exit(status=1, message="iact-estimator terminated abnormally.")

        plotting_options = config["plotting_options"]
        use_seaborn = config["use_seaborn"]
        if use_seaborn:
            import seaborn as sns

            seaborn_options = config["seaborn_options"]
            sns.set_theme(**seaborn_options)

        # Basic observability checks
        target_source = FixedTarget.from_name(source_name)
        observer = Observer.at_site("Roque de los Muchachos")

        crab = FixedTarget.from_name("Crab")

        constraints = define_constraints(config)

        from datetime import datetime

        start_datetime = (
            Time(config["observation"]["start_datetime"])
            if config["observation"]["start_datetime"] is not None
            else Time(datetime.now(tz=observer.timezone))
        )
        end_datetime = (
            Time(config["observation"]["end_datetime"])
            if config["observation"]["end_datetime"] is not None
            else start_datetime + 1 * u.day
        )
        logger.info("Observation starts at %s", start_datetime)
        logger.info("Observation ends at %s", end_datetime)

        time_range = [start_datetime, end_datetime]
        time_grid_resolution = (
            u.Quantity(config["observation"]["time_resolution"])
            if config["observation"]["time_resolution"]
            else 1 * u.h
        )

        ever_observable, best_months = check_observability(
            constraints, observer, [target_source], time_range, time_grid_resolution
        )

        _ = plot_observability_constraints_grid(
            source_name,
            config,
            observer,
            target_source,
            start_datetime,
            end_datetime,
            time_grid_resolution,
            constraints,
            ax=None,
            savefig=True,
            output_path=output_path,
        )

        if not ever_observable:
            logger.info("The source is never observable from this location!")
            parser.exit(0)

        logger.info(f"The best months to observe the target source are {best_months}.")

        with quantity_support():
            plot_transit(
                config,
                source_name,
                target_source,
                observer,
                start_datetime,
                merge_profiles=config["plotting_options"]["merge_horizon_profiles"],
                plot_crab=True if (crab.coord == target_source.coord) else False,
                style_kwargs=None,
                savefig=True,
                output_path=output_path,
            )

            plot_altitude_airmass(
                config,
                source_name,
                target_source,
                observer,
                start_datetime,
                brightness_shading=True,
                airmass_yaxis=True,
                savefig=True,
                output_path=output_path,
            )

        logger.info("Initializing assumed model")
        assumed_spectrum = initialize_model(config)

        plot_energy_bounds = [
            u.Quantity(plotting_options["min_energy"]),
            u.Quantity(plotting_options["max_energy"]),
        ]

        logger.info("Producing plot of the assumed source model")
        with quantity_support():
            plot_spectrum(
                config,
                plot_energy_bounds,
                assumed_spectrum,
                source_name,
                plotting_options,
                savefig=True,
                output_path=None,
                label=source_name,
            )

        energy_bins, gamma_rate, background_rate = prepare_data(config)

        en, sed, dsed, sigmas, detected = calculate(
            energy_bins, gamma_rate, background_rate, config, assumed_spectrum
        )

        combined_significance = source_detection(
            sigmas, u.Quantity(config["observation"]["time"])
        )

        with quantity_support():
            plot_sed(
                config,
                sigmas,
                combined_significance,
                source_name,
                assumed_spectrum,
                en,
                sed,
                dsed,
                detected,
                savefig=True,
                output_path=output_path,
            )

        logger.info("All expected operations have been perfomed succesfully.")

        if config["plotting_options"]["show"]:
            plt.show()


if __name__ == "__main__":
    main()

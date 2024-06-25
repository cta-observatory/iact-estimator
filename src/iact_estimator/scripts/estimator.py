"""Script to estimate telescopes performance."""

import argparse
import logging
from pathlib import Path

from astroplan import FixedTarget, Observer
from astropy.time import Time
import astropy.units as u
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt

from ..io import read_yaml, load_performance_ecsv
from ..core import (
    setup_logging,
    initialize_model,
    check_input_configuration,
    prepare_data,
    source_detection,
    calculate,
)
from ..plots import plot_spectrum, plot_sed, plot_transit, plot_altitude_airmass

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config", required=True, type=str, help="Path to configuration file."
)
parser.add_argument(
    "--source-name",
    default="test_source",
    type=str,
    help="Name of the source to estimate.",
)
parser.add_argument(
    "--performance",
    default="",
    type=str,
    help="Custom performance data.",
)
parser.add_argument(
    "--output-path",
    default=None,
    type=str,
    help="Path where output will be saved (defaults to current working directory)",
)
parser.add_argument(
    "--log-level",
    default="INFO",
    type=str.upper,
    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    help="Logging level.",
)
parser.add_argument(
    "--overwrite", action="store_true", help="Overwrite any generated output."
)


def main():
    args = parser.parse_args()

    output_path = Path(args.output_path) if args.output_path is not None else Path.cwd()
    source_name = args.source_name

    previous_output = len([file for file in output_path.rglob(f"**/{source_name}*")])
    if (previous_output > 0) and (not args.overwrite):
        raise ValueError(
            "Previous results are present and --overwrite option was not used."
        )

    logger = setup_logging(args.log_level, source_name)

    logger.info("Loading configuration file")
    config = read_yaml(args.config)

    performance_data = None
    if args.performance:
        performance_data_file = Path(args.performance).resolve()
        logger.info("Loading performance data from %s", performance_data_file)
        performance_data = load_performance_ecsv(performance_data_file)

    logging.info("Validating input configuration")
    if not check_input_configuration(config, performance_data):
        logging.critical("One or more invalid configuration settings have been found.")
        parser.exit(status=1, message="iact-estimator terminated abnormally.")

    plotting_options = config["plotting_options"]
    use_seaborn = config["use_seaborn"]
    if use_seaborn:
        import seaborn as sns

        seaborn_options = config["seaborn_options"]
        sns.set_theme(**seaborn_options)

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

    if not performance_data:
        energy_bins, gamma_rate, background_rate = prepare_data(config)
    else:
        energy_bins, gamma_rate, background_rate = prepare_data(performance_data)

    en, sed, dsed, sigmas, detected = calculate(
        energy_bins, gamma_rate, background_rate, config, assumed_spectrum
    )

    combined_significance = source_detection(
        sigmas, u.Quantity(config["observation_time"])
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

    target_source = FixedTarget.from_name(source_name)
    observer = Observer.at_site("Roque de los Muchachos")
    time = Time(config["observation_datetime"])

    crab = FixedTarget.from_name("Crab")

    with quantity_support():
        plot_transit(
            config,
            source_name,
            target_source,
            observer,
            time,
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
            time,
            brightness_shading=True,
            airmass_yaxis=True,
            savefig=True,
            output_path=output_path,
        )

    plt.show()


if __name__ == "__main__":
    main()

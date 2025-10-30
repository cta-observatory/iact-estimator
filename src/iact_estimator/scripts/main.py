"""Script to estimate telescopes performance."""

import argparse
import logging
from pathlib import Path
from datetime import datetime
import shutil

from astroplan import FixedTarget, Observer
from astropy.coordinates.name_resolve import NameResolveError
from astropy.time import Time
import astropy.units as u
from astropy.visualization import quantity_support
import matplotlib.pyplot as plt

from .. import __version__
from ..io import read_yaml, save_fits_hdu, load_performance_ecsv
from ..core import (
    setup_logging,
    load_target_source_coordinates,
    initialize_model,
    check_input_configuration,
    prepare_data,
    source_detection,
    calculate,
)
from ..plots.physics import plot_spectrum, plot_sed
from ..plots.observability import (
    plot_transit,
    plot_altitude_airmass,
    plot_observability_constraints_grid,
    create_observability_heatmap,
)
from ..plots.multi_wavelength import plot_from_skyview_survey
from ..observability import (
    define_constraints,
    check_observability,
    get_days_in_this_year,
)
from ..plots.wobble_skymap import plot_skymap_with_wobbles, load_wobbles
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
parser.add_argument(
    "--performance",
    default="",
    type=str,
    help="Custom performance data.",
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

        logging.info("Loading configuration file")
        config = read_yaml(args.config)

        source_name = (
            config["target_source"]["name"]
            if config["target_source"]["name"]
            else "test_source"
        )

        previous_output = len(
            [file for file in output_path.rglob(f"**/{source_name}*")]
        )
        if (previous_output > 0) and (not args.overwrite):
            raise ValueError(
                "Previous results are present and --overwrite option was not used."
            )

        logger = setup_logging(args.log_level, source_name)

        performance_data = None
        if args.performance:
            performance_data_file = Path(args.performance).resolve()
            logger.info("Loading performance data from %s", performance_data_file)
            performance_data = load_performance_ecsv(performance_data_file)

        logging.info("Validating input configuration")
        if not check_input_configuration(config, performance_data):
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
        if source_name:
            try:
                target_source = FixedTarget.from_name(source_name)
            except NameResolveError:
                target_source = load_target_source_coordinates(config)
        elif (
            source_name != "test_source"
            and config["target_source"]["coordinates"]["force"]
        ):
            target_source = load_target_source_coordinates(config)

        if config["observer"]["auto"]:
            observer = Observer.at_site(config["observer"]["auto"])
        else:
            obs_cfg = config["observer"]["manual"]
            observer = Observer(
                timezone=obs_cfg["timezone"],
                name=obs_cfg["name"],
                latitude=u.Quantity(obs_cfg["latitude"]),
                longitude=u.Quantity(obs_cfg["longitude"]),
                elevation=u.Quantity(obs_cfg["elevation"]),
            )

        logger.debug("Defining observation constraints")
        constraints = define_constraints(config)

        start_datetime = (
            Time(config["observation"]["start_datetime"])
            if config["observation"]["start_datetime"] is not None
            else Time(datetime.now(tz=observer.timezone))
        )
        year_days = get_days_in_this_year()
        end_datetime = (
            Time(config["observation"]["end_datetime"])
            if config["observation"]["end_datetime"] is not None
            else start_datetime + year_days
        )
        logger.info("Observation starts at %s", start_datetime)
        logger.info("Observation ends at %s", end_datetime)

        time_range = [start_datetime, end_datetime]
        time_grid_resolution = (
            u.Quantity(config["observation"]["time_resolution"])
            if config["observation"]["time_resolution"]
            else 1 * u.h
        )

        logger.debug("Checking observability")
        ever_observable, best_months = check_observability(
            constraints, observer, [target_source], time_range, time_grid_resolution
        )

        logger.debug("Producing observability constraints grid")
        obs_grid_time_res = (
            1 * u.h
            if (end_datetime - start_datetime).to("day").value <= 1
            else 1 * u.day
        )
        _ = plot_observability_constraints_grid(
            source_name,
            config,
            observer,
            target_source,
            start_datetime,
            end_datetime,
            obs_grid_time_res,
            constraints,
            ax=None,
            savefig=True,
            output_path=output_path,
        )

        if not ever_observable:
            logger.info("The source is never observable from this location!")
            parser.exit(0)

        logger.info(f"The best months to observe the target source are {best_months}.")

        logger.debug("Producing observability heatmap")
        create_observability_heatmap(
            target_source,
            observer,
            constraints,
            start_datetime,
            end_datetime,
            time_resolution=1 * u.hour,
            cmap="YlGnBu",
            sns_plotting_context="paper",
            sns_axes_style="whitegrid",
            savefig=True,
            output_path=None,
            save_format="png",
        )

        with quantity_support():
            plot_transit(
                config,
                source_name,
                target_source,
                observer,
                start_datetime,
                merge_profiles=config["plotting_options"]["merge_horizon_profiles"],
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

        for survey in config["skyview"]["surveys"]:
            survey_name = survey["name"]
            fig, ax = plt.subplots()
            ax, hdu = plot_from_skyview_survey(
                target_source,
                survey_name=survey_name,
                fov_radius=u.Quantity(survey["fov_radius"]),
                log=survey["log"],
                ax=ax,
                reticle=survey["reticle"],
                style_kwargs=survey["style_kwargs"],
                reticle_style_kwargs=survey["reticle_style_kwargs"],
            )
            ax.set_title(f"{source_name} - {survey_name}")
            output_path = output_path if output_path is not None else Path.cwd()
            fig.savefig(
                output_path
                / f"{source_name}_{survey_name}.{config['plotting_options']['file_format']}",
                bbox_inches=config["plotting_options"]["bbox_inches"],
            )
            if config["skyview"]["save_hdus"]:
                save_fits_hdu(
                    hdu,
                    output_path
                    / f"{source_name}_skyview_image_from_{survey_name.replace(' ', '')}.fits",
                    overwrite=args.overwrite,
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

        if not performance_data:
            energy_bins, gamma_rate, background_rate = prepare_data(config)
        else:
            energy_bins, gamma_rate, background_rate = prepare_data(
                config, performance_data
            )

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

        instrument_fov = u.Quantity(config["fov"])
        wobble_offsets, wobble_angles = load_wobbles(config["wobbles"])
        plot_skymap_with_wobbles(
            target_source,
            observer,
            instrument_fov,
            wobble_angles,
            wobble_offsets,
            config,
        )

        logger.info("All expected operations have been perfomed succesfully.")

        logger.info("All output can be found at %s", output_path)

        if config["plotting_options"]["show"]:
            plt.show()


if __name__ == "__main__":
    main()

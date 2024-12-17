"""Script to estimate telescopes performance."""

import argparse
import logging
from pathlib import Path
import shutil

from astroplan import FixedTarget, Observer
from astropy.time import Time
import astropy.units as u
from astropy.visualization import quantity_support
from gammapy.irf import load_irf_dict_from_file

from gammapy.datasets import SpectrumDataset
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SpectrumDatasetMaker,
    WobbleRegionsFinder,
)
from gammapy.modeling.models import Models, PointSpatialModel
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
from ..gammapy_interface import (
    load_real_observations,
    create_simulated_observation,
    has_radmax_2d,
    define_on_region_geometry,
    estimate_sed,
    fake_onoff_from_fake_observation,
    fake_onoff_from_real_observations,
    load_energy_axis_from_config,
    get_wstat_table,
)
from ..plots import (
    plot_spectrum,
    plot_sed,
    plot_transit,
    plot_altitude_airmass,
    plot_gammapy_sed,
)
from .. import RESOURCES_PATH

parser = argparse.ArgumentParser()

parser.add_argument("--version", action="store_true")

subparsers = parser.add_subparsers(dest="command")

parser_config = subparsers.add_parser(
    "config",
    help="Copy the configuration files to your new project directory.",
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
    "--use-gammapy", action="store_true", help="Use also the new gammapy interface."
)

parser_run.add_argument(
    "--config", required=True, type=str, help="Path to configuration file."
)
parser_run.add_argument(
    "--models",
    required=False,
    type=str,
    help="Path to gammapy model file (models.yml).",
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
        shutil.copy(RESOURCES_PATH / "models.yml", config_output_path)

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

        simulated_observation_livetime = u.Quantity(config["observation_time"])

        combined_significance = source_detection(sigmas, simulated_observation_livetime)

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

        target_source_coordinates = target_source.coord

        if args.use_gammapy:
            gammapy_config = config["gammapy_interface"]

            gammapy_analysis_type = gammapy_config["datasets"]["type"]

            if gammapy_analysis_type == "1d":
                required_irfs = "point-like"
            else:
                raise NotImplementedError(
                    "Only 1D analyses from point-like simulations are available at the moment with the gammapy interface."
                )

            gammapy_input = Path(gammapy_config["input"])

            sky_models = Models.from_dict(read_yaml(args.models))

            if gammapy_input.is_dir():
                logger.info(
                    "Creating data store from real observations at %s", gammapy_input
                )
                observations = load_real_observations(gammapy_input, required_irfs)
            else:
                logger.info(
                    "Loading instrument response functions from %s", gammapy_input
                )
                irfs = load_irf_dict_from_file(gammapy_input)

                observations = create_simulated_observation(
                    target_source_coordinates,
                    observer,
                    simulated_observation_livetime,
                    irfs,
                )

            reco_energy_axis = load_energy_axis_from_config(gammapy_config, "energy")
            true_energy_axis = load_energy_axis_from_config(
                gammapy_config, "energy_true"
            )

            on_region_config = gammapy_config["datasets"]["on_region"]
            offset = on_region_config["offset"]
            on_region_radius = on_region_config["radius"] or None
            on_region_geometry = define_on_region_geometry(
                target_source_coordinates,
                observations[0],
                offset,
                reco_energy_axis,
                on_region_radius=on_region_radius,
            )

            empty_spectrum_dataset = SpectrumDataset.create(
                geom=on_region_geometry,
                energy_axis_true=true_energy_axis,
            )

            logger.info("Angular cut depends on reconstructed energy.")

            is_point_like = has_radmax_2d(observations[0])

            containment_correction = gammapy_config["datasets"][
                "containment_correction"
            ]
            if is_point_like:
                logger.warning(
                    "Input observation is point-like, so containment correction will be disabled."
                )
                containment_correction = False

            n_off_regions = config["n_off_regions"]

            if len(observations) >= 1 and has_radmax_2d(observations[0]):
                region_finder = WobbleRegionsFinder(n_off_regions=n_off_regions)
                background_maker = ReflectedRegionsBackgroundMaker(
                    region_finder=region_finder
                )

                spectrum_dataset_maker = SpectrumDatasetMaker(
                    containment_correction=containment_correction,
                    selection=["counts", "exposure", "edisp"],
                )

                empty_spectrum_dataset = SpectrumDataset.create(
                    geom=on_region_geometry,
                    energy_axis_true=true_energy_axis,
                )

                spectrum_datasets_on_off = fake_onoff_from_real_observations(
                    observations,
                    empty_spectrum_dataset,
                    simulated_observation_livetime,
                    spectrum_dataset_maker,
                    background_maker,
                    sky_models,
                )

            else:
                spectrum_dataset_maker = SpectrumDatasetMaker(
                    containment_correction=containment_correction,
                    selection=["exposure", "edisp", "background"],
                )

                if gammapy_config["datasets"]["type"] == "1d":
                    sky_models[
                        target_source.name
                    ].spatial_model = PointSpatialModel.from_position(
                        on_region_geometry.region.center
                    )

                spectrum_datasets_on_off = fake_onoff_from_fake_observation(
                    observations[0],
                    empty_spectrum_dataset,
                    spectrum_dataset_maker,
                    sky_models,
                    n_off_regions,
                )

            spectrum_datasets_on_off_stacked = spectrum_datasets_on_off.stack_reduce()
            spectrum_datasets_on_off_stacked.models = sky_models
            wstat_table = get_wstat_table(spectrum_datasets_on_off_stacked)
            logging.info("Computing WSTAT on stacked spectrum ON/OFF dataset")
            print(wstat_table)

            flux_points = estimate_sed(
                target_source,
                reco_energy_axis,
                spectrum_datasets_on_off,
                **gammapy_config["estimators"]["flux_points"] or {},
            )

            plot_gammapy_sed(
                simulated_observation_livetime,
                flux_points,
                assumed_spectrum,
                plot_ts_profiles=True,
                energy_flux_unit="TeV cm-2 s-1",
                source_name=target_source.name,
                grid={"which": "both", "axis": "both", "alpha": 0.5, "ls": "dashed"},
                savefig={"format": "pdf"},
                seaborn={"context": config["seaborn_options"]["context"]},
            )

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

        if config["plotting_options"]["show"]:
            plt.show()


if __name__ == "__main__":
    main()

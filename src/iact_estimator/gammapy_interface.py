"""Interface to gammapy related to computation."""

import logging

import astropy.units as u
from astropy.coordinates import Angle
from astropy.table import Table
import numpy as np

from gammapy.data import DataStore, FixedPointingInfo, Observation, Observations
from gammapy.datasets import Datasets, SpectrumDatasetOnOff

from gammapy.maps import MapAxis, RegionGeom

from gammapy.estimators import FluxPointsEstimator
from gammapy.stats import wstat

from regions import PointSkyRegion, CircleSkyRegion

from .core import log_and_raise

logger = logging.getLogger(__name__)


def load_real_observations(data_store_path, required_irfs):
    data_store = DataStore.from_dir(data_store_path)
    observations = data_store.get_observations(required_irf=required_irfs)

    return observations


def create_simulated_observation(
    pointing_coordinates,
    observer,
    observation_livetime,
    irfs,
):
    simulated_observation = Observation.create(
        pointing=FixedPointingInfo(
            fixed_icrs=pointing_coordinates.icrs,
        ),
        livetime=observation_livetime,
        irfs=irfs,
        location=observer.location,
    )
    observations = Observations([simulated_observation])
    return observations


def has_radmax_2d(observation):
    rad_max = observation.rad_max

    return True if (rad_max is not None and rad_max.has_offset_axis) else False


def load_energy_axis_from_config(config, which):
    axis_config = config["datasets"]["geom"]["axes"][which]
    axis = MapAxis.from_energy_bounds(
        u.Quantity(axis_config["min"]),
        u.Quantity(axis_config["max"]),
        u.Quantity(axis_config["nbins"]),
        per_decade=axis_config["per_decade"],
        name=which,
    )

    return axis


def define_on_region_geometry(
    target_source_coordinates,
    observation,
    offset,
    reconstructed_energy_axis,
    on_region_radius=None,
):
    if has_radmax_2d(observation):
        on_region = PointSkyRegion(target_source_coordinates)
    else:
        if not on_region_radius:
            log_and_raise(
                logger,
                "Your IRFs do not have angular energy dependence and you did not define a fixed ON region radius.",
                exc_type=ValueError,
            )

        on_region_radius = Angle(on_region_radius)

        center = target_source_coordinates.directional_offset_by(
            position_angle=0 * u.deg, separation=Angle(offset)
        )

        on_region = CircleSkyRegion(center=center, radius=on_region_radius)

    on_region_geometry = RegionGeom.create(
        region=on_region, axes=[reconstructed_energy_axis]
    )

    return on_region_geometry


def fake_onoff_from_real_observations(
    observations,
    empty_spectrum_dataset,
    observation_livetime,
    spectrum_dataset_maker,
    background_maker,
    sky_models,
):
    """Produce ON/OFF spectrum datasets from real observations.

    This is the use case for experiments like MAGIC which provide DL3
    files with run-wise IRFs and real OFF counts in place of a background model."""

    logger.debug("Faking ON / OFF spectrum dataset from real observations.")

    fake_spectrum_datasets_on_off = Datasets()

    for observation in observations:
        spectrum_dataset = spectrum_dataset_maker.run(
            empty_spectrum_dataset.copy(name=str(observation.obs_id)), observation
        )

        spectrum_dataset.models = sky_models

        spectrum_dataset.exposure *= (
            observation_livetime / observation.observation_live_time_duration.to("h")
        )
        spectrum_dataset.fake()

        # Make a real ON/OFF dataset because we need the real OFF counts
        # since these IRFs do not have simulated background
        fake_spectrum_dataset_on_off = background_maker.run(
            spectrum_dataset, observation
        )
        fake_spectrum_dataset_on_off.models = sky_models
        # then fake it using its OFF counts
        fake_spectrum_dataset_on_off.fake(
            npred_background=fake_spectrum_dataset_on_off.npred_background()
        )
        fake_spectrum_datasets_on_off.append(fake_spectrum_dataset_on_off)

    return fake_spectrum_datasets_on_off


def fake_onoff_from_fake_observation(
    simulated_observation,
    empty_spectrum_dataset,
    spectrum_dataset_maker,
    sky_models,
    n_off_regions,
    target_source,
):
    """Produce ON/OFF spectrum datasets from a simulated observation.

    This is the use case for CTA which provide IRFs with a background model
    but no real events (for now!).
    """

    spectrum_dataset = spectrum_dataset_maker.run(
        empty_spectrum_dataset, simulated_observation
    )

    sky_models[target_source.name].spatial_model = None
    spectrum_dataset.models = sky_models
    spectrum_dataset.fake()

    spectrum_dataset_on_off = SpectrumDatasetOnOff.from_spectrum_dataset(
        dataset=spectrum_dataset, acceptance=1, acceptance_off=n_off_regions
    )
    spectrum_dataset_on_off.fake(npred_background=spectrum_dataset.npred_background())

    return Datasets(spectrum_dataset_on_off)


def estimate_sed(target_source, energy_axis, spectrum_datasets_on_off, **kwargs):
    flux_points_estimator = FluxPointsEstimator(
        source=target_source.name,
        energy_edges=energy_axis.edges,
        selection_optional="all",
        **kwargs,
    )

    flux_points = flux_points_estimator.run(spectrum_datasets_on_off)

    return flux_points


def get_wstat_table(spectrum_dataset_on_off):
    table = Table()
    table["mu_sig"] = spectrum_dataset_on_off.npred_signal().data.flatten()
    table["n_on"] = spectrum_dataset_on_off.counts.data.flatten()
    table["n_off"] = spectrum_dataset_on_off.counts_off.data.flatten()
    table["alpha"] = spectrum_dataset_on_off.alpha.data.flatten()
    table["wstat"] = wstat(
        n_on=table["n_on"],
        n_off=table["n_off"],
        alpha=table["alpha"],
        mu_sig=table["mu_sig"],
    )
    table["li_ma_sigma"] = np.sqrt(
        wstat(
            n_on=table["n_on"],
            n_off=table["n_off"],
            alpha=table["alpha"],
            mu_sig=table["mu_sig"],
        )
    )
    return table

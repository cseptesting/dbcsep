#!/usr/bin/env python3

# Copyright (C) 2024:
#   Swiss Seismological Service, ETH Zurich, Zurich, Switzerland
#   Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences, Potsdam, Germany
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero
# General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import numpy

from csep.models import EvaluationResult
from csep.core.poisson_evaluations import _poisson_likelihood_test
from csep.utils.plots import plot_likelihood_test

from dbcsep.database import SpatiaLiteCSEPDatabase


def conditional_likelihood_test(
    db: SpatiaLiteCSEPDatabase,
    plot: bool = False,
    num_simulations: int = 1000,
    seed: int = None,
    random_numbers: numpy.ndarray = None,
):
    """
    Perform Conditional Likelihood Test on a spatial database containing
     a forecast and an observed Catalogs.
    Args:
        db:
        plot:
        num_simulations:
        seed:
        random_numbers:

    Returns:

    """
    magic_array = db.spatial_magnitude_counts(1)
    observation = magic_array[:, -1]
    forecast = magic_array[:, -2]

    qs, obs_ll, simulated_ll = _poisson_likelihood_test(
        forecast,
        observation,
        num_simulations=num_simulations,
        seed=seed,
        use_observed_counts=True,
        random_numbers=random_numbers,
        verbose=True,
        normalize_likelihood=False,
    )

    result = EvaluationResult(
        test_distribution=simulated_ll,
        name="Poisson CL-Test",
        observed_statistic=obs_ll,
        quantile=qs,
        sim_name="sim",
        obs_name="cat",
        status="normal",
        min_mw=4.95,
    )

    if plot:
        plot_likelihood_test(result)

    return result


#
# def spatial_test(db: SpatiaLiteCSEPDatabase,
#                  plot: bool = False,
#                  num_simulations: int = 1000,
#                  seed: int = None,
#                  random_numbers: numpy.ndarray = None):
#     """
#     Performs the Spatial Test on the Forecast using the Observed Catalogs.
#
#     Note: The forecast and the observations should be scaled to the same time period before calling this function. This increases
#     transparency as no assumptions are being made about the length of the forecasts. This is particularly important for
#     gridded forecasts that supply their forecasts as rates.
#
#     Args:
#         gridded_forecast: csep.core.forecasts.GriddedForecast
#         observed_catalog: csep.core.catalogs.Catalog
#         num_simulations (int): number of simulations used to compute the quantile score
#         seed (int): used fore reproducibility, and testing
#         random_numbers (numpy.ndarray): random numbers used to override the random number generation. injection point for testing.
#
#     Returns:
#         evaluation_result: csep.core.evaluations.EvaluationResult
#     """
#
#     gridded_catalog_data = observed_catalog.spatial_counts()
#
#     # simply call likelihood test on catalog data and forecast
#     qs, obs_ll, simulated_ll = _poisson_likelihood_test(
#         gridded_forecast.spatial_counts(), gridded_catalog_data,
#         num_simulations=num_simulations,
#         seed=seed,
#         random_numbers=random_numbers,
#         use_observed_counts=True,
#         verbose=verbose,
#         normalize_likelihood=True)

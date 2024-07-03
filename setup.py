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

from setuptools import setup, find_packages

tests_require = ["pytest"]
linters_require = ["black>=20.8b1", "pylint", "flake8"]

setup(
    name="databaselib",
    description="Database library for handling Spatialite and PostGIS databases",
    license="AGPLv3+",
    install_requires=["packaging", "obspy", "mercantile", "shapely", "numpy", "pycsep"],
    extras_require={
        "tests": tests_require,
        "linters": linters_require,
    },
    packages=find_packages(),
    python_requires=">=3.11",
)

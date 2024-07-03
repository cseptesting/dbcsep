#!/usr/bin/env python3

# Copyright (C) 2024:
#   Swiss Seismological Service, Zurich, ETH Zurich.
#   Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences, Potsdam, Germany.
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

from abc import abstractmethod, ABC

from database import AbstractDatabase


class AbstractManager(AbstractDatabase, ABC):
    """
    Todo: Fix docstring
    The `AbstractExposure` class contains all functions properties that are common to
    the derived exposure classes `SpatiaLiteExposure` and `PostGISExposure`.
    """

    def __init__(self):
        super().__init__()


class SpatiaLiteManager(SpatiaLiteDatabase, AbstractManager):
    """
    Todo: Fix docstring
    The `SpatiaLiteExposure` class represents a database for exposure models.
    It is derived from the generic `SpatiaLiteDatabase` class for providing functionality
    related to SpatiaLite databases and the `AbstractExposure` class for database-type
    independent functionality regarding handling of exposure data.
    """

    def __init__(self, database_filepath, spatialite_filepath="mod_spatialite"):
        """
        Initializes the SpatiaLite database.

        Args:
            database_filepath(str):
                Filepath of the exposure database file.
            spatialite_filepath(str):
                Filepath of the SpatiaLite extension.
        """

        super().__init__(database_filepath, spatialite_filepath=spatialite_filepath)

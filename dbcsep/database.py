#!/usr/bin/env python3

# Copyright (C) 2024:
#   Swiss Seismological Service, ETH Zurich, Zurich, Switzerland.
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

import abc
import logging
import sqlite3

logger = logging.getLogger()


class AbstractDatabase(abc.ABC):
    """
    The AbstractDatabase class represents any database. It manages the database connection and
    cursor.

    Attributes:
        connection:
            Database connection.
        cursor:
            Database cursor.
    """

    def __init__(self):
        self.connection = None
        self.cursor = None

    @abc.abstractmethod
    def connect(self):
        """
        Abstract connection to the database.
        """

        pass

    def commit_and_close(self):
        """
        Commits the executed statements and closes the database connection.
        """

        self.connection.commit()
        self.close()

    def close(self):
        """
        Closes the database connection.
        """

        self.connection.close()


class SpatiaLiteDatabase(AbstractDatabase):
    """
    The `SpatiaLiteDatabase` class represents a SpatiaLite database. It manages the database
    connection and cursor. Upon connection, the SpatiaLite extension is loaded.

    Args:
        database_filepath (str):
            Filepath for the SpatiaLite database file.
        spatialite_filepath (str, default: 'mod_spatialite'):
            Filepath of the SpatiaLite extension.

    Attributes:
        database_filepath (str):
            SpatiaLite database filepath.
        spatialite_filepath (str):
            Filepath to the SpatiaLite extension.
        connection (sqlite3.Connection):
            Database connection.
        cursor (sqlite3.Cursor):
            Database cursor.
    """

    def __init__(self, database_filepath, spatialite_filepath="mod_spatialite"):
        super().__init__()
        self.database_filepath = database_filepath
        self.spatialite_filepath = spatialite_filepath

    def connect(self, init_spatial_metadata=True, **kwargs):
        """
        Connects to a database, loads the SpatiaLite extension and initializes it if requested.
        If the database file does not exist, it will be automatically created.

        Args:
            init_spatial_metadata (Bool):
                Defines if the spatial metadata should be initialized. This is not necessary if
                an existing SpatiaLite database is opened.
        """

        # Connect (if exists) or create SQLite database.
        logger.debug(f"Connecting to database at {self.database_filepath}.")
        try:
            self.connection = sqlite3.connect(self.database_filepath, **kwargs)
        except Exception:
            logger.exception("Error connecting to the database file.")
            raise
        logger.debug("Connection to database established.")

        # Load and initialize the SpatiaLite extension.
        self.connection.enable_load_extension(True)
        sql_statement = f"SELECT load_extension('{self.spatialite_filepath}');"
        try:
            self.connection.execute(sql_statement)
        except sqlite3.OperationalError:
            logger.exception("Error loading spatial extension.")
            raise
        logger.debug("SpatiaLite extension loaded.")

        if init_spatial_metadata:
            try:
                self.connection.execute("SELECT InitSpatialMetaData(1);")
            except Exception:
                logger.exception("Error initializing spatial metadata.")
                raise
            logger.debug("Spatial metadata initialized.")

        # Debug output.
        try:
            versions = self.connection.execute("SELECT sqlite_version(), spatialite_version()")
        except Exception:
            logger.exception("Error querying database versions.")
            raise
        for row in versions:
            logger.debug(f"SQLite version: {row[0]}.")
            logger.debug(f"SpatiaLite version: {row[1]}.")

        try:
            self.connection.execute("PRAGMA foreign_keys = ON")
        except Exception:
            logger.exception("Error with enabling foreign keys.")
            raise

        self.cursor = self.connection.cursor()

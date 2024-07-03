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
# /home/danijel/dev/pmc-lib
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import abc
import logging
import sqlite3
import csv
import numpy
import mercantile

import dbcsep.utils as utils

import time
from obspy import UTCDateTime, Catalog
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException

logger = logging.getLogger()


class AbstractDatabase(abc.ABC):
    """
    The AbstractDatabase class represents any database. It manages the
    database connection and cursor.

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
    The `SpatiaLiteDatabase` class represents a SpatiaLite database. It manages the
    database connection and cursor. Upon connection, the SpatiaLite extension is loaded.

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
            Database connection
        cursor (sqlite3.Cursor):
            Database cursor
    """

    def __init__(self, database_filepath, spatialite_filepath="mod_spatialite"):
        super().__init__()
        self.database_filepath = database_filepath
        self.spatialite_filepath = spatialite_filepath

    def connect(self, init_spatial_metadata=True, **kwargs):
        """
        Connects to a database, loads the SpatiaLite extension and
        initializes it if requested. If the database file does not exist,
        it will be automatically created.

        Args:
            init_spatial_metadata (Bool):
                Defines if the spatial metadata should be initialized. This is not
                necessary if an existing SpatiaLite database is opened.
        """

        # Connect (if exists) or create SQLite database
        logger.debug("Connecting to database at %s" % self.database_filepath)
        try:
            self.connection = sqlite3.connect(self.database_filepath, **kwargs)
        except Exception:
            logger.exception("Error connecting to the database file")
            raise
        logger.debug("Connection to database established")

        # Load and initialize the SpatiaLite extension
        self.connection.enable_load_extension(True)
        sql_statement = "SELECT load_extension('%s');" % self.spatialite_filepath
        try:
            self.connection.execute(sql_statement)
        except sqlite3.OperationalError:
            logger.exception("Error loading spatial extension")
            raise
        logger.debug("SpatiaLite extension loaded")

        if init_spatial_metadata:
            try:
                self.connection.execute("SELECT InitSpatialMetaData(1);")
            except Exception:
                logger.exception("Error initializing spatial metadata")
                raise
            logger.debug("Spatial metadata initialized")

        # Debug output
        try:
            versions = self.connection.execute("SELECT sqlite_version(), spatialite_version()")
        except Exception:
            logger.exception("Error querying database versions")
            raise
        for row in versions:
            logger.debug("SQLite version: %s" % row[0])
            logger.debug("SpatiaLite version: %s" % row[1])

        try:
            self.connection.execute("PRAGMA foreign_keys = ON")
        except Exception:
            logger.exception("Error with enabling foreign keys")
            raise

        self.cursor = self.connection.cursor()


class SpatiaLiteCSEPDatabase(SpatiaLiteDatabase):
    """
    The `SpatiaLiteCatalogDatabase` class represents a SpatiaLite database to store and provide
    catalog data for the FTLS.

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
        super().__init__(database_filepath, spatialite_filepath)

    def create_tables(self):
        """
        Creates all necessary tables in the database. These are:
            Events  : Stores event parameters.
            Metadata: Stores metadata tags as key-value pairs.
        """

        # Create table `Catalogs`.
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Catalogs (
                id               INTEGER PRIMARY KEY AUTOINCREMENT,
                start_timestamp  REAL,
                end_timestamp    REAL,
                description      TEXT
            );
            """
        self.connection.executescript(sql_statement)
        logger.debug("Table `Catalogs` created.")

        # Create table `Events`.
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Events (
                id           INTEGER PRIMARY KEY AUTOINCREMENT,
                catalog_id   INTEGER,
                resource_id  TEXT,
                magnitude    REAL,
                origin_time  INTEGER,
                UNIQUE(resource_id)
            );
            SELECT AddGeometryColumn('Events', 'geom', 4979, 'POINT', 'XYZ');
            SELECT AddGeometryColumn('Events', 'centroid', 4326, 'POINT', 'XY');
            SELECT CreateSpatialIndex('Events', 'geom');
            SELECT CreateSpatialIndex('Events', 'centroid');
            CREATE INDEX IF NOT EXISTS idx_events_catalogid ON Events(catalog_id);
            CREATE INDEX IF NOT EXISTS idx_events_origintime ON Events(origin_time);
            """
        self.connection.executescript(sql_statement)
        logger.debug("Table `Events` created.")

        # Create table `Metadata`.
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Metadata (
                tag_key      TEXT,
                tag_value    TEXT,
                tag_comment  TEXT,
                UNIQUE(tag_key)
            )
            """
        self.connection.executescript(sql_statement)
        logger.debug("Table `Metadata` created.")

        # Create table `Models`
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Models (
                id           INTEGER PRIMARY KEY AUTOINCREMENT,
                grid_id      INTEGER,
                description  TEXT
            );"""
        self.connection.executescript(sql_statement)
        logger.debug("Table `Models` created.")

        # Create table `Forecasts`
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Forecasts (
                id           INTEGER PRIMARY KEY AUTOINCREMENT,
                model_id  INTEGER,
                start_timestamp  REAL,
                end_timestamp REAL
            );"""
        self.connection.executescript(sql_statement)
        logger.debug("Table `Forecasts` created.")

        # Create table `Grids`
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Grids (
                id           INTEGER PRIMARY KEY AUTOINCREMENT,
                description  TEXT
            );"""
        self.connection.executescript(sql_statement)
        logger.debug("Table `Grids` created.")

        # Create table `Cells`
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Cells (
                id           INTEGER PRIMARY KEY AUTOINCREMENT,
                grid_id      INTEGER,
                quadkey      TEXT
            );
            SELECT AddGeometryColumn('Cells', 'geom', 4326, 'POLYGON', 'XY');
            CREATE INDEX IF NOT EXISTS idx_cells_quadkey ON Cells(quadkey COLLATE NOCASE);
            SELECT CreateSpatialIndex('Cells', 'geom');
            """
        self.connection.executescript(sql_statement)
        logger.debug("Table `Cells` created")

        # Create table `Bins`
        sql_statement = """
            CREATE TABLE IF NOT EXISTS Bins (
                id             INTEGER PRIMARY KEY AUTOINCREMENT,
                grid_id        INTEGER,
                min_magnitude  REAL,
                max_magnitude  REAL
                );
            """
        self.connection.executescript(sql_statement)
        logger.debug("Table `Bins` created.")

        # Create table `Forecasts`
        sql_statement = """
            CREATE TABLE IF NOT EXISTS ForecastRates (
                id              INTEGER PRIMARY KEY AUTOINCREMENT,
                forecast_id  INTEGER,
                cell_id         INTEGER,
                bin_id          INTEGER,
                rate     REAL
            );
            CREATE INDEX IF NOT EXISTS idx_forecasts_forecastid ON ForecastRates(forecast_id);
            CREATE INDEX IF NOT EXISTS idx_forecasts_cellid ON ForecastRates(cell_id);
            CREATE INDEX IF NOT EXISTS idx_forecasts_binid ON ForecastRates(bin_id);
            """
        self.connection.executescript(sql_statement)
        logger.debug("Table `ForecastRates` created.")

    def insert_catalog(self, dd):
        """"""

        return True

    def insert_event(
        self,
        resource_id: str,
        magnitude: float,
        origin_time: float,
        longitude: float,
        latitude: float,
        elevation: float,
    ) -> int:
        """
        Inserts an event dataset into the `Events` table and returns its ID. The dataset
        consists of the resource ID, magnitude, timestamp of the origin time, and the longitude,
        latitude and elevation of the hypocenter. The vertical coordinate of events is usually
        given as depth and need to be given as negative values (for elevation). In case the
        UNIQUE constraint on the field `resource_id` is violated, the function raises a
        `SystemExit` exception and reports the resource ID causing the exception.

        Args:
            resource_id (str):
                Resource ID of the event.
            magnitude (float):
                Magnitude of the event.
            origin_time (float):
                Timestamp of the origin time of the event.
            longitude (float):
                Longitude of the event hypocenter in degrees (WGS 84).
            latitude (float):
                Latitude of the event hypocenter in degrees (WGS 84).
            elevation (float):
                Elevation of the event hypocenter in meters. Depths are given with negative
                values.

        Returns:
            ID of the event.
        """

        utils.check_is_number_between(magnitude, -10, 10)
        utils.check_is_number_between(longitude, -180, 180)
        utils.check_is_number_between(latitude, -90, 90)
        # Elevation/depth is between center of the Earth and the highest peak.
        utils.check_is_number_between(elevation, -6738000, 8848)
        sql_statement = f"""
            INSERT INTO Events (resource_id, magnitude, origin_time, geom, centroid)
            VALUES (
                '{resource_id}',
                {magnitude},
                {origin_time},
                GeomFromText('POINTZ({longitude} {latitude} {elevation})', 4979),
                GeomFromText('POINT({longitude} {latitude})', 4326)
            )
            RETURNING id
            """
        try:
            self.cursor.execute(sql_statement)
        except sqlite3.IntegrityError:
            raise SystemExit(
                f"Resource ID '{resource_id}' is already present in the `Events` table. "
                f"Adding this event is prevented by the UNIQUE constraint on `resource_id`. "
                f"Exiting."
            )
        return self.cursor.fetchone()[0]

    def get_event(self, event_id: int) -> tuple:
        """
        Retrieves an event dataset from the `Events` table given the ID of the event.

        Args:
            event_id (int):
                ID of the event.

        Returns:
            A tuple containing the ID, network, station, and location code of the station, the
            installation and de-installation timestamp of the station, the longitude, latitude,
            elevation and geometry in WKT of the station.
        """

        sql_statement = f"""
            SELECT id, resource_id, magnitude, origin_time, ST_X(geom), ST_Y(geom), ST_Z(geom),
                ST_AsText(geom)
            FROM Events
            WHERE id = {event_id}
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()

    def get_event_by_resource_id(self, resource_id: str) -> tuple:
        """
        Retrieves an event dataset from the `Events` table given the resource ID of the event.
        Resource IDs are part of the QuakeMl definition and are defined as globally unique for
        each event.

        Args:
            resource_id (str):
                Globally unique resource ID of the event, see the QuakeML documentation for
                further information.

        Returns:
            A tuple containing the ID, network, station, and location code of the station, the
            installation and de-installation timestamp of the station, the longitude, latitude,
            elevation and geometry in WKT of the station.
        """

        sql_statement = f"""
            SELECT id, resource_id, magnitude, origin_time, ST_X(geom), ST_Y(geom), ST_Z(geom),
                ST_AsText(geom)
            FROM Events
            WHERE resource_id = '{resource_id}'
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()

    def get_number_of_events(self) -> int:
        """"""

        sql_statement = "SELECT COUNT(id) FROM Events"
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()[0]

    def upsert_tag(self, key: str, value: str, comment: str = "", overwrite: bool = False):
        """
        Inserts or updates a tag as key-value pair into the `Metadata` table. If the key does
        not yet exist, the key-value pair is inserted. Otherwise, the value of the existing key
        will be updated if `overwrite` is set to True. If `overwrite` is set to False, the
        function will keep the existing key-value pair.

        Args:
            key (str):
                The key describing the metadata topic.
            value (str):
                The metadata value of the given topic.
            comment (str, optional, default: Empty string):
                Optional comment to be added to the key-value pair.
            overwrite (bool, optional, default: True):
                If set, an existing metadata value for the given key will be overwritten. If not
                set, the function will keep the existing key-value pair. In case the key does
                not yet exist in the table, this option has no effect.
        """

        if not overwrite:
            sql_statement = f"SELECT * FROM Metadata WHERE tag_key = '{key}' LIMIT 1"
            self.cursor.execute(sql_statement)
            if self.cursor.fetchone() is not None:
                logger.warning(f"Metadata key {key} already exists.")
                return

        sql_statement = f"""
            INSERT INTO Metadata (tag_key, tag_value, tag_comment)
            VALUES ('{key}', '{value}', '{comment}')
            ON CONFLICT (tag_key) DO UPDATE SET
                tag_value = EXCLUDED.tag_value,
                tag_comment = EXCLUDED.tag_comment
            """
        self.cursor.execute(sql_statement)

    def get_tag(self, key: str) -> tuple | None:
        """
        Retrieves the value of a tag (key-value pair) from the `Metadata` table given its key.

        Args:
            key (str):
                The key describing the metadata topic for which the value should be retrieved.

        Returns
            A tuple containing the metadata value and comment of the given key.
        """

        sql_statement = f"""
            SELECT tag_value, tag_comment FROM Metadata
            WHERE tag_key = '{key}'
            """
        self.cursor.execute(sql_statement)
        result = self.cursor.fetchone()
        if result is None:
            logger.debug(f"No tag with key {key} exists.")
            return None
        else:
            return result

    def download_catalog_from_fdsn_service(
        self,
        base_url: str,
        start_time: int | float | str,
        end_time: int | float | str,
        min_magnitude: float,
        min_longitude: float,
        max_longitude: float,
        min_latitude: float,
        max_latitude: float,
        max_attempts: int = 3,
    ) -> int | None:
        """
        Downloads an earthquake catalog from a service of the International Federation of
        Digital Seismograph Networks (FDSN) and stores the event data in the `Events` table. The
        dataset to be downloaded is specified by a start and end time, the minimum magnitude and
        a bounding box given minimum and maximum longitudes and latitudes. The service to
        download from is specified in the `base_url` parameter. After a successful download, the
        download parameters are stored in as metadata tags in the `Metadata` table to be used
        for the updating of the catalog at later times. In case no catalog data exists for the
        given download parameters or the maximum number of download attempts has been reached,
        the function returns `False`, otherwise `True`.

        Args:
            base_url (str):
                URL of a FDSN web service or a shortcut name (e.g. 'USGS'). A full list of
                available providers can be found on https://docs.obspy.org.
            start_time (int|float|str):
                Start time of the catalog to be downloaded. The time can be given as a UNIX
                timestamp (as `int` or `float`) or a ISO8601:2004 string (e.g.
                '2010-01-01T12:00:00+01:00').
            end_time (int|float|str):
                End time of the catalog to be downloaded. The time can be given as a UNIX
                timestamp (as `int` or `float`) or a ISO8601:2004 string (e.g.
                '2010-01-01T12:00:00+01:00').
            min_magnitude (float):
                Minimum magnitude of the catalog to be downloaded.
            min_longitude (float):
                Minimum longitude of the catalog to be downloaded.
            max_longitude (float):
                Maximum longitude of the catalog to be downloaded.
            min_latitude (float):
                Minimum latitude of the catalog to be downloaded.
            max_latitude (float):
                Maximum latitude of the catalog to be downloaded.
            max_attempts (int, optional, default: 3):
                Maximum number of attempts to download the catalog.

        Returns:
            Todo: Fix
            `True` if the catalog download was successful, `False` if no data is matching the
            download parameters or the maximum number of unsuccessful attempts is reached.
        """


        if self.get_number_of_events() > 0:
            raise Exception("Catalog is not empty. Cannot use the download function.")

        start_time = UTCDateTime(start_time)
        end_time = UTCDateTime(end_time)
        utils.check_is_number_between(min_magnitude, -10, 10)
        utils.check_is_number_between(min_longitude, -180, 180, "left")
        utils.check_is_number_between(max_longitude, -180, 180, "right")
        utils.check_is_number_between(min_latitude, -90, 90, "left")
        utils.check_is_number_between(max_latitude, -90, 90, "right")

        client = Client(base_url)
        client.timeout = 180

        # Repeat until the request has been successfully completed or the maximum number of
        # attempts is reached.
        for attempts in range(0, max_attempts):
            try:
                logger.debug(f"Downloading from FDSN service (attempt no. {attempts})...")
                catalog: Catalog = client.get_events(
                    starttime=start_time,
                    endtime=end_time,
                    minmagnitude=min_magnitude,
                    minlongitude=min_longitude,
                    maxlongitude=max_longitude,
                    minlatitude=min_latitude,
                    maxlatitude=max_latitude,
                    eventtype="earthquake",
                )
                last_event_timestamp = start_time.timestamp
                logger.debug("Importing data from FDSN service to the database...")
                for event in catalog:
                    origin = event.preferred_origin()
                    # Insert the event into the `Events` table.
                    self.insert_event(
                        event.resource_id,
                        event.preferred_magnitude().mag,
                        origin.time.timestamp,
                        origin.longitude,
                        origin.latitude,
                        -1 * origin.depth,  # Depths are negative in EPSG:4979.
                    )
                    if origin.time.timestamp > last_event_timestamp:
                        last_event_timestamp = origin.time.timestamp
                self.upsert_tag("Catalog.BaseURL", base_url)
                self.upsert_tag("Catalog.MinMagnitude", str(min_magnitude))
                self.upsert_tag("Catalog.MinLongitude", str(min_longitude))
                self.upsert_tag("Catalog.MaxLongitude", str(max_longitude))
                self.upsert_tag("Catalog.MinLatitude", str(min_latitude))
                self.upsert_tag("Catalog.MaxLatitude", str(max_latitude))
                self.upsert_tag("Catalog.LastEventTimestamp", str(last_event_timestamp))
                return True
            except FDSNNoDataException:
                logger.error("No data available for the parameters provided.")
                return False
            except Exception as e:
                logger.error(e)
                time.sleep(3)
        return False

    def update_catalog_from_fdsn_service(self, max_attempts: int = 3) -> bool:
        """
        Downloads new data to the already existing catalog of events in the `Events` table from
        the respective FDSN web service. All relevant download parameters are taken from the
        `Metadata` table where the `download_catalog_from_fdsn_service()` function has stored
        them after a successful initial download of the catalog. One of the parameters stored is
        the timestamp of the last event in the catalog. This timestamp is used to define the
        time from which the catalog should be updated. In case no catalog data exists from this
        time on or the maximum number of download attempts has been reached, the function
        returns `False`, otherwise `True`.

        Args:
            max_attempts (int, optional, default: 3):
                Maximum number of attempts to download the catalog.

        Returns:
            `True` if the catalog download was successful, `False` if no data is matching the
            download parameters or the maximum number of unsuccessful attempts is reached.

        """

        base_url = self.get_tag("Catalog.BaseURL")[0]
        update_from_timestamp = float(self.get_tag("Catalog.LastEventTimestamp")[0])
        min_magnitude = float(self.get_tag("Catalog.MinMagnitude")[0])
        min_longitude = float(self.get_tag("Catalog.MinLongitude")[0])
        max_longitude = float(self.get_tag("Catalog.MaxLongitude")[0])
        min_latitude = float(self.get_tag("Catalog.MinLatitude")[0])
        max_latitude = float(self.get_tag("Catalog.MaxLatitude")[0])

        client = Client(base_url)
        client.timeout = 180

        # Repeat until the request has been successfully completed or the maximum number of
        # attempts is reached.
        for attempts in range(0, max_attempts):
            try:
                logger.debug(f"Downloading from FDSN service (attempt no. {attempts})...")
                catalog = client.get_events(
                    updatedafter=update_from_timestamp,
                    minmagnitude=min_magnitude,
                    minlongitude=min_longitude,
                    maxlongitude=max_longitude,
                    minlatitude=min_latitude,
                    maxlatitude=max_latitude,
                    eventtype="earthquake",
                )
                last_event_timestamp = update_from_timestamp
                logger.debug("Importing data from FDSN service to the database...")
                for event in catalog:
                    origin = event.preferred_origin()
                    # Insert the event into the `Events` table.
                    self.insert_event(
                        event.resource_id,
                        event.preferred_magnitude().mag,
                        origin.time.timestamp,
                        origin.longitude,
                        origin.latitude,
                        -1 * origin.depth,  # Depths are negative in EPSG:4979.
                    )
                    if origin.time.timestamp > last_event_timestamp:
                        last_event_timestamp = origin.time.timestamp
                self.upsert_tag(
                    "Catalog.LastEventTimestamp", str(last_event_timestamp), overwrite=True
                )
                return True
            except FDSNNoDataException:
                logger.error("No data available for the parameters provided.")
                return False
            except Exception as e:
                logger.error(e)
                time.sleep(3)
        return False

    def insert_grid(self, description):
        """ """

        sql_statement = f"""
            INSERT INTO Grids (description)
            VALUES ('{description}')
            RETURNING id
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()[0]

    def insert_cell(
        self, grid_id: int, min_lon: float, min_lat: float, max_lon: float, max_lat: float
    ) -> int:
        """"""

        cell_wkt = f"""
            POLYGON((
                {min_lon} {min_lat},
                {min_lon} {max_lat},
                {max_lon} {max_lat},
                {max_lon} {min_lat},
                {min_lon} {min_lat}
            ))
            """
        sql_statement = f"""
            INSERT INTO Cells (grid_id, quadkey, geom)
            VALUES (
                {grid_id},
                NULL,
                GeomFromText('{cell_wkt}', 4326)
            )
            RETURNING id
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()[0]

    def get_cell_id(self, grid_id: int, min_lon: float, min_lat: float) -> int | None:
        """
        Todo: Fix docstring, Make the ST_MinX, etc. functions compatible with PostGIS

        Args:
            grid_id:
            min_lon:
            min_lat:

        Returns:
            ID of the grid cell.
        """

        sql_statement = f"""
            SELECT id
            FROM Cells
            WHERE grid_id = {grid_id}
                AND ST_MinX(geom) = {min_lon}
                AND ST_MinY(geom) = {min_lat}
            """
        self.cursor.execute(sql_statement)
        result = self.cursor.fetchone()
        if result is None:
            return None
        else:
            return result[0]

    def insert_bin(self, grid_id: int, min_magnitude: float, max_magnitude: float) -> int:
        """ """

        sql_statement = f"""
            INSERT INTO Bins (grid_id, min_magnitude, max_magnitude)
            VALUES ({grid_id}, {min_magnitude}, {max_magnitude})
            """
        self.cursor.execute(sql_statement)

    def get_bin_id(self, grid_id: int, min_magnitude: float) -> int | None:
        """"""

        sql_statement = f"""
            SELECT id FROM Bins
            WHERE grid_id = {grid_id} AND min_magnitude = {min_magnitude} 
            """
        self.cursor.execute(sql_statement)
        result = self.cursor.fetchone()
        if result is None:
            return None
        else:
            return result[0]

    def insert_computation(self, grid_id: int, description: str = "") -> int:
        """"""

        sql_statement = f"""
            INSERT INTO Computations (grid_id, description)
            VALUES ({grid_id}, '{description}')
            RETURNING id
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()[0]

    def insert_model(self, description: str, grid_id: int) -> int:
        """"""

        sql_statement = f"""
            INSERT INTO Models (description, grid_id)
            VALUES ('{description}', {grid_id})
            RETURNING id
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()[0]

    def insert_forecast(self, model_id, start_time, end_time):
        """"""

        start_time = UTCDateTime(start_time)
        end_time = UTCDateTime(end_time)
        sql_statement = f"""
            INSERT INTO Forecasts (model_id, start_timestamp, end_timestamp)
            VALUES ({model_id}, {start_time.timestamp}, {end_time.timestamp})
            RETURNING id
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchone()[0]

    def get_grid_id(self, computation_id: int) -> int | None:
        """"""

        sql_statement = f"""
            SELECT grid_id FROM Computations WHERE id = {computation_id}
            """
        self.cursor.execute(sql_statement)
        result = self.cursor.fetchone()
        if result is None:
            return None
        else:
            return result[0]

    def load_csep_grid(self, grid_filepath: str, description: str = "") -> int:
        """
        Todo: Fix docstring

        Use a standard CSEP forecast to get the cells for the grid.
        Args:
            grid_filepath:
            description:

        Returns:
            ID of the grid.

        """

        grid_id = self.insert_grid(description)
        logger.info("Loading CSEP grid...")
        with open(grid_filepath) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter="\t")
            min_lon = max_lon = min_lat = max_lat = 0
            for row in csv_reader:
                # Test if the current line is describing the same cell as the previous line.
                if not (
                    (min_lon == float(row[0]))
                    and (max_lon == float(row[1]))
                    and (min_lat == float(row[2]))
                    and (max_lat == float(row[3]))
                ):
                    min_lon = float(row[0])
                    max_lon = float(row[1])
                    min_lat = float(row[2])
                    max_lat = float(row[3])
                    self.insert_cell(grid_id, min_lon, min_lat, max_lon, max_lat)
                min_mag = float(row[6])
                max_mag = float(row[7])
                if self.get_bin_id(grid_id, min_mag) is None:
                    self.insert_bin(grid_id, min_mag, max_mag)
            self.connection.commit()
        logger.info("CSEP grid loaded.")
        return grid_id

    def create_grid(self, min_lon, max_lon, min_lat, max_lat, zoom_level, description=""):
        """
        depth in meter
        """

        grid_id = self.insert_grid(description)
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, zoom_level)
        for tile in tiles:
            bb = mercantile.bounds(tile)
            tile_wkt = f"""
                POLYGON(({bb.west} {bb.south}, {bb.west} {bb.north}, {bb.east} {bb.north},
                    {bb.east} {bb.south}, {bb.west} {bb.south}))
                """
            sql_statement = f"""
                INSERT INTO Cells (grid_id, quadkey, geom)
                VALUES (
                    {grid_id},
                    '{mercantile.quadkey(tile)}',
                    GeomFromText('{tile_wkt}', 4326)
                )
                """
            self.cursor.execute(sql_statement)
        self.connection.commit()
        return grid_id

    def get_nodes(self, grid_id):
        """ """

        sql_statement = f"""
            SELECT quadkey, ST_X(geom), ST_Y(geom), ST_Z(geom) FROM Node
            WHERE grid_id = {grid_id}
            """
        self.cursor.execute(sql_statement)
        return self.cursor.fetchall()

    def insert_forecast_rate(self, forecast_id, cell_id, bin_id, rate):
        """ """

        sql_statement = f"""
            INSERT INTO ForecastRates (forecast_id, cell_id, bin_id, rate)
            VALUES ({forecast_id}, {cell_id}, {bin_id}, {rate})
            """
        self.cursor.execute(sql_statement)

    def load_csep_forecast(self, forecast_id: int, grid_id: int, forecast_filepath: str):
        """ """
        if grid_id is None:
            raise Exception("No grid given.")

        with (open(forecast_filepath) as csv_file):
            csv_reader = csv.reader(csv_file, delimiter="\t")
            min_lon = min_lat = 0
            for row in csv_reader:
                if not ((min_lon == float(row[0])) and (min_lat == float(row[2]))):
                    min_lon = float(row[0])
                    min_lat = float(row[2])
                    cell_id = self.get_cell_id(grid_id, min_lon, min_lat)
                min_mag = float(row[6])
                bin_id = self.get_bin_id(grid_id, min_mag)
                rate = float(row[8])
                self.insert_forecast_rate(forecast_id, cell_id, bin_id, rate)
            self.connection.commit()

    def spatial_magnitude_counts(self, forecast_id):

        sql_statement = f"""
            SELECT
                B.id,
                rate AS forecast,
                COUNT(Events.id) AS observation
            FROM
            (
                SELECT 
                    ForecastRates.id,
                    Cells.geom,
                    Bins.min_magnitude,
                    Bins.max_magnitude,
                    ForecastRates.rate
                FROM ForecastRates
                INNER JOIN Cells ON ForecastRates.cell_id = Cells.id
                INNER JOIN Bins ON ForecastRates.bin_id = Bins.id
                WHERE forecast_id = {forecast_id}
            ) AS B
            LEFT OUTER JOIN Events 
                ON (ST_Contains(B.geom, Events.centroid) 
                    AND Events.magnitude >= B.min_magnitude
                    AND Events.magnitude < B.max_magnitude)
                    AND ST_X(Events.centroid) < ST_MaxX(B.geom)
                    AND ST_Y(Events.centroid) < ST_MaxY(B.geom) 
            GROUP BY B.id
            """

        # sql_statement = """
        #     SELECT
        #         Bins.id,
        #         binid,
        #         min_magnitude,
        #         probability AS forecast,
        #         COUNT(Events.id) AS observation
        #     FROM
        #     (
        #             SELECT
        #             Node.id,
        #             NodeProbability.id as binid,
        #             Node.geom, NodeProbability.min_magnitude,
        #             NodeProbability.max_magnitude, probability
        #         FROM Node
        #         INNER JOIN NodeProbability ON Node.id = NodeProbability.tile_id
        #     ) AS Bins
        #     LEFT OUTER JOIN Events
        #         ON (ST_Contains(Bins.geom, Events.centroid)
        #             AND Events.magnitude >= Bins.min_magnitude
        #             AND Events.magnitude < Bins.max_magnitude)
        #             AND ST_X(Events.centroid) < ST_MaxX(Bins.geom)
        #             AND ST_Y(Events.centroid) < ST_MaxY(Bins.geom)
        #     GROUP BY binid
        #     """
        self.cursor.execute(sql_statement)
        return numpy.array(self.cursor.fetchall())

    def spatial_counts(self):

        sql_statement = """
            SELECT MainQuery.id, SUM(rate) AS forecast, SUM(events_per_bin) AS observation
            FROM
            (            
                SELECT B.id, rate, COUNT(Events.id) AS events_per_bin
                FROM
                (
                    SELECT 
                        ForecastRates.id,
                        Cells.geom,
                        Bins.min_magnitude,
                        Bins.max_magnitude,
                        ForecastRates.rate
                    FROM ForecastRates
                    INNER JOIN Cells ON ForecastRates.cell_id = Cells.id
                    INNER JOIN Bins ON ForecastRates.bin_id = Bins.id
                    WHERE forecast_id = 1
                ) AS B
                LEFT OUTER JOIN Events 
                    ON (ST_Contains(B.geom, Events.centroid) 
                        AND Events.magnitude >= B.min_magnitude
                        AND Events.magnitude < B.max_magnitude)
                        AND ST_X(Events.centroid) < ST_MaxX(B.geom)
                        AND ST_Y(Events.centroid) < ST_MaxY(B.geom) 
                GROUP BY B.id
            ) AS MainQuery
            GROUP BY MainQuery.id
            """
        self.cursor.execute(sql_statement)
        return numpy.array(self.cursor.fetchall())

    def magnitude_counts(self):

        sql_statement = """
            SELECT min_magnitude, SUM(probability) AS forecast, SUM(eventsperbin) AS observation FROM
            (            
                SELECT B.id, rate, COUNT(Events.id) AS events_per_bin
                FROM
                (
                    SELECT 
                        ForecastRates.id,
                        Cells.geom,
                        Bins.min_magnitude,
                        Bins.max_magnitude,
                        ForecastRates.rate
                    FROM ForecastRates
                    INNER JOIN Cells ON ForecastRates.cell_id = Cells.id
                    INNER JOIN Bins ON ForecastRates.bin_id = Bins.id
                    WHERE forecast_id = 1
                ) AS B
                LEFT OUTER JOIN Events 
                    ON (ST_Contains(B.geom, Events.centroid) 
                        AND Events.magnitude >= B.min_magnitude
                        AND Events.magnitude < B.max_magnitude)
                        AND ST_X(Events.centroid) < ST_MaxX(B.geom)
                        AND ST_Y(Events.centroid) < ST_MaxY(B.geom) 
                GROUP BY B.id
) AS MainQuery
GROUP BY MainQuery.min_magnitude    """
        self.cursor.execute(sql_statement)
        return numpy.array(self.cursor.fetchall())

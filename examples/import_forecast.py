#!/usr/bin/env python3

# Copyright (C) 2024:
#   Swiss Seismological Service, ETH Zurich, Zurich, Switzerland
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

import os
import sys
import logging

from dbcsep import SpatiaLiteCSEPDatabase
from dbcsep import conditional_likelihood_test

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler(sys.stdout))


def main():
    """ """

    os.remove("test.db")
    db = SpatiaLiteCSEPDatabase("test.db")
    db.connect()
    db.create_tables()
    grid_id = db.load_csep_grid("helmstetter_et_al.hkj-fromXML.dat")
    model_id = db.insert_model("Helmstetter et al.", grid_id)
    forecast_id = db.insert_forecast(model_id, "2006-01-01", "2011-01-01")
    db.load_csep_forecast(forecast_id, grid_id, "helmstetter_et_al.hkj-fromXML.dat")
    db.download_catalog_from_fdsn_service("USGS", "2023-01-01", "2024-01-01", 4.95, -126, -111, 30, 43)
    conditional_likelihood_test(db, plot=True)


    db.commit_and_close()


main()

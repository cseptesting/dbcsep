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


from dbcsep import SpatiaLiteCSEPDatabase


def main():
    """ """

    db = SpatiaLiteCSEPDatabase("test.db")
    db.connect()
    db.create_tables()
    db.create_grid(-126, -112, 32, 43, 12, description="test-grid")
    # db.download_catalog_from_fdsn_service("USGS", "2024-01-01", "2024-03-01", 0, -126, -111, 30, 43)
    # db.update_catalog_from_fdsn_service()
    db.commit_and_close()


main()

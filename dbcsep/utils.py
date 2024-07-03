#!/usr/bin/env python3

# Copyright (C) 2024:
#   Swiss Seismological Service, ETH Zurich
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


def check_is_number_between(
    number: float, min_value: float, max_value: float, inclusive="both"
) -> bool:
    """
    Checks whether the given number `number` is within the given interval specified by a minimum
    and a maximum value. The parameter `inclusive` defines which of the interval limits are
    inclusive (possible values are: `both`, `neither`, `left`, right`). In case of any error,
    the function raises an exception. If all checks are successful, the function returns `True`.

    Args:
        number (float or int):
            Number to be checked if it is a floating point or integer value and within the
            given interval.
        min_value (float or int):
            Minimum value of the interval.
        max_value (float or int):
            Maximum value of the interval.
        inclusive (`both`, `neither`, `left`, `right`):
            Defines which interval limit is considered inclusive.

    Returns:
        True.
    """

    # Check if the `inclusive` parameter is correct.
    if inclusive not in ["both", "neither", "left", "right"]:
        raise ValueError(
            "Parameter `inclusive` can be of `both`, `neither`, `left`, `right`"
            f" but `{inclusive}` was given."
        )

    # Check for lower bound of interval.
    if number < min_value:
        raise ValueError(f"Number `{number}` is smaller than minimum value `{min_value}`.")
    if inclusive in ["neither", "right"] and number == min_value:
        raise ValueError(f"Number `{number}` is equal to minimum value `{min_value}`.")

    # Check for upper bound of interval.
    if number > max_value:
        raise ValueError(f"Number `{number}` is larger than maximum value `{max_value}`.")
    if inclusive in ["neither", "left"] and number == max_value:
        raise ValueError(f"Number `{number}` is equal to maximum value `{max_value}`.")

    return True

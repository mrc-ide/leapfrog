from typing import Any

import numpy as np


def assert_keys(data: dict[str, Any], keys: list[str]):
    missing_keys = [key for key in keys if key not in data.keys()]
    if missing_keys:
        msg = f"Missing required keys: {', '.join(missing_keys)}"
        raise ValueError(msg)


def assert_fortran_order(data: dict[str, np.ndarray]):
    for key, value in data.items():
        if value.flags["F_CONTIGUOUS"] is False:
            msg = f"Array '{key}' is not in Fortran order"
            raise ValueError(msg)

import numpy as np


def assert_keys(data: dict[str, any], keys: list[str]):
    missing_keys = [key for key in keys if key not in data.keys()]
    if missing_keys:
        raise ValueError(f"Missing required keys: {', '.join(missing_keys)}")


def assert_fortran_order(data: dict[str, np.ndarray]):
    for key, value in data.items():
        if value.flags['F_CONTIGUOUS'] is False:
            raise ValueError(f"Array '{key}' is not in Fortran order")

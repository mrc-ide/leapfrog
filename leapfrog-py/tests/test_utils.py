import numpy as np
import pytest

from leapfrog_py.utils import assert_fortran_order, assert_keys


def test_assert_keys():
    test_input = {"basepop": np.arange(10)}
    with pytest.raises(
        ValueError,
        match="Missing required keys: survival_probability, net_migration",
    ):
        assert_keys(
            test_input, ["basepop", "survival_probability", "net_migration"]
        )


def test_assert_fortran_order():
    test_input = {"basepop": np.zeros(shape=(2, 3))}
    with pytest.raises(
        ValueError, match="Array 'basepop' is not in Fortran order"
    ):
        assert_fortran_order(test_input)

    # No error raised
    test_input_f = {
        "basepop": np.zeros(shape=(2, 3), order="F"),
        "1d_array": np.arange(10),
    }
    assert_fortran_order(test_input_f)

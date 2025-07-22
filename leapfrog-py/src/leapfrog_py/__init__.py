from typing import Union, cast
from h5py._hl.dataset import Dataset
import numpy as np
import h5py

from leapfrog_py._core import run_base_model, run_base_model_from_state, run_base_model_single_year


def read_h5_file(file_path: str) -> dict[str, Union[np.ndarray, int, float]]:
    ret = {}
    f = h5py.File(file_path, "r")
    for key in list(f.keys()):
        dset = cast(Dataset, f[key])
        if dset.size == 1:
            val = dset[0]
            if isinstance(val, bytes):
                ret[key] = val.decode("utf-8")
            else:
                ret[key] = val
        else:
            ret[key] = np.transpose(np.array(f[key]))
    f.close()
    return ret


def save_h5_file(dat: dict[str, Union[np.ndarray, int, float]], file_path: str) -> None:
    f = h5py.File(file_path, "w")
    for key in list(dat.keys()):
        arr = dat[key]
        if isinstance(arr, np.ndarray):
            f.create_dataset(key, data=np.transpose(arr))
        else:
            f[key] = dat[key]
    f.close()

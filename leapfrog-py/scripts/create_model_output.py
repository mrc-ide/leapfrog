#!/usr/bin/env python3

"""Run leapfrog model and save output to specified dir
Usage:
  run_model <output-dir>

Arguments:
  <output-dir>  Path to save output to.

Options:
  -h --help                  Show this screen.
"""

from docopt import docopt
import os

from leapfrog_py import read_h5_file, run_model, save_h5_file

if __name__ == '__main__':
    args = docopt(__doc__)
    output_dir = args["<output-dir>"]
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


    parameters = read_h5_file("../r-package/tests/testthat/testdata/adult_parms_full.h5")
    ret = run_model(parameters)

    save_h5_file(ret, os.path.join(output_dir, "py-output.h5"))

# src/radar_lwc/preprocessing/radar.py

import numpy as np
import sys
from pathlib import Path
sys.path.append(str(Path("../src").resolve()))

from radar_lwc.preprocessing.chirps import combine_chirps

def get_reflectivity(ds):
    """
    Return combined reflectivity in dBZ.
    """

    ze = combine_chirps(ds, "ZE")

    ze = ze.where(ze != -999)

    ze = 10 * np.log10(ze)

    return ze
"""
NIRCam Coronagraphy
"""
import re

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

from whippot import whippot_tools
from whippot import whippot_plots

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):

    def _prefilter_apertures(self, aperture_list) -> list:
        """
        Function that returns a filtered list of apertures that can be selected
        Written this way so that it can be overriden by subclasses
        """
        # NRCA only 
        apernames = list(filter(lambda i: i.startswith("NRCA"), aperture_list))
        # just the simple apertures
        apernames = list(filter(lambda i: i.count("_") == 2, aperture_list))
        # Bar or Round masks
        apernames = list(filter(lambda i: ( i.endswith('R') or i.endswith("WB") ), aperture_list))
        return list(apernames)


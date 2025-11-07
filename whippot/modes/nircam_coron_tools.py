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
        # apernames = [i for i in  aperture_list if 'MASK' in i.upper()]
        apernames = [i for i in aperture_list if (i.endswith('R') or i.endswith("WB"))]
        return list(apernames)


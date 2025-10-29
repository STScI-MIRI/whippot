"""
Miscellaneous plot methods
"""
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

from whippot import whippot_tools
from whippot import whippot_plots

# trace size reference: Andreea Petric, personal communication
aper = whippot_tools.Siaf("MIRI")['MIRIM_FULL']
TRACE_UP = 100 * aper.YSciScale
TRACE_DOWN = 300 * aper.YSciScale

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):

    def aperture_filter(self):
        """
        Function that returns a filtered list of apertures that can be selected
        Written this way so that it can be overriden by subclasses
        """
        apernames = [i for i in  whippot_tools.Siaf(self.parameter_values['instr']).apernames if 'MIRIFU_CHANNEL' in i]
        return list(apernames)

    def plot_scene(self, *args) -> mpl.figure.Figure:
        # copy the docstring
        super().plot_scene.__doc__

        fig = super().plot_scene(self, *args)
        # get the IDL and Sky axes
        idl_ax, sky_ax = fig.get_axes()

        # for the chosen channel, overlay the slit footprints on the axes
        channelid = self.aperture.AperName[-2:]
        slit_apernames = [i for i in self.instr.apernames if f'MIRIFU_{channelid}' in i]
        for apername in slit_apernames:
            new_aper = self.instr[apername]
            # idl
            footprint = whippot_tools.transform_aper_footprint(
                new_aper, self.aperture, 'idl',
                label=apername, zorder=-1
            )
            idl_ax.add_patch(footprint)
            # sky
            footprint = whippot_tools.transform_aper_footprint(
                new_aper, self.aperture, 'sky',
                label=apername, zorder=-1
            )
            sky_ax.add_patch(footprint)

        return fig

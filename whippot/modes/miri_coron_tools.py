"""
MIRI Coronagraphy
"""
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

from whippot import whippot_tools, whippot_plots

aper = whippot_tools.Siaf("MIRI")['MIRIM_ILLUM']
TRACE_UP = 100 * aper.YSciScale
TRACE_DOWN = 300 * aper.YSciScale

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):
    def plot_scene(self, *args) -> mpl.figure.Figure:
        # copy the docstring
        super().plot_scene.__doc__

        fig = super().plot_scene(self, *args)
        if self._SELF_TA:
            idl_ax, sky_ax = fig.get_axes()
            idl_axes, sky_axes = [idl_ax], [sky_ax]
        else:
            axes = fig.get_axes()
            idl_axes, sky_axes = axes[::2], axes[1::2]

        # show the ILLUM and FULL apertures
        # also show the SLITLESSUPPER and LOWER apertures
        coron_id = self.aperture.AperName[-4:]
        for apername in ['MIRIM_CORON'+coron_id]:
            for idl_ax in idl_axes:
                new_aper = self.instr[apername]
                footprint = whippot_plots.transform_aper_footprint(new_aper, self.aperture, 'idl', label=apername)
                idl_ax.add_patch(footprint)
            for sky_ax in sky_axes:
                footprint = whippot_plots.transform_aper_footprint(new_aper, self.aperture, 'sky', label=apername)
                sky_ax.add_patch(footprint)

        return fig

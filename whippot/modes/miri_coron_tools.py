"""
Miscellaneous plot methods
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
        idl_ax, sky_ax = fig.get_axes()

        # show the ILLUM and FULL apertures
        # also show the SLITLESSUPPER and LOWER apertures
        which = self.aperture.AperName[-4:]
        for apername in ['MIRIM_MASK'+which, 'MIRIM_CORON'+which]:
            new_aper = self.instr[apername]
            footprint = whippot_plots.transform_aper_footprint(new_aper, self.aperture, 'idl', label=apername)
            idl_ax.add_patch(footprint)
            footprint = whippot_plots.transform_aper_footprint(new_aper, self.aperture, 'sky', label=apername)
            sky_ax.add_patch(footprint)

        whippot_plots.include_patches_in_axes(idl_ax, idl_traces, invert_ra_axis=False)
        whippot_plots.include_patches_in_axes(sky_ax, sky_traces, invert_ra_axis=True)
        # add the sky trace to the sky axes
        # add the idl trace to the idl axis
        for it, st in zip(idl_traces, sky_traces):
            idl_ax.add_patch(it)
            sky_ax.add_patch(st)

        return fig

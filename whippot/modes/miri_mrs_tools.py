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
        idl_ax, sky_ax = fig.get_axes()
        # idl_axes, sky_axes = fig.get_axes()

        # for the chosen channel, show the slits
        channelid = self.aperture.AperName[-2:]
        slit_apernames = [i for i in self.instr.apernames if f'MIRIFU_{channelid}' in i]
        print(channelid, slit_apernames)
        for apername in slit_apernames:
            new_aper = self.instr[apername]
            # footprint = whippot_tools.transform_aper_footprint(new_aper, self.aperture, 'idl', label=apername)
            # idl_ax.add_patch(footprint)
            # footprint = whippot_tools.transform_aper_footprint(new_aper, self.aperture, 'sky', label=apername)
            # sky_ax.add_patch(footprint)

    #     trace_properties = dict(
    #         facecolor='C0', alpha=0.5, zorder=-1, linestyle='none'
    #     )

    #     # for each source, add its trace to the detector and sky axes
    #     idl_traces, sky_traces = [], []
    #     for i, (k, coord) in enumerate(self.idl_coords_after_slew.items()):
    #         # plot each trace as a Rectangle, defining the height, width, and bottom corner
    #         height, width = (TRACE_UP + TRACE_DOWN), 1
    #         # ll -> lower left corner
    #         ll = (coord[0]-width/2, coord[1]-TRACE_DOWN)
    #         idl_traces.append(mpl.patches.Rectangle(ll, width, height, **trace_properties))

    #         # transform the trace idl vertices to sky coordinates
    #         sky_traces.append(whippot_tools.transform_patch_footprint(
    #             idl_traces[-1], self.aperture, 'idl', 'sky', **trace_properties
    #         ))
    #     # for some reason you have to find the axis limits *before* you add the patches to the plot,
    #     # perhaps because the act of adding them changes their vertices
    #     whippot_plots.include_patches_in_axes(idl_ax, idl_traces, invert_ra_axis=False)
    #     whippot_plots.include_patches_in_axes(sky_ax, sky_traces, invert_ra_axis=True)
    #     # add the sky trace to the sky axes
    #     # add the idl trace to the idl axis
    #     for it, st in zip(idl_traces, sky_traces):
    #         idl_ax.add_patch(it)
    #         sky_ax.add_patch(st)
    #     return fig

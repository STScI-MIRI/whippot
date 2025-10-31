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


def plot_wfss_traces(
    cp : whippot_tools.ComputePositions,
    ax : mpl.axes.Axes = None,
    title : str = '',
    plot_full : bool = False,
    show_mirim_illum : bool = True
)-> mpl.figure.Figure:
    """
    Compute the idl coordinates overlaid on the aperture of interest when the
    SCI target is at the reference position. Add the WFSS traces as well.

    Parameters
    ----------
    cp : a whippot_tools.ComputePositions instance
    aper : usually the return value of utils.ComputeOffsets().get_aper()
    ax : axis to plot on
    title : title for the axis
    show_aper_title : if True, display the text MIRIM_ILLUM
    plot_full : Show (True) or hide (False) the FULL detector array boundaries
    show_mirim_illum : Show (True) or hide (False) the words "MIRIM_ILLUM" at the center of the aperture
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    ax.set_title(title)

    aper = cp.get_aper()
    # trace size reference: Andreea Petric, personal communication
   
    for i, (k, coord) in enumerate(cp.idl_coords_after_slew.items()):
        # plot each trace as a Rectangle, defining the height, width, and bottom corner
        height, width = (TRACE_UP + TRACE_DOWN), 1
        # ll -> lower left corner
        ll = (coord[0]-width/2, coord[1]-TRACE_DOWN)
        trace = mpl.patches.Rectangle(ll, width, height, facecolor=f'C0', alpha=0.5)
        ax.add_patch(trace)
        ax.scatter(*coord, c=f'C{i}', marker='x',label=k)
    aper.plot(ax=ax, fill=False, mark_ref=False, frame='idl', label=show_mirim_illum, zorder=-1)

    if plot_full:
        full_aper = whippot_tools.Siaf("MIRI")['MIRIM_FULL']
        # conver the corners to ILLUM IDL
        corners = aper.convert(*full_aper.corners("tel"), from_frame="tel", to_frame='idl')
        full_rect = mpl.patches.Rectangle(
            xy = np.min(corners, axis=1),
            width = np.diff(corners[0]).max(),
            height = np.diff(corners[1]).max(),
            fill=False, ec='gray'
        )
        ax.add_patch(full_rect)
    ax.legend(loc=(1.05, 0.3), title='Sources')
    return fig

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):
    def plot_scene(self, *args) -> mpl.figure.Figure:
        # copy the docstring
        super().plot_scene.__doc__

        fig = super().plot_scene(self, *args)
        idl_ax, sky_ax = fig.get_axes()

        # for each source, add its trace to the detector and sky axes
        idl_traces, sky_traces = [], []
        trace_properties = dict(
            facecolor='C0', alpha=0.5, zorder=-1, linestyle='none'
        )
        for i, (k, coord) in enumerate(self.idl_coords_after_slew.items()):
            # plot each trace as a Rectangle, defining the height, width, and bottom corner
            height, width = (TRACE_UP + TRACE_DOWN), 1
            # ll -> lower left corner
            ll = (coord[0]-width/2, coord[1]-TRACE_DOWN)
            idl_traces.append(mpl.patches.Rectangle(ll, width, height, **trace_properties))

            # transform the trace idl vertices to sky coordinates
            sky_traces.append(whippot_plots.transform_patch_footprint(
                idl_traces[-1], self.aperture, 'idl', 'sky', **trace_properties
            ))
        # for some reason you have to find the axis limits *before* you add the patches to the plot,
        # perhaps because the act of adding them changes their vertices
        whippot_plots.include_patches_in_axes(idl_ax, idl_traces, invert_ra_axis=False)
        whippot_plots.include_patches_in_axes(sky_ax, sky_traces, invert_ra_axis=True)
        # add the sky trace to the sky axes
        # add the idl trace to the idl axis
        for it, st in zip(idl_traces, sky_traces):
            idl_ax.add_patch(it)
            sky_ax.add_patch(st)

        return fig

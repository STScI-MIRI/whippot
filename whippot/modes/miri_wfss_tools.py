"""
MIRI Wide-Field Slitless Spectroscopy
"""
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

from whippot import whippot_tools, whippot_plots

aper = whippot_tools.Siaf("MIRI")['MIRIM_ILLUM']
TRACE_UP = 100 * aper.YSciScale
TRACE_DOWN = 300 * aper.YSciScale

# shared plotting properties for traces
trace_properties = dict(
    facecolor='C0', alpha=0.5, zorder=-1, linestyle='none'
)

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):

    def make_trace_patch(
        self, idl_coord, frame='idl', **trace_properties
    ) -> mpl.patches.PathPatch:
        """
        Parametrize the trace in the IDL coordinate frame to compute where it
        falls on the detector. Then, use the points to make a matplotlib
        PathPatch that can be transformed to the coordinate system specified by
        `frame`. Return the transformed patch.
        """
        height, width = (TRACE_UP + TRACE_DOWN), 1
        # ll -> lower left corner
        ll = (idl_coord[0]-width/2, idl_coord[1]-TRACE_DOWN)
        patch = mpl.patches.Rectangle(ll, width, height)

        # convert this to a pathpatch
        xverts, yverts = patch.get_verts().T
        codes = patch.get_path().codes
        transf_verts = np.array(self.aperture.convert(xverts, yverts, 'idl', frame)).T
        trace = mpl.patches.PathPatch(
            mpl.path.Path(list(transf_verts), codes=codes),
            **trace_properties,
        )
        return trace

    def plot_scene(self, *args) -> mpl.figure.Figure:
        # copy the docstring
        super().plot_scene.__doc__

        fig = super().plot_scene(self, *args)
        idl_ax, sky_ax = fig.get_axes()

        # show the ILLUM and FULL apertures
        # also show the SLITLESSUPPER and LOWER apertures
        for apername in ['MIRIM_FULL', 'MIRIM_ILLUM']:
            new_aper = self.instr[apername]
            footprint = whippot_plots.transform_aper_footprint(new_aper, self.aperture, to_frame=frame, label=apername)
            whippot_plots.include_patches_in_axes(idl_ax, [footprint], invert_ra_axis=False)
            idl_ax.add_patch(footprint)
            footprint = whippot_plots.transform_aper_footprint(new_aper, self.aperture, to_frame='sky', label=apername)
            whippot_plots.include_patches_in_axes(sky_ax, [footprint], invert_ra_axis=True)
            sky_ax.add_patch(footprint)

        # for each source, add its trace to the detector and sky axes
        idl_traces, sky_traces = [], []
        for i, (k, coord) in enumerate(self.idl_coords_after_slew.items()):
            idl_traces.append(self.make_trace_patch(coord, 'idl', **trace_properties))
            sky_traces.append(self.make_trace_patch(coord, 'sky', **trace_properties))

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

def plot_on_data(
        cp : ComputePositions,
        img : np.ndarray,
        hdr : dict | None = None,
) -> mpl.figure.Figure:
    """
    Overlay the scene on provided data. You have to be careful here about handling the aperture name

    Parameters
    ----------
    cp : ComputePositions
    img : np.ndarray
      A 2-D exposure that corresponds to the SUBARRAY keyword from the data header
    hdr : dict = None
      the PRI header corresponding to the data. Can be NONE if the APERNAME keyword corresponds to a subarray

    Output
    ------
    fig : mpl.figure.Figure hande

    """
    # get an aperture defined on the detector
    idl_pos = cp.idl_coords_after_slew
    # if your aperture can go straight between idl and sci, great
    # if not, go through tel
    aper_prefix = cp.aperture.AperName.split("_")[0]
    subarray_aper = whippot_tools.Siaf(hdr['INSTRUME'])[aper_prefix + "_" + hdr['SUBARRAY']]
    try:
        det_pos = {k: cp.aperture.idl_to_sci(*v) for k, v in idl_pos.items()}
    except TypeError:
        tel_pos = {k: cp.aperture.idl_to_tel(*v) for k, v in idl_pos.items()}
        det_pos = {k: subarray_aper.tel_to_sci(*v) for k, v in tel_pos.items()}
    dims = img.shape
    xcoords = np.arange(dims[1]+1) + 0.5
    ycoords = np.arange(dims[0]+1) + 0.5
    fig, ax = plt.subplots(1, 1)
    vmin, vmax = np.nanquantile(img, [0.01, 0.99])
    ax.pcolormesh(
        xcoords, ycoords, img,
        vmin=vmin, vmax=vmax,
        cmap=mpl.cm.cividis,
    )
    ax.set_aspect("equal")
    for label, pos in det_pos.items():
        ax.scatter(*pos, label=label, marker='x', s=100)
    ax.legend(loc=(1.1, 0.5))

    traces = []
    properties = trace_properties.copy()
    properties.update({'zorder': 10})
    for i, (k, coord) in enumerate(det_pos.items()):
        # plot each trace as a Rectangle, defining the height, width, and bottom corner
        height, width = (TRACE_UP + TRACE_DOWN) / subarray_aper.YSciScale, 1 / subarray_aper.YSciScale
        # ll -> lower left corner
        ll = (coord[0]-width/2, coord[1]-TRACE_DOWN/subarray_aper.YSciScale)
        idl_trace = mpl.patches.Rectangle(ll, width, height, **properties)
        traces.append(idl_trace)
    # for some reason you have to find the axis limits *before* you add the patches to the plot,
    # perhaps because the act of adding them changes their vertices
    whippot_plots.include_patches_in_axes(ax, traces, invert_ra_axis=False)
    for t in traces:
        ax.add_patch(t)

    return fig



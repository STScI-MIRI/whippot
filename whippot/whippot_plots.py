"""
More advanced plotting tools for use in notebooks
Importable into different observing mode plots.
These methods can be overloaded or wrapper.
"""
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import patches, path

from astropy.coordinates import SkyCoord
from astropy import units

import pysiaf
from pysiaf import Siaf
from whippot import list_of_masks as lom


def plot_aper_to_frame(
    aper : pysiaf.aperture.JwstAperture,
    star_positions : dict[str, np.ndarray] = {},
    frame_from : str = 'idl',
    frame_to : str = 'idl',
    offset : list | np.ndarray = [0., 0.],
    ax = None,
    title : str = '',
    show_legend : bool = True,
    idl_mask = None,
    show_diffraction_spikes : bool = False,
) -> mpl.figure.Figure:
    """
    Make a plot (aperture, masks, sources, etc.) in an arbitrary telescope frame

    Parameters
    ----------
    aper : pysiaf.aperture.JwstAperture
      the aperture to plot in IDL
    star_positions : dict[str, np.ndarray] = {}
      dictionary of star coordinates in IDL
    offset : list | np.ndarray = [0., 0.]
      any additional offset to aplpy to the star positions
    ax = None
      the axes object to plot on
    title = ''
      a title for the axis
    show_legend : bool = True
      if True, show the legend
    idl_mask = None
      a matplotlib Patch object to overlay on the aperture, in IDL coordinates
    frame_from = 'idl'
    frame_to = 'idl'
      options are: tel (telescope), det (detector), sci (aperture)

    Output
    ------
    fig : a handle to the matplotlib figure
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout='constrained')
    else:
        fig = ax.get_figure()
    if title != '':
        ax.set_title(title)
    aper.plot(ax=ax, label=False, frame=frame_to, c='gray', mark_ref=True, fill=False)
    if idl_mask is not None:
        if frame_to == 'idl':
            ax.add_artist(idl_mask)
        else:
            ax.add_artist(
                transform_patch_footprint(
                    idl_mask, aper,
                    frame_from='idl', frame_to=frame_to,
                    **lom.mask_kwargs
                )
            )

    offset = np.array(offset)
    # star_positions = star_positions.copy()
    star_positions = {
        label: aper.convert(*(np.array(pos)+np.array(offset)), frame_from, frame_to)
        for label, pos in star_positions.items()
    }
    if 'ACQ' in star_positions.keys():
        acq_pos = star_positions['ACQ']
        ax.scatter(*acq_pos,
                   c='k',
                   label=f"ACQ",
                   marker='x',
                   s=100)
    if 'SCI' in star_positions.keys():
        sci_pos = star_positions["SCI"]
        ax.scatter(*sci_pos,
                   label=f"SCI",
                   marker="*",
                   s=100,
                   c='k')
    for star, position in star_positions.items():
        if star in ['SCI','ACQ']:
            continue
        ax.scatter(*position,
                   # c='k',
                   label=star,
                   marker='.',
                   s=100)

    if show_diffraction_spikes:
        draw_diffraction_spikes(
            ax,
            aper,
            star_positions,
            source_frame=frame_to
        )

    if show_legend:
        ax.legend(loc=(1.05, 0.3))
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    if frame_to == 'sky':
        ax.invert_xaxis()
    return fig


def include_patches_in_axes(ax, list_of_patches, invert_ra_axis=False) -> None:
    """
    Adjust axis limits to include patches

    invert_ra_axis : False
      if plotting RA on the x-axis, invert it
    """
    (xmin, xmax), (ymin, ymax) = ax.get_xlim(), ax.get_ylim()
    if invert_ra_axis:
        xmin, xmax = xmax, xmin

    verts = np.concatenate([p.get_verts() for p in list_of_patches])
    vxmin, vymin = np.min(verts, axis=0)
    vxmax, vymax = np.max(verts, axis=0)
    xmin = np.min([xmin, vxmin])
    ymin = np.min([ymin, vymin])
    xmax = np.max([xmax, vxmax])
    ymax = np.max([ymax, vymax])

    # add a 10% buffer on either side of the mins and maxes
    dx = np.abs(xmax-xmin)
    dy = np.abs(ymax-ymin)
    ax.set_ylim(ymin - 0.1*dy, ymax + 0.1*dy)
    ax.set_xlim(xmin - 0.1*dx, xmax + 0.1*dy)

    if invert_ra_axis:
        ax.invert_xaxis()

    return 

def draw_diffraction_spikes(
    ax,
    aperture,
    sources={},
    source_frame = 'idl',
    show_inner_diff_spikes : bool = False,
) -> None:
    """
    Draw diffraction spikes on the sources at the given position
    Code written originally by M. Perrin

    parameters
    ----------
    source_frame : the frame in which the source coordinates are defined
    """
    spike_length = 4 # arcsec
    # if source_frame == 'sky':
    #     spike_length = 4 * units.arcsec.to(units.degree)
    #     print(spike_length)
    outer_spikelen = spike_length
    inner_spikelen = 0.5
    for angle in range(6):
        # Big outer diffraction spikes, from the individual hexagons
        for label, pos in sources.items():
            # compute the spike positions in TEL frame
            tel_pos = aperture.convert(*pos, source_frame, 'tel')
            ang_rad = np.deg2rad(angle * 60)
            spike_pts = (
                tel_pos[0] - np.asarray([inner_spikelen, outer_spikelen]) * np.sin(ang_rad),
                tel_pos[1] + np.asarray([inner_spikelen, outer_spikelen]) * np.cos(ang_rad)
            )
            # convert back to  the source frame
            spike_pts = aperture.convert(*spike_pts, 'tel', source_frame)
            ax.plot(
                # pos[0] - np.asarray([inner_spikelen, outer_spikelen]) * np.sin(ang_rad),
                # pos[1] + np.asarray([inner_spikelen, outer_spikelen]) * np.cos(ang_rad),
                *spike_pts,
                color='red', lw=1, marker='none',
            )
    #         # smaller inner diffraction spikes, from overall primary
    #         if show_inner_diff_spikes:
    #             ang_rad = np.deg2rad(angle * 60 + 30 + aperture.V3IdlYAngle)
    #             ax.plot(
    #                 [pos[0], pos[0] - np.sin(ang_rad) * inner_spikelen],
    #                 [pos[1], pos[1] + np.cos(ang_rad) * inner_spikelen],
    #                 color='black', lw=2, marker='none', zorder=-1, alpha=0.5,
    #             )
    #         # Extra horizontal spikes from the +V3 SM strut
    # for angle in range(2):
    #     ang_rad = np.deg2rad(angle * 180 + 90 + aperture.V3IdlYAngle)
    #     spikelen = 1.5
    #     for label, pos in sources.items():
    #         ax.plot(
    #             [pos[0], pos[0] - np.sin(ang_rad) * outer_spikelen],
    #             [pos[1], pos[1] + np.cos(ang_rad) * outer_spikelen],
    #             color='black', lw=1, marker='none'
    #         )
    return


def transform_aper_footprint(
    aper_from,
    aper_to,
    to_frame='idl',
    **patch_kwargs,
) -> mpl.patches.Patch:
    """
    Translate an aperture footprint to the IDL, SCI, or DET frame of another aperture
    """
    vertices = aper_from.closed_polygon_points("tel", rederive=False)
    vertices  = list(zip(*aper_to.convert(
        *vertices,
        from_frame='tel',
        to_frame=to_frame
    )))
    default_kwargs = dict(fill=False, ec='gray')
    default_kwargs.update(patch_kwargs)
    footprint = mpl.patches.PathPatch(
        mpl.patches.Path(vertices=vertices, closed=True),
        **default_kwargs
    )
    return footprint

def transform_patch_footprint(
    patch : mpl.patches.Patch,
    aper : pysiaf.aperture.JwstAperture,
    frame_from : str = 'idl',
    frame_to : str = 'sky',
    **patch_kwargs,
) -> mpl.patches.Patch :
    """
    Take a patch and transform its vertices to a different coordinate system

    Parameters
    ----------
    patch : mpl.patches.Patch
    aper : pysiaf.aperture.JwstAperture
    frame_from : str = 'idl'
    frame_to : str = 'sky'

    Output
    ------
    transf_patch : mpl.patches.Patch

    """
    p = patch.get_path()
    codes = patch.get_path().codes
    x, y = patch.get_verts().T
    transf_verts = np.array(aper.convert(x, y, frame_from, frame_to)).T
    transf_patch = patches.PathPatch(
        path.Path(list(transf_verts), codes=codes),
        **patch_kwargs,
    )
    return transf_patch

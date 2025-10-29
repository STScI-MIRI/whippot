"""
More advanced plotting tools for use in notebooks
Importable into different observing mode plots.
These methods can be overloaded or wrapper.
"""
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon

from astropy.coordinates import SkyCoord
from astropy import units

import pysiaf
from pysiaf import Siaf
from whippot import list_of_masks as lom

def plot_aper_idl(
        aper : pysiaf.aperture.JwstAperture,
        star_positions : dict[str, np.ndarray] = {},
        offset : list | np.ndarray = [0., 0.],
        ax = None,
        title = '',
        show_legend : bool = True,
        idl_mask = None
):
    """
    Plot the scene on the detector when you're pointed at the science target
    
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
    """
    # plot 1 : POV of the detector
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout='constrained')
    else:
        fig = ax.get_figure()
    frame = 'idl' # options are: tel (telescope), det (detector), sci (aperture)
    if title != '':
        ax.set_title(title)
    aper.plot(ax=ax, label=False, frame=frame, c='gray', mark_ref=True)
    if idl_mask is not None:
        ax.add_artist(idl_mask)

    offset = np.array(offset)
    star_positions = star_positions.copy()
    if 'ACQ' in star_positions.keys():
        acq_pos = star_positions.pop('ACQ')
        ax.scatter(acq_pos[0] + offset[0], acq_pos[1] + offset[1],
                   c='k',
                   label=f"ACQ",
                   marker='x',
                   s=100)
    if 'SCI' in star_positions.keys():
        sci_pos = star_positions.pop("SCI")
        ax.scatter(*(sci_pos+ offset),
                   label=f"SCI",
                   marker="*",
                   c='k')
    for star, position in star_positions.items():
        ax.scatter(*(position + offset),
                   # c='k',
                   label=star,
                   marker='.',
                   s=50)

    if show_legend:
        ax.legend(loc=(1.05, 0.3))
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    return fig


def plot_aper_sky(
        aper : pysiaf.aperture.JwstAperture,
        star_positions : dict[str, np.ndarray] = {},
        offset : list | np.ndarray = [0., 0.],
        ax = None,
        title = '',
        show_legend : bool = True,
        idl_mask = None,
):
    """
    Plot the scene on the sky when you're pointed at the science target
    
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
    """
    # plot 1 : POV of the detector
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout='constrained')
    else:
        fig = ax.get_figure()
    frame = 'sky' # options are: tel (telescope), det (detector), sci (aperture)
    if title != '':
        ax.set_title(title)
    aper.plot(ax=ax, label=False, frame=frame, c='gray', mark_ref=True)
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")

    if idl_mask is not None:
        vertices = np.array(aper.idl_to_sky(*idl_mask.get_verts().T)).T
        properties = idl_mask.properties()
        properties.pop('path')
        sky_mask = mpl.patches.PathPatch(mpl.path.Path(vertices, closed=True))
        sky_mask.set(**lom.mask_kwargs)
        ax.add_artist(sky_mask)

    offset = np.array(offset)
    star_positions = star_positions.copy()
    if 'ACQ' in star_positions.keys():
        acq_pos = star_positions.pop('ACQ')
        acq_pos = np.array(acq_pos) + np.array(offset)
        acq_pos = aper.idl_to_sky(*acq_pos)
        ax.scatter(*acq_pos,
                   c='k',
                   label=f"ACQ",
                   marker='x',
                   s=100)
    if 'SCI' in star_positions.keys():
        sci_pos = star_positions.pop("SCI")
        sci_pos = np.array(sci_pos) + np.array(offset)
        sci_pos = aper.idl_to_sky(*sci_pos)
        ax.scatter(*sci_pos,
                   label=f"SCI",
                   marker="*",
                   c='k')
    for star, position in star_positions.items():
        position = np.array(position) + np.array(offset)
        position = aper.idl_to_sky(*position)
        ax.scatter(*position,
                   # c='k',
                   label=star,
                   marker='.',
                   s=50)

    if show_legend:
        ax.legend(loc=(1.05, 0.3))
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
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

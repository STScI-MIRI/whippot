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
    aper.plot(ax=ax, label=False, frame=frame, c='black', mark_ref=True)
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
        print(sky_mask.get_verts())
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

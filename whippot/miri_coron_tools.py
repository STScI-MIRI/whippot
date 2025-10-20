"""
Custom tools for the MIRI coronagraphs

The main thing we want to add here is to show the 4QPM and Lyot masks on the plots
"""
import itertools
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import patches 

from whippot import whippot_tools
from whippot import whippot_plots


def lyot_fixtures(aperture, kwargs={}):
    """Generate a polygon to show the Lyot mask and support bar"""

    y_angle = np.deg2rad(aperture.V3IdlYAngle)
    corners_x, corners_y = aperture.corners(to_frame='idl')
    min_x, min_y = np.min(corners_x), np.min(corners_y)
    max_x, max_y = np.max(corners_x), np.max(corners_y)

    width = 0.7 # bar width in arcsec
    # draw a rectangle and then rotate it around the fiducial point
    ll = [-width/2, aperture.corners('idl')[1].min()]
    height = aperture.YSciSize * aperture.YSciScale
    bar = patches.Rectangle(ll, width=width, height=height, rotation_point='center', angle=-aperture.V3IdlYAngle, **kwargs)
    spot = patches.Circle((0, 0), radius=2.16, **kwargs)
    # verts = np.concatenate([x_verts[:, np.newaxis], y_verts[:, np.newaxis]], axis=1)
    # mask = Polygon(verts, alpha=0.5, **kwargs)
    mask = mpl.collections.PatchCollection([bar, spot], **kwargs)
    return mask

def quad_boundaries(aperture, kwargs={}):
    """
    Generate a polygon to plot the 4QPM quadrant boundaries. Stolen from the JWST
    coronagraphic visibility tool

    Parameters
    ----------
    aperture: a pysiaf.Siaf aperture for the 1065, 1140, or 1550 coronagraph
    kwargs : {} arguments to pass to Polygon

    Output
    ------
    mask : matplotlib.patches.Polygon object
    """

    y_angle = np.deg2rad(aperture.V3IdlYAngle)
    corners_x, corners_y = aperture.corners(to_frame='idl')
    min_x, min_y = np.min(corners_x), np.min(corners_y)
    max_x, max_y = np.max(corners_x), np.max(corners_y)

    width_arcsec = 0.33
    x_verts0 = np.array([
        min_x,
        -width_arcsec,
        -width_arcsec,
        width_arcsec,
        width_arcsec,
        max_x,
        max_x,
        width_arcsec,
        width_arcsec,
        -width_arcsec,
        -width_arcsec,
        min_x
    ])
    y_verts0 = np.array([
        width_arcsec,
        width_arcsec,
        max_y,
        max_y,
        width_arcsec,
        width_arcsec,
        -width_arcsec,
        -width_arcsec,
        min_y,
        min_y,
        -width_arcsec,
        -width_arcsec
    ])
    x_verts = np.cos(y_angle) * x_verts0 + np.sin(y_angle) * y_verts0
    y_verts = -np.sin(y_angle) * x_verts0 + np.cos(y_angle) * y_verts0

    verts = np.concatenate([x_verts[:, np.newaxis], y_verts[:, np.newaxis]], axis=1)
    mask = patches.Polygon(verts, alpha=0.5, **kwargs)
    return mask

# wrap whippot_plots.plot_aper_idl and plot_aper_sky to show the 4QPM and Lyot masks
# doing it this way avoids having to copy over all the parameters
def mask_decorator(func):
    def wrapper(*args, **kwargs):
        aper = args[0]
        if aper.AperName[-4:] in ['1065', '1140', '1550']:
            mask = quad_boundaries(aper, kwargs={'fc': 'grey'})
        elif aper.AperName[-4:] == 'LYOT':
            mask = lyot_fixtures(aper, kwargs={'fc': 'grey'})
        else:
            mask = None
        kwargs.update({"idl_mask": mask})
        fig = whippot_plots.plot_aper_idl(*args, **kwargs)
        return fig
    return wrapper
# make a new plotting function with the wrapper
@mask_decorator
def plot_aper_idl(*args, **kwargs):
    pass

# same as the built-in object _plot_scene() function from whippot_tools,
# but with masks generated for the apertures
def plot_scene(cp : whippot_tools.ComputePositions, *args) -> mpl.figure.Figure:
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10), layout='constrained')
    fig = plot_aper_idl(
        cp.get_aper(),
        cp.idl_coords_after_ta,
        ax = axes[0, 0],
        title='Before slew (IDL)',
        show_legend = True,
    )
    fig = plot_aper_idl(
        cp.get_aper(),
        cp.idl_coords_after_slew,
        ax = axes[0, 1],
        title='After slew (IDL)',
        show_legend = False,
    )
    fig = whippot_plots.plot_aper_sky(
        cp.get_aper(),
        cp.idl_coords_after_ta,
        ax = axes[1, 0],
        title='Before slew (Sky)',
        show_legend = False,
    )
    fig = whippot_plots.plot_aper_sky(
        cp.get_aper(),
        cp.idl_coords_after_slew,
        ax = axes[1, 1],
        title='After slew (Sky)',
        show_legend = False,
    )
    return fig

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):
    def _plot_scene(self, *args) -> mpl.figure.Figure:
        return plot_scene(self, *args)

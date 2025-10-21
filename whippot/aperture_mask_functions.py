"""
Define matplotlib patches to represent stuff in the optical path for different apertures

File format: this file is divided into two halves. The bottom half contains
functions that define the masks. The top half contains a dictionary where those
masks are assigned.

"""
import numpy as np
from matplotlib import patches
from pysiaf import Siaf

# MIRI
def miri_lyot_mask(aperture, kwargs={}):
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

def miri_4qpm_mask(aperture, kwargs={}):
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

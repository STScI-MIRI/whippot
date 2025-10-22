"""
Define matplotlib patches to represent stuff in the optical path for different apertures

File format: this file is divided into two halves. The bottom half contains
functions that define the masks. The top half contains a dictionary where those
masks are assigned.

"""
import numpy as np
import matplotlib as mpl
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

def nrc_coron_mask(aperture, kwargs={}):
    """
    Define the masks for the NIRCam coronagraphic apertures
    "Borrowed" from the JWST Coronagraphic Visibility Tool
    """
    mask_artists = []
    aperture_name = aperture.AperName
    arcsec_per_pixel = np.average([aperture.XSciScale, aperture.YSciScale])
    x_sci_size, y_sci_size = aperture.XSciSize, aperture.YSciSize
    if aperture_name[-1] == 'R':
        if '210R' in aperture_name:
            radius_arcsec = 0.40
        elif '335R' in aperture_name:
            radius_arcsec = 0.64
        elif '430R' in aperture_name:
            radius_arcsec = 0.82
        else:
            raise RuntimeError("Invalid mask!")
        # make a circle
        mask_artists.append(patches.Circle((0, 0), radius=radius_arcsec, alpha=0.5))
    else:
        x_verts = x_sci_size / 2 * np.array([-1, 1, 1, -1])
        if 'LWB' in aperture_name:
            thin_extent_arcsec = 0.58 * (2 / 4)
            thick_extent_arcsec = 0.58 * (6 / 4)
            # x_verts *= -1  # flip LWB left to right
        elif 'SWB' in aperture_name:
            thin_extent_arcsec = 0.27 * (2 / 4)
            thick_extent_arcsec = 0.27 * (6 / 4)
        else:
            raise RuntimeError("Invalid mask!")
        
        y_verts = np.array([
            thin_extent_arcsec / arcsec_per_pixel,
            thick_extent_arcsec / arcsec_per_pixel,
            -thick_extent_arcsec / arcsec_per_pixel,
            -thin_extent_arcsec / arcsec_per_pixel
        ])
        x_idl_verts, y_idl_verts = aperture.sci_to_idl(x_verts + aperture.XSciRef, y_verts + aperture.YSciRef)
        verts = np.concatenate([x_idl_verts[:, np.newaxis], y_idl_verts[:, np.newaxis]], axis=1)
        patch = patches.Polygon(verts, alpha=0.5)
        mask_artists.append(patch)
    mask = mpl.collections.PatchCollection(mask_artists, **kwargs)
    return mask

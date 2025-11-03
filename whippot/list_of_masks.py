"""
Define matplotlib patches to represent stuff in the optical path for different apertures

This file contains a dictionary indexed by the SIAF aperture name, where the
entry is function that a matplotlib patch showing various features of the field of view.
The function takes the aperture object

"""
import numpy as np
from matplotlib import patches
from pysiaf import Siaf

from whippot import aperture_mask_functions as amf

# use this function to generate the mask

mask_kwargs = {'fc': 'grey'}

def mask_maker(mask_func, aper):
    def func():
        return mask_func(aper, kwargs=mask_kwargs)
    return func

# Format: {aperture_name: mask}
list_of_masks = {}


# MIRI coronagraphs
# use the CORON apertures for all masks
list_of_masks['MIRIM_CORONLYOT'] = mask_maker(amf.miri_lyot_mask, Siaf("MIRI")['MIRIM_CORONLYOT'])
list_of_masks['MIRIM_MASKLYOT'] = list_of_masks['MIRIM_CORONLYOT']
list_of_masks['MIRIM_FULLLYOT'] = list_of_masks['MIRIM_CORONLYOT']
list_of_masks['MIRIM_CORON1065'] = mask_maker(amf.miri_4qpm_mask, Siaf("MIRI")['MIRIM_CORON1065'])
list_of_masks['MIRIM_MASK1065'] = list_of_masks['MIRIM_CORON1065']
list_of_masks['MIRIM_FULL1065'] = list_of_masks['MIRIM_CORON1065']
list_of_masks['MIRIM_CORON1140'] = mask_maker(amf.miri_4qpm_mask, Siaf("MIRI")['MIRIM_CORON1140'])
list_of_masks['MIRIM_MASK1140'] = list_of_masks['MIRIM_CORON1140']
list_of_masks['MIRIM_FULL1140'] = list_of_masks['MIRIM_CORON1140']
list_of_masks['MIRIM_CORON1550'] = mask_maker(amf.miri_4qpm_mask, Siaf("MIRI")['MIRIM_CORON1550'])
list_of_masks['MIRIM_MASK1550'] = list_of_masks['MIRIM_CORON1550']
list_of_masks['MIRIM_FULL1550'] = list_of_masks['MIRIM_CORON1550']

# NIRCam coronagraphs
list_of_masks['NRCA2_MASK210R'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCA2_MASK210R'])
list_of_masks['NRCA5_MASK335R'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCA5_MASK335R'])
list_of_masks['NRCA5_MASK430R'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCA5_MASK430R'])
list_of_masks['NRCA4_MASKSWB'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCA4_MASKSWB'])
list_of_masks['NRCA5_MASKLWB'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCA5_MASKLWB'])
list_of_masks['NRCB1_MASK210R'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCB1_MASK210R'])
list_of_masks['NRCB5_MASK335R'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCB5_MASK335R'])
list_of_masks['NRCB5_MASK430R'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCB5_MASK430R'])
list_of_masks['NRCB3_MASKSWB'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCB3_MASKSWB'])
list_of_masks['NRCB5_MASKLWB'] = mask_maker(amf.nrc_coron_mask, Siaf("NIRcam")['NRCB5_MASKLWB'])

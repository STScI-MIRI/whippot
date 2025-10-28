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
make_mask = lambda value: None if value is None else value[0](value[1], kwargs=mask_kwargs)

# Format: {aperture_name: mask}
list_of_masks = {}


# MIRI coronagraphs
# use the CORON apertures for all masks
list_of_masks['MIRIM_CORONLYOT'] = (amf.miri_lyot_mask, Siaf("MIRI")['MIRIM_CORONLYOT'])
list_of_masks['MIRIM_MASKLYOT'] = list_of_masks['MIRIM_CORONLYOT']
list_of_masks['MIRIM_FULLLYOT'] = list_of_masks['MIRIM_CORONLYOT']
list_of_masks['MIRIM_CORON1065'] = (amf.miri_4qpm_mask, Siaf("MIRI")['MIRIM_CORON1065'])
list_of_masks['MIRIM_MASK1065'] = list_of_masks['MIRIM_CORON1065']
list_of_masks['MIRIM_FULL1065'] = list_of_masks['MIRIM_CORON1065']
list_of_masks['MIRIM_CORON1140'] = (amf.miri_4qpm_mask, Siaf("MIRI")['MIRIM_CORON1140'])
list_of_masks['MIRIM_MASK1140'] = list_of_masks['MIRIM_CORON1140']
list_of_masks['MIRIM_FULL1140'] = list_of_masks['MIRIM_CORON1140']
list_of_masks['MIRIM_CORON1550'] = (amf.miri_4qpm_mask, Siaf("MIRI")['MIRIM_CORON1550'])
list_of_masks['MIRIM_MASK1550'] = list_of_masks['MIRIM_CORON1550']
list_of_masks['MIRIM_FULL1550'] = list_of_masks['MIRIM_CORON1550']

# NIRCam coronagraphs
list_of_masks['NRCA2_MASK210R'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCA2_MASK210R'])
list_of_masks['NRCA5_MASK335R'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCA5_MASK335R'])
list_of_masks['NRCA5_MASK430R'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCA5_MASK430R'])
list_of_masks['NRCA4_MASKSWB'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCA4_MASKSWB'])
list_of_masks['NRCA5_MASKLWB'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCA5_MASKLWB'])
list_of_masks['NRCB1_MASK210R'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCB1_MASK210R'])
list_of_masks['NRCB5_MASK335R'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCB5_MASK335R'])
list_of_masks['NRCB5_MASK430R'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCB5_MASK430R'])
list_of_masks['NRCB3_MASKSWB'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCB3_MASKSWB'])
list_of_masks['NRCB5_MASKLWB'] = (amf.nrc_coron_mask, Siaf("NIRcam")['NRCB5_MASKLWB'])

"""
Define matplotlib patches to represent stuff in the optical path for different apertures

This file contains a dictionary indexed by the SIAF aperture name, where the
entry is a matplotlib patch showing various features of the field of view.

"""
import numpy as np
from matplotlib import patches
from pysiaf import Siaf

from whippot import aperture_mask_functions as amf


# Format: {aperture_name: mask}
list_of_masks = {}


# MIRI coronagraphs
list_of_masks['mirim_coronlyot'] = amf.miri_lyot_mask(Siaf("miri")["MIRIM_CORONLYOT"], kwargs={'fc': 'grey'})
list_of_masks['mirim_coron1065'] = amf.miri_4qpm_mask(Siaf("miri")["MIRIM_CORON1065"], kwargs={'fc': 'grey'})
list_of_masks['mirim_coron1140'] = amf.miri_4qpm_mask(Siaf("miri")["MIRIM_CORON1140"], kwargs={'fc': 'grey'})
list_of_masks['mirim_coron1550'] = amf.miri_4qpm_mask(Siaf("miri")["MIRIM_CORON1550"], kwargs={'fc': 'grey'})
list_of_masks['mirim_masklyot'] = list_of_masks['mirim_coronlyot']
list_of_masks['mirim_mask1065'] = list_of_masks['mirim_coron1065']
list_of_masks['mirim_mask1140'] = list_of_masks['mirim_coron1140']
list_of_masks['mirim_mask1550'] = list_of_masks['mirim_coron1550']
list_of_masks['mirim_fulllyot'] = list_of_masks['mirim_coronlyot']
list_of_masks['mirim_full1065'] = list_of_masks['mirim_coron1065']
list_of_masks['mirim_full1140'] = list_of_masks['mirim_coron1140']
list_of_masks['mirim_full1550'] = list_of_masks['mirim_coron1550']

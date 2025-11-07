"""
NIRSpec IFU
"""
import re

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

from whippot import whippot_tools
from whippot import whippot_plots

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):

    def _prefilter_apertures(self, aperture_list) -> list:
        """
        Function that returns a filtered list of apertures that can be selected
        Written this way so that it can be overriden by subclasses
        """
        apernames = [i for i in  aperture_list if 'IFU' in i]
        return list(apernames)

    def plot_scene(self, *args, frame='idl') -> mpl.figure.Figure:
        # copy the docstring
        super().plot_scene.__doc__

        fig = super().plot_scene(self, *args, frame=frame)
        # get the IDL and Sky axes
        idl_ax, sky_ax = fig.get_axes()

        # channel_apernames = [f'MIRIFU_CHANNEL{i}{band}' for i in [1, 2, 3, 4]]
        slit_apernames = sorted([i for i in self.instr.apernames if 'NRS_IFU_SLICE' in i])

        idl_patches, sky_patches = [], []
        for i, apername in enumerate(slit_apernames):
            slit_params = dict(label=apername, zorder=-1, alpha=0.3, ec='k')
            new_aper = self.instr[apername]
            # idl
            idl_footprint = whippot_plots.transform_aper_footprint(
                new_aper, self.aperture, 'idl',
                **slit_params,
            )
            # sky
            sky_footprint = whippot_plots.transform_aper_footprint(
                new_aper, self.aperture, 'sky',
                **slit_params,
            )

            idl_patches.append(idl_footprint)
            sky_patches.append(sky_footprint)

        # adjust the axes
        whippot_plots.include_patches_in_axes(idl_ax, idl_patches)
        whippot_plots.include_patches_in_axes(sky_ax, sky_patches, invert_ra_axis=True)
        # add the footprints to the axes
        for patch in idl_patches:
            idl_ax.add_patch(patch)
        for patch in sky_patches:
            sky_ax.add_patch(patch)


        return fig

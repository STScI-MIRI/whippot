"""
Miscellaneous plot methods
"""
import re

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

from whippot import whippot_tools
from whippot import whippot_plots

# trace size reference: Andreea Petric, personal communication
aper = whippot_tools.Siaf("MIRI")['MIRIM_FULL']
TRACE_UP = 100 * aper.YSciScale
TRACE_DOWN = 300 * aper.YSciScale

# Make a new class that overrides whippot_tools.ComputePositions.plot_scene()
# with the one defined above
class ComputePositions(whippot_tools.ComputePositions):

    # def filter_aperture_options(self):
    #     """
    #     Function that returns a filtered list of apertures that can be selected
    #     Written this way so that it can be overriden by subclasses
    #     """
    #     apernames = [i for i in  whippot_tools.Siaf(self.parameter_values['instr']).apernames if 'MIRIFU_CHANNEL' in i]
    #     return list(apernames)

    def plot_scene(self, *args) -> mpl.figure.Figure:
        # copy the docstring
        super().plot_scene.__doc__

        fig = super().plot_scene(self, *args, frame='idl')
        # get the IDL and Sky axes
        idl_ax, sky_ax = fig.get_axes()

        channelid, band = re.search("[1-4][A-C]", self.aperture.AperName).group()

        # overlay the channel footprints as well as the slits for the selected channel

        channel_apernames = [f'MIRIFU_CHANNEL{i}{band}' for i in [1, 2, 3, 4]]
        slit_apernames = [i for i in self.instr.apernames if f'MIRIFU_{channelid}{band}' in i]

        idl_patches, sky_patches = [], []
        for apername in slit_apernames:
            slit_params = dict(label=apername, zorder=-1, alpha=0.3, ec='k')
            new_aper = self.instr[apername]
            # idl
            footprint = whippot_plots.transform_aper_footprint(
                new_aper, self.aperture, 'idl',
                **slit_params,
            )
            idl_patches.append(footprint)
            # sky
            footprint = whippot_plots.transform_aper_footprint(
                new_aper, self.aperture, 'sky',
                **slit_params,
            )
            sky_patches.append(footprint)

        for apername in channel_apernames:
            channel_params = dict(label=apername, zorder=-1, alpha=0.2, ec='k', lw=2)
            new_aper = self.instr[apername]
            # idl
            footprint = whippot_plots.transform_aper_footprint(
                new_aper, self.aperture, 'idl',
                **channel_params,
            )
            idl_patches.append(footprint)
            # sky
            footprint = whippot_plots.transform_aper_footprint(
                new_aper, self.aperture, 'sky',
                **channel_params,
            )
            sky_patches.append(footprint)

        # adjust the axes
        whippot_plots.include_patches_in_axes(idl_ax, idl_patches)
        whippot_plots.include_patches_in_axes(sky_ax, sky_patches, invert_ra_axis=True)
        # add the footprints to the axes
        for patch in idl_patches:
            idl_ax.add_patch(patch)
        for patch in sky_patches:
            sky_ax.add_patch(patch)


        return fig

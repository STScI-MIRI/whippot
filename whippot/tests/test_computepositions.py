import pytest

from astropy import units

from whippot import whippot_tools
from whippot import list_of_masks as lom
from whippot.modes import (
    miri_lrs_slitless_tools,
    miri_wfss_tools,
    miri_mrs_tools,
    miri_coron_tools,
    nirspec_ifu_tools,
)

instr_abbr = {
    'MIR': 'MIRI',
    'NRC': 'NIRCam',
    'NRS': 'NIRSpec',
    'NIS': 'NIRISS',
}

sources = {
    'SCI': whippot_tools.SkyCoord("05 24 20.7552 -70 05 1.60", frame='icrs', unit=("hourangle","degree")),
    'a': whippot_tools.SkyCoord("05 24 26.33969 -70 05 22.3545", frame='icrs', unit=("hourangle","degree")),
    'b': whippot_tools.SkyCoord("05 24 28.67861 -70 05 24.4484", frame='icrs', unit=("hourangle","degree")),
    'c': whippot_tools.SkyCoord("05 24 36.22460 -70 05 28.1876", frame='icrs', unit=("hourangle","degree")),
    'd': whippot_tools.SkyCoord("05 24 25.608 -70 05 01.66", frame='icrs', unit=("hourangle","degree")),
}

default_init = {
    'instr': 'miri',
    'sci_aper': 'mirim_full', 
    'pa': 290.,
    'sci_ra': sources['SCI'].ra.deg, 'sci_dec': sources['SCI'].dec.deg,
    'other_stars': '\n'.join(
        f"{k}: ({v.ra.deg}, {v.dec.deg})" for k, v in sources.items() if k != 'SCI'
    ),
    'filter_apertures': False,
    'show_diffraction_spikes': False,
}
# add a multi-line string of the other stars, copied from the cell above
default_init['other_stars'] = "\n".join(f"{k}: ({v.ra.deg}, {v.dec.deg})" for k, v in sources.items())

# mode-specific initializations

# test that things run without crashing
def test_ComputePositions():
    cp = whippot_tools.ComputePositions(initial_values=default_init)
    fig = cp.plot_scene()
    whippot_tools.plt.close(fig)
    return

def test_MiriLrsSlitless_ComputePositions():
    config = default_init.copy()
    config.update({'sci_aper': 'MIRIM_SLITLESS'})
    cp = miri_lrs_slitless_tools.ComputePositions(initial_values=config)
    fig = cp.plot_scene()
    whippot_tools.plt.close(fig)
    return

def test_MiriWFSS_ComputePositions():
    config = default_init.copy()
    config.update({'sci_aper': 'MIRIM_ILLUM'})
    cp = miri_wfss_tools.ComputePositions(initial_values=config)
    fig = cp.plot_scene()
    whippot_tools.plt.close(fig)
    return

def test_MiriMRS_ComputePositions():
    config = default_init.copy()
    config.update({
        'sci_aper': 'MIRIFU_CHANNEL2C',
        'show_diffraction_spikes': True,
    })
    cp = miri_mrs_tools.ComputePositions(initial_values=config)
    fig = cp.plot_scene()
    whippot_tools.plt.close(fig)
    return

@pytest.mark.parametrize(
    'sci_aper',
    ['MIRIM_MASK1550','MIRIM_MASKLYOT'],
)
def test_MiriCoron_ComputePositions(sci_aper):
    config = default_init.copy()
    config.update({
        # update the aperture
        'sci_aper': sci_aper,
        # add a TA star
        'acq_ra': config['sci_ra'] + 60*units.arcsec.to(units.deg),
        'acq_dec': config['sci_dec'] - 30*units.arcsec.to(units.deg),
        # show diffraction spikes
        'show_diffraction_spikes': True,
    })
    cp = miri_coron_tools.ComputePositions(initial_values=config)
    fig = cp.plot_scene()
    whippot_tools.plt.close(fig)
    return


def test_NirspecIFU_ComputePositions():
    config = default_init.copy()
    config.update({
        'instr': 'nirspec',
        'sci_aper': 'nrs_full_ifu',
        'show_diffraction_spikes': True,
    })
    cp = nirspec_ifu_tools.ComputePositions(initial_values=config)
    fig = cp.plot_scene()
    whippot_tools.plt.close(fig)
    return



@pytest.mark.parametrize(
    'apername',
    sorted(lom.list_of_masks.keys())
)
def test_masks(apername):
    config = default_init.copy()
    config.update({
        'instr': instr_abbr[apername[:3]],
        'sci_aper': apername
    })
    print(config)
    cp = whippot_tools.ComputePositions(initial_values=config)
    fig = cp.plot_scene()
    whippot_tools.plt.close(fig)
    return

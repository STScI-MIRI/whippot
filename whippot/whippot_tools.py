import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units

import pysiaf
from pysiaf import Siaf

import ipywidgets as widgets

from whippot import whippot_plots
from whippot import list_of_masks as lom

#------------------------------------------------------#
#----------------- Offset Computation -----------------#
#------------------------------------------------------#

# define user interface
instruments = ['NIRCAM', 'NIRSPEC', 'NIRISS', 'MIRI', 'FGS']

class ComputePositions():
    def __init__(self, initial_values={}):
        """
        Generate the Position and Offset Computation GUI.
        Can initialize from a dictionary containing the entries (e.g.):
        initial_values={
          'instr': 'MIRI',
          'sci_aper': 'MIRIM_CORON1550',
          'pa': 180.,
          'acq_ra' : 90., 'acq_dec' : 90.,
          'sci_ra' : 91., 'sci_dec' : 89.,
          'other_stars' : '',
          'exclude_roi' : True,
          'show_diffraction_spikes': False,
        }
        """
        self._initial_values = initial_values.copy()
        # if no acq star coordinates are given, copy the sci star
        self._initial_values['acq_ra'] = initial_values.get("acq_ra", initial_values.get("sci_ra", 0))
        self._initial_values['acq_dec'] = initial_values.get("acq_dec", initial_values.get("sci_dec", 0))
        self.parameter_values = self.default_parameters()
        self.parameter_values.update(self._initial_values)
        self.instr = None
        self.aperture = None
        # are we doing TA on the same star or a different star?
        self.SELF_TA = True
        self.ui = self._make_ui()
        # if an initial dictionary is provided, run the computations.
        if initial_values != {}:
            self.compute_positions()

    def default_parameters(self) -> dict:
        defaults = {
          'instr': 'MIRI',
          'sci_aper': 'MIRIM_CORON1550',
          'pa': 0.,
          'acq_ra' : 0., 'acq_dec' : 0.,
          'sci_ra' : 0., 'sci_dec' : 0.,
          'other_stars' : '',
          'exclude_roi' : True,
          'show_diffraction_spikes' : False,
        }

        return defaults

    def filter_aperture_options(self):
        apertures = Siaf(self._instr_picker.value).apertures
        apernames = [name.upper() for name, aper in apertures.items()]
        if self._exclude_roi_chkbx.value:
            apernames = [name.upper() for name, aper in apertures.items() if aper.AperType != 'ROI']
        return apernames

    def _update_aperture_options(self, *args):
        self._sci_apers = self.filter_aperture_options()#[i for i in Siaf(self._instr_picker.value).apernames]
        self._sci_aper_picker.options = self._sci_apers
        curr_aper = self.parameter_values['sci_aper'].upper()
        if curr_aper in self._sci_apers:
            self._sci_aper_picker.value = curr_aper
        else:
            self._sci_aper_picker.value = self._sci_apers[0]

    def _update_parameter_dict(self):
        # read in the GUI fields to the parameter dictionary
        self.parameter_values['instr'] = self._instr_picker.value
        self.parameter_values['sci_aper'] = self._sci_aper_picker.value
        self.parameter_values['pa'] = self._PA_setter.value
        self.parameter_values['final_idl_x'] = self._slew_to_this_idl.children[1].children[0].value
        self.parameter_values['final_idl_y'] = self._slew_to_this_idl.children[1].children[1].value
        self.parameter_values['acq_ra'] = self._acq_pos_widget.children[1].children[0].value
        self.parameter_values['acq_dec'] = self._acq_pos_widget.children[1].children[1].value
        self.parameter_values['sci_ra'] = self._sci_pos_widget.children[1].children[0].value
        self.parameter_values['sci_dec'] = self._sci_pos_widget.children[1].children[1].value
        self.parameter_values['other_stars'] = self._other_stars_widget.value
        self.parameter_values['exclude_roi'] = self._exclude_roi_chkbx.value
        self.parameter_values['show_diffraction_spikes'] = self._show_spikes_chkbx.value

        # update various object attributes from the parameter dictionary
        self.instr = Siaf(self.parameter_values['instr'])
        self.aperture = self.instr[self.parameter_values['sci_aper']]
        # set the self-ta flag
        acq_ra, acq_dec = self.parameter_values['acq_ra'], self.parameter_values['acq_dec'] 
        sci_ra, sci_dec = self.parameter_values['sci_ra'], self.parameter_values['sci_dec'] 
        if (acq_ra == sci_ra) and (acq_dec == sci_dec):
            self.SELF_TA = True
        else:
            self.SELF_TA = False

    def _update_widgets_from_parameters(self):
        self._instr_picker.value = self.parameter_values['instr']
        self._sci_aper_picker.value = self.parameter_values['sci_aper']
        self._PA_setter.value = self.parameter_values['pa']
        self._slew_to_this_idl.children[1].children[0].value = self.parameter_values['final_idl_x']
        self._slew_to_this_idl.children[1].children[1].value = self.parameter_values['final_idl_y']
        self._acq_pos_widget.children[1].children[0].value = self.parameter_values['acq_ra']
        self._acq_pos_widget.children[1].children[1].value = self.parameter_values['acq_dec']
        self._sci_pos_widget.children[1].children[0].value = self.parameter_values['sci_ra']
        self._sci_pos_widget.children[1].children[1].value = self.parameter_values['sci_dec']
        self._other_stars_widget.value = self.parameter_values['other_stars']
        self._exclude_roi_chkbx.value = self.parameter_values['exclude_roi']
        self._show_spikes_chkbx.value = self.parameter_values['show_diffraction_spikes']

        # update various object attributes from the parameter dictionary
        self.instr = Siaf(self.parameter_values['instr'])
        self.aperture = self.instr[self.parameter_values['sci_aper']]
        # set the self-ta flag
        acq_ra, acq_dec = self.parameter_values['acq_ra'], self.parameter_values['acq_dec'] 
        sci_ra, sci_dec = self.parameter_values['sci_ra'], self.parameter_values['sci_dec'] 
        if (acq_ra == sci_ra) and (acq_dec == sci_dec):
            self.SELF_TA = True
        else:
            self.SELF_TA = False

    def get_aper(self):
        if hasattr(self, "aperture"):
            aper = self.aperture
        else:
            aper = Siaf(self._instr_picker.value)[self._sci_aper_picker.value]
        return aper

    def _make_starpos_widget(self, title, initial_ra=0., initial_dec=0.):
        """Make a widget for getting a star's RA and Dec"""
        star_widget = widgets.VBox([
            widgets.Label(value=title, layout = widgets.Layout(display='flex', justify_content='center')),
            widgets.VBox([
                widgets.FloatText(value=initial_ra, description='RA [deg]', disabled=False),
                widgets.FloatText(value=initial_dec, description='Dec [deg]', disabled=False)
            ])
        ])
        return star_widget

    def _make_final_idl_widget(self, title, initial_x=0., initial_y=0.):
        """Make a widget for getting a star's RA and Dec"""
        star_widget = widgets.VBox([
            widgets.Label(value=title, layout = widgets.Layout(display='flex', justify_content='center')),
            widgets.VBox([
                widgets.FloatText(value=initial_x, description='Final IDL X', disabled=False),
                widgets.FloatText(value=initial_y, description='Final IDL Y', disabled=False)
            ])
        ])
        return star_widget


    def _parse_other_stars(self) -> list:
        """Parse the entries in the other_stars_widget and return a list"""
        other_stars = []
        if self.parameter_values['other_stars'] != '':
            stars = self.parameter_values['other_stars'].strip().split("\n")
            stars = [dict(zip(['label', 'position'], i.split(":"))) for i in stars]
            for star in stars:
                position = SkyCoord(
                    *[float(i) for i in star['position'].strip()[1:-1].split(",")],
                    frame='icrs', unit='deg',
                )
                star['position'] = position
            other_stars = stars
        return other_stars

    def _make_widgets(self, widget_values={}):
        """
        Container method for making and initializing the widgets
        """
        self._instr_picker = widgets.Dropdown(
            options = instruments,
            value = widget_values.get('instr', instruments[0]).upper(),
            description='Instrument'
        )
        # run update_apers when instr_picker changes
        self._instr_picker.observe(self._update_aperture_options)
        # set self.sci_apers
        self._sci_aper_picker = widgets.Dropdown(description='Aperture')
        self._exclude_roi_chkbx = widgets.Checkbox(
            value = widget_values.get('exclude_roi', True),
            description='Exclude ROI apers',
            disabled=False,
            indent=False
        )
        self._exclude_roi_chkbx.observe(self._update_aperture_options)
        # initialize the apers
        self._update_aperture_options()
        # to show or not to show the full aperture list


        self._show_spikes_chkbx = widgets.Checkbox(
            value = widget_values.get('show_diffraction_spikes', False),
            description='Show diffraction spikes (WIP)',
            disabled=False,
            indent=False
        )

        # Position Angle
        self._PA_setter = widgets.BoundedFloatText(
            value=widget_values.get("pa", 0),
            min=0.,
            max=360.,
            step=0.1,
            description='PA (deg):',
            disabled=False
        )
        self._slew_to_this_idl = self._make_final_idl_widget(
            "Slew to this IDL",
            initial_x = widget_values.get("final_idl_x", 0.),
            initial_y = widget_values.get("final_idl_y", 0.),
        )
        # star positions
        self._acq_pos_widget = self._make_starpos_widget(
            "ACQ target position",
            initial_ra = widget_values.get("acq_ra", 0.),
            initial_dec = widget_values.get("acq_dec", 0.),
        )
        self._sci_pos_widget = self._make_starpos_widget(
            "SCI target position",
            initial_ra = widget_values.get("sci_ra", 0.),
            initial_dec = widget_values.get("sci_dec", 0.)
        )
        self._other_stars_widget = widgets.Textarea(
            value=widget_values.get("other_stars", '').strip(),
            placeholder='label : (ra.deg, dec.deg)',
            description='Other stars: ',
            disabled=False,
        )

        self._compute_positions_button = widgets.Button(
            description='Compute positions',
            disabled=False,
            button_style='success'
        )
        self._compute_positions_button.on_click(self.compute_positions)
        self._plot_scene_button = widgets.Button(
            description = 'Plot scenes',
            disabled=False,
            button_style = 'info'
        )
        self._plot_scene_button.on_click(self.plot_scene)

        self._output_offset = widgets.Output()
        self._output_before = widgets.Output()
        self._output_after = widgets.Output()

    def _make_ui(self):
        self._make_widgets(self.parameter_values)
        grid = widgets.GridspecLayout(
            n_rows=11, n_columns=3,
            style=dict(background='white')
        )
        grid[0, :] = widgets.Label(
            value="IDL Coordinate and Offset TA Calculator".upper(),
            layout = widgets.Layout(display='flex', justify_content='center'),
        )
        # grid[1, 0] = self._instr_picker
        grid[1, 0] = self._instr_picker
        grid[2, 0] = self._sci_aper_picker
        grid[3, 0] = widgets.Label(value='Options', layout = widgets.Layout(display='flex', justify_content='center'),)
        # toggles
        grid[4:10, 0] = widgets.VBox(
            [self._exclude_roi_chkbx,
             self._show_spikes_chkbx],
            layout=widgets.Layout(justify_content='flex-start'),
        )
        # position column
        grid[1:4, 1] = self._acq_pos_widget
        grid[4:7, 1] = self._sci_pos_widget
        grid[7:10, 1] = self._slew_to_this_idl

        # other stars text entry
        grid[1:-1, 2] = self._other_stars_widget

        # place the buttons at the bottom of the central column
        grid[-1, 1] = widgets.HBox([self._compute_positions_button, self._plot_scene_button])
        output_grid = widgets.GridspecLayout(
            n_rows=1, n_columns=3,
        )
        output_grid[0, 0] = self._output_offset
        output_grid[0, 1] = self._output_before
        output_grid[0, 2] = self._output_after
        ui = widgets.VBox([
            grid, output_grid
            # widgets.HBox([, self.output_before, self.output_after]),
        ])
        return ui

    def plot_scene(self, *args, frame='idl') -> mpl.figure.Figure:
        """
        Plot the detector- and sky-oriented scenes with the footprint of the selected aperture

        Parameters
        ----------
        *args is a dummy argument used to make the function compatible with ipywidgets calls
        """
        nrows = 1 if self.SELF_TA else 2
        ncols = 2
        fig, axes = plt.subplots(
            nrows=nrows,
			ncols=ncols,
			figsize=(5*ncols, 5*nrows),
			layout='constrained',
			squeeze=False
        )

        aper = self.get_aper()

        fig.suptitle(aper.AperName)

        # the function to generate mask patches
        mask_func = lom.list_of_masks.get(aper.AperName.upper(), None)

        i = 0
        if not self.SELF_TA:
            fig = whippot_plots.plot_aper_to_frame(
                self.get_aper(),
                self.idl_coords_after_ta,
                frame_from='idl',
                frame_to=frame,
                ax = axes[i, 0],
                title='Detector-oriented view\n(ACQ star centered)',
                show_legend = False,
                idl_mask=lom.make_mask(mask_func),
            )
            fig = whippot_plots.plot_aper_to_frame(
                self.get_aper(),
                self.idl_coords_after_ta,
                frame_from='idl',
                frame_to='sky',
                ax = axes[i, 1],
                title='Sky-oriented view\n(ACQ star centered)',
                show_legend = False,
                idl_mask=lom.make_mask(mask_func),
            )
            i += 1
        fig = whippot_plots.plot_aper_to_frame(
            self.get_aper(),
            self.idl_coords_after_slew,
            frame_from='idl',
            frame_to=frame,
            ax = axes[i, 0],
            title='Detector-oriented view',
            show_legend = True,
            idl_mask=lom.make_mask(mask_func),
        )
        fig = whippot_plots.plot_aper_to_frame(
            self.get_aper(),
            self.idl_coords_after_slew,
            frame_from='idl',
            frame_to='sky',
            ax = axes[i, 1],
            title='Sky-oriented view',
            show_legend = False,
            idl_mask=lom.make_mask(mask_func),
        )
        # fig = whippot_plots.plot_aper_sky(
        #     self.get_aper(),
        #     self.idl_coords_after_slew,
        #     ax = axes[i, 1],
        #     title='Sky-oriented view',
        #     show_legend = False,
        #     idl_mask=lom.make_mask(mask_func),
        # )

        return fig

    def show_ui(self):
        return self.ui

    def compute_positions(self, *args, update_params_from_widgets=True):
        # first, update either the values config dictionary or the widgets
        if update_params_from_widgets:
            self._update_parameter_dict()
        acq_ra = self.parameter_values['acq_ra']
        acq_dec = self.parameter_values['acq_dec']
        sci_ra = self.parameter_values['sci_ra']
        sci_dec = self.parameter_values['sci_dec']

        acq_pos = {
            'label': 'ACQ',
            'position': SkyCoord(
                acq_ra, acq_dec,
                frame='icrs', unit='deg',
            ),
        }
        sci_pos = {
            'label': 'SCI',
            'position': SkyCoord(
                sci_ra, sci_dec,
                frame='icrs', unit='deg',
            )
        }

        v3pa = self.parameter_values['pa']

        other_stars = self._parse_other_stars()
        # any extra IDL slews after TA is complete
        slew_to_idl = np.array([
            self.parameter_values['final_idl_x'],
            self.parameter_values['final_idl_y'],
        ])

        idl_coords = compute_idl_after_ta(
            acq_pos, sci_pos, v3pa, self.aperture,
            other_stars = other_stars,
        )
        self.idl_coords_after_ta = {i['label']: i['position'] for i in idl_coords}
        self.offset_to_sci = -self.idl_coords_after_ta['SCI'] + slew_to_idl
        self.idl_coords_after_slew = {
            k: v + self.offset_to_sci
            for k, v in self.idl_coords_after_ta.items()
        }
        # set the aperture attitude matrix to the SCI position
        create_attmat(sci_pos['position'], self.aperture, v3pa, set_matrix=True)
        # if ACQ and SCI stars are the same, remove the SCI star
        if self.SELF_TA == True:
            self.idl_coords_after_ta.pop("ACQ")
            self.idl_coords_after_slew.pop("ACQ")

        self._output_offset.clear_output()
        self._output_before.clear_output()
        self._output_after.clear_output()
        outputstr = "Special Requirement -> Offset values\n"
        outputstr += "-"*(len(outputstr)-1) + "\n"
        outputstr += f"Offset X [arcsec]: {self.offset_to_sci[0]:+0.4f}\nOffset Y [arcsec]: {self.offset_to_sci[1]:+0.4f}"
        self._output_offset.append_stdout(outputstr)
        outputstr = f"IDL positions of stars after TA:\n"
        outputstr += "-"*(len(outputstr)-1) + "\n"
        outputstr += "\n".join(f"{k}\t{v[0]:>+10.4f}\t{v[1]:>+10.4f}" for k, v in self.idl_coords_after_ta.items())
        self._output_before.append_stdout(outputstr)
        outputstr = f"IDL positions of stars after slew:\n"
        outputstr += "-"*(len(outputstr)-1) + "\n"
        outputstr += "\n".join(f"{k}\t{v[0]:>+10.4f}\t{v[1]:>+10.4f}" for k, v in self.idl_coords_after_slew.items())
        self._output_after.append_stdout(outputstr)


#------------------------------------------------------#
#----------------- Offset Computation -----------------#
#------------------------------------------------------#
def compute_idl_after_ta(
    slew_from: dict,
    slew_to: dict,
    v3pa: float,
    sci_aper : pysiaf.aperture.JwstAperture,
    other_stars : list = [],
) -> list[dict]:
    """
    Compute the IDL positions of all the given stars whent he ACQ target is at the reference position.

    Parameters
    ----------
    slew_from: dict
      A dictionary containing the label and position of the TA target, set by
      the user in compute_offsets.py
    slew_to: dict
      A dictionary containing the label and position of the science target, set by
      the user in compute_offsets.py
    v3pa: float
      The PA_V3 angle of the telescope for this observation
    coron_id : str
      Identifier for the desired coronagraphic subarray. Must be one of '1065',
      '1140', '1550', and 'LYOT'
    other_stars : list
      A list of dicts of other stars in the field, in the same format as slew_from/slew_to
    verbose : int = 1
      print diagnostics and offsets to screen.
      0 : nothing is printed
      1 : all output is printed
      2 : only the final IDl coordinates of the stars are printed

    Output
    ------
    Returns a list of dicts of floats

    """
    # make sure coron_id is valid
    star_positions = [slew_from, slew_to] + other_stars
    idl_coords = sky_to_idl(star_positions,
                            sci_aper,
                            v3pa)
    return idl_coords

def apply_offset(
    idl_coords : dict,
) -> dict:
    sci_idl = idl_coords.pop('SCI')
    offset_to_sci = -1 * sci_idl
    final_idl = {l: c + offset_to_sci for l, c in idl_coords.items()}
    return final_idl

def compute_offsets(
        slew_from: dict,
        slew_to: dict,
        v3pa: float,
        sci_aper : pysiaf.aperture.JwstAperture,
        other_stars : list = [],
        verbose : int = 1,
        return_offsets : bool = False,
) -> np.ndarray :
    """
    Compute the slews for the TA sequences, print the offsets, and show the plots if requested.
    How it works:
    - Point the coronagraphic aperture (MASK or CORON) at the TA star by
      setting an attitude matrix for the V3PA value
    - The offset is the negative of the IDL coordinates of the SCI star
    - The rest of the machinery is basically just making verification plots

    Parameters
    ----------
    slew_from: dict
      A dictionary containing the label and position of the TA target, set by
      the user in compute_offsets.py
    slew_to: dict
      A dictionary containing the label and position of the science target, set by
      the user in compute_offsets.py
    v3pa: float
      The PA_V3 angle of the telescope for this observation
    coron_id : str
      Identifier for the desired coronagraphic subarray. Must be one of '1065',
      '1140', '1550', and 'LYOT'
    other_stars : list
      A list of dicts of other stars in the field, in the same format as slew_from/slew_to
    verbose : int = 1
      print diagnostics and offsets to screen.
      0 : nothing is printed
      1 : all output is printed
      2 : only the final IDl coordinates of the stars are printed
    show_plots : bool = True
      If True, display the diagnostic plots. If False, only print the offsets.
    return_offsets : bool = False
      If True, return an array of dx and dy offsets

    Output
    ------
    Prints offsets and shows plots. Returns a dict of floats

    """
    # make sure coron_id is valid
    star_positions = [slew_from, slew_to] + other_stars
    labels = {'ACQ': slew_from['label'],
              'SCI': slew_to['label']}


    # Offsets to science target
    sep = star_positions[0]['position'].separation(
        star_positions[1]['position']
    ).to(units.arcsec)
    pa = star_positions[0]['position'].position_angle(
        star_positions[1]['position']
    ).to(units.deg)

    idl_coords = sky_to_idl(star_positions,
                            sci_aper,
                            v3pa)
    # The offset you apply is as if you were moving the science target - i.e.
    # the negative of its position
    offset = -1*np.array(idl_coords[1]['position'])

    if verbose == 1:
        len_label = max(len(star['label']) for star in idl_coords)
        print("IDL coordinates of all stars after slew:")
        for star in idl_coords:
            label = star['label']
            idl = star['position']
            print(f"{label:{len_label}s}:\t{idl[0]+offset[0]:+0.10f}, {idl[1]+offset[1]:+0.10f}")
        print("")

    if return_offsets == True:
        return offset

def create_attmat(
        position : SkyCoord,
        aper : pysiaf.aperture.JwstAperture,
        pa : float,
        idl_offset : tuple[float, float] = (0., 0.),
        set_matrix : bool = True,
) -> np.ndarray:
    """
    Create an attitude matrix for JWST when the reference point of a particular
    aperture is pointed at a given position for a specified PA

    Parameters
    ----------
    position : SkyCoord
      skycoord position on the sky
    aper : pysiaf.aperture.JwstAperture
      pySIAF-defined aperture object
    pa : float
      PA angle with respect to the V3 axis measured at the aperture reference
      point. This corresponds to the V3PA field in the APT PA range special
      requirement, and the ROLL_REF keyword in the data. This is *not* the
      PA_V3 keyword value, which is the PA angle of the V3 axis measured at the
      telescope boresight.
    idl_offset : tuple[float, float] = (0.0, 0.0)
      allows you to specify an arbitrary position in IDL coordinates that
      corresponds to the position
    set_matrix : bool = True
      if True, also set the matrix on the current aperture in addition to returning it
    Output
    ------
    attmat : np.ndarray
      matrix that pySIAF can use to specify the attitude of the telescope
    """
    v2, v3 = aper.idl_to_tel(idl_offset[0], idl_offset[1])
    # v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3,
                                                    ra=position.ra.deg,
                                                    dec=position.dec.deg,
                                                    pa=pa)
    if set_matrix == True:
        aper.set_attitude_matrix(attmat)
    return attmat


def sky_to_idl(
        stars : list[dict],
        aper : pysiaf.aperture.JwstAperture,
        pa : float,
        idl_offset : tuple[float, float] = (0., 0.)
) -> list[dict]:
    """
    Convert RA and Dec positions of a TA star and its target with an offset
    into a detector position (measured from the reference point, in arcsec)
    for a given PA of the V3 axis w.r.t. North.
    Assume the TA star is centered on the aperture (idl = (0, 0))

    Parameters
    ----------
    stars : list of {"label": label, "position": pos} elements to plot
      the first element is the ACQ target. The aperture is centerd on this target.
    aper : SIAF object for the aperture being used (e.g. MIRIM_MASK1550)
    pa : the PA in degrees of the V3 axis of the telescope (measured eastward of North) at the observation epoch

    Output
    ------
    idl_coords : dict {label: idl_tuple} of IDL coordinates for each target
      a dictionary of IDL x, y coordinates for the TA star and the science target
      the TA star coordinates should be very close to 0
    """
    acq_pos = stars[0]['position']
    attmat = create_attmat(acq_pos, aper, pa, idl_offset)
    aper.set_attitude_matrix(attmat)
    idl_coords = []
    # ta star - should be close to 0
    for star in stars:
        label = star['label']
        pos = star['position']
        idl_coords.append({'label': label, 'position': np.array(aper.sky_to_idl(pos.ra.deg, pos.dec.deg))})
    return idl_coords



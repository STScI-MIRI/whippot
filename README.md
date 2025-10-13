# Offset TA Tools

Author: Jonathan Aguilar (jaguilar@stsci.edu)

Last update: Oct 13, 2025

## Requirements ##
- ipywidgets
- numpy
- matplotlib
- astropy
- pySIAF

## Installation ##

There are two ways to use these tools, after downloading the repository: 
1. Navigate to the top-level directory, and run `pip install .`.
   - This enables import via `from whippot import whippot_tools`.
2. Nagivate to the `whippot` subfolder, and copy the notebook(s) and the python scripts to whatever directory you plan to be working in. 
   - In this case, the statement to import `whippot_tools` must be changed to `import whippot_tools`.
   - These files can also be copied directly from the online repository, without downloading the rest.

## All-Instrument Interface Usage ##

For the purposes of this README, the word "aperture" is used in the same sense as it is in the SIAF (Science Instrument Aperture File) - any defined region of the telescope that can be used to command the telescope pointing. Some of these "apertures" also correspond to subarrays that are read out as data.

The "All-Instrument Interface" is much easier to use than the previous script-based version, and is recommended. Some highlights, while the README is slowly updated:
1. This is GUI-based, instead of script-based.
2. Any aperture can be chosen from any instrument, no longer just MIRI. Unfortunately, no help is given to the user to figure out which apertures are the correct ones to use. Some options are:
   - Guess based on the names;
   - Check the APERNAME keyword of data from the same observing configuration you plan to use;
   - Poke around on JDox for any mentions of the SIAF;
   - Ask an experienced JWST user (_very_ experienced).
4. An initialization dictionary can be provided. If none is, the GUI will initialize with default values. Initialization dictionaries are useful for repeatability.
5. The SCI and ACQ target positions are intended for use with offset target acquisition (that is, TA is performed on a different target from the science target).
  - If performing self-TA, use the same coordinates for the ACQ and SCI position fields.
  - the `acq_ra` and `acq_dec` fields can be omitted from the initialization dictionary, in which case the SCI target positions will be copied over.
5. If you want your SCI star to land somewhere in the aperture other than the reference position, enter this in arcsec into the "Final IDL X" and "Final IDL Y" positions. Read these as, "This is the final position (in IDL X and Y) where I want my SCI target to end up"., enter this in arcsec into the "Final IDL X" and "Final IDL Y" positions. Read these fields as, "This is the final position (in IDL X and Y) where I want my SCI target to end up".


!! Attention !!
It is not recommended to use a ComputePositions instance for multiple calculations. Instead, it is recommended to create a new ComputePositions object for each set of parameters you wish to compare or manipulate.

### Explanation of fields ###

- Instrument : one of the JWST instruments ('NIRCAM', 'NIRSPEC', 'NIRISS', 'MIRI', or 'FGS')
- SCI aperture : the SIAF-defined name of the part of the telescope whose reference position will be pointed at the SCI target. Choose from the drop-down menu.
- PA (deg) : the PA angle of the V3 axis of the telescope. Not to be confused with ROLL_REF.
- ACQ target position RA/Dec : the RA and Dec positions of the target used for Target Acquisition, in decimal degrees
- SCI target position RA/Dec : the RA and Dec positions of the final science target, in decimal degrees
- Other stars : a mutliline string containing a label and decimal-degree coordinates for any other targets whose positions you wish to know. The format is: "label: (ra, dec)". 
- Final IDL X/Y: After TA, this is where you want your SCI target to end up (in IDL X/Y arcsec, corresponding to APT's "Special Requirements ->  Offsets" option).

### Initialization dictionary keywords ###

You may prefer to pass an initialization dictionary to the ComputePositions object, instead of setting each field by hand. These are the available keywords, and the field above to which they correspond:
- instr -> Instrument 
- sci\_aper -> SCI aperture 
- pa -> PA
- sci_ra -> SCI target position / RA [deg]
- sci_dec -> SCI target position / Dec [deg]
- acq_ra -> ACQ target position / RA [deg]
- acq_dec -> ACQ target position / Dec [deg]
- other_stars -> Other Stars
- final\_idl\_x -> Final IDL X
- final\_idl\_y -> Final IDL Y

## (Old format) USAGE ##

This library is meant to help users compute offsets in the detector frame of
reference to enter into APT as Offset Special Requirements, in the case that
they need to perform TA on a target other than their science target. Users
should edit the `if __name__ == '__main__:` section at the bottom of
`compute_offsets.py` with the coordinates of their TA and Science targets, as
well as the position angle of the V3 axis of the telescope. The instructions for
editing each section are written in comments in the script. This will print out
appropriate X and Y offsets that can be entered into APT, along with some
diagnostic information. Alternatively, the file `example.py` shows how to use
`compute_offsets.py` as an imported module.

It is suggested that users make a new copy of `compute_offsets.py` or
`example.py` for each variation of TA star and science target.

It can also be imported as a module after downloading the source code from github:
```
cd $directory
pip install .
```

Then from a python terminal:
```
> from miri_coro_offset_ta import compute_offsets
```

### User-specified variables ###

#### `slew_to`, `slew_from` ####

`slew_to` and `slew_from` are dictionaries used to specify the science (SCI) and
acquisition (ACQ) targets, respectively.
- The `label` keyword is used for annotating output text and figures.
- The `position` keyword is used to store an astropy `SkyCoord` object
  ([SkyCoord
  documentation](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html)).
  This gives considerable flexibility in specifying the coordinates; for
  example, by providing a distance and proper motions, the user can propagate
  the positions using `SkyCoord.apply_space_motion()` to compute offsets for
  multiple epochs.


#### `v3pa` ####

`v3pa` refers to the position angle of the telescope's V3 axis *at the reference
position of the aperture used for the observation*. If that sounds confusing,
the short version is it corresponds to the angle in APT's `Special Requirements
-> PA -> PA Range` menu if the `V3PA` radio button is selected. It also
corresponds to the `ROLL_REF` header keyword (see the [JWST Keyword
Dictionary](https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/) for
keyword definitions). If you are using this library to plan observations in APT,
the `V3PA` field should match your `v3pa` variable.

This is not to be confused with the `PA_APER` header keyword, which corresponds
to the `Aperture PA Range` radio button and refers to the amount by which the
detector-aligned coordinate system is rotated with respect to the `V3` axis. It
also is not to be confused with the `PA_V3` header keyword, which refers to the
V3 position angle at the position of the telescope boresight. Due to spherical
trigonometric effects, the PA of the V3 axis varies across the telescope's focal
plane and varies strongly at high and low latitudes.

Offset slews are specified along the detector axes, in units of arcsec (see
https://jwst-docs.stsci.edu/jppom/special-requirements/general-special-requirements).
In order to convert between the detector coordinate system and two positions on
the sky, pySIAF requires information about the orientation of the telescope.
Here, we provide this information using a combination of the coronagraph used
(see `coron_id`), and position angle of the v3 axis of the telescope, measured
at the chosen coronagraph's reference position.

More details about the different coordinate systems used in describing positions
in the telescope can be found here:
https://jwst-docs.stsci.edu/jwst-observatory-characteristics/jwst-observatory-coordinate-system-and-field-of-regard/
.

## FAQs for creating your APT program ##

### How do I choose my acquisition target? ###

To choose an acquisition target, you should consider the brightness, separation,
and position angle:
- Brightness: it should be bright enough to achieve high SNR in the TA filter
  without saturating (see the [ETC](https://jwst.etc.stsci.edu/)).
- Separation : it should be closer than the [visit-splitting
  distance](https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview/additional-jwst-apt-functionality/apt-visit-splitting),
  which ranges between 30"-80" depending on the availability of guide stars for
  a particular target.
- Position angle: The acquisition target must be clear of diffraction spikes
  from nearby sources. This is especially important if the science target is
  very bright. [STPSF](https://stpsf.readthedocs.io/) can be used to
  determine if your TA target will be clear of the diffraction spikes. You can place the field stars into your simulated imaging by providing the IDL coordinates computed by WHIPPOT for a given PA.

### How do I choose dates and V3PA angles? ###

To see available dates and V3PA angles, go to the Visit Planner window in APT
and find the `Reports` menu at the bottom. Select a visit, and then select
`Total Roll Analysis for Visit`. This will give you a plot of available V3 PA
angles against dates, as well as a table that can be read into a script.

### Do I need to calculate a separate offset for each roll? ###

Yes, unless your roll angle is very small or your acquisition target is very
close.

### How to I place my SCI target at the position of one of the field targets in a separate observation?

This situation may arise, for example, in a high contrast imaging scenario where you have a nearby bright star (let's call it, BS) in field of view of your science observation, and you wish to create a PSF reference observation to subtract it out. Let's call the science observation with the interloper Obs 1, and the PSF reference observation, Obs 2.

For Obs 1, you would create a `ComputePositions` object and save it to a variable, like so: `ob1 = ComputePositions()`. Fill the ACQ and SCI RA/Dec fields as normal, and add "BS: (bs.ra, bs.dec)" to the "Other stars" entry. Press "Compute positions" to get the IDL coordinates of BS in your science observation.

Then, create a new `ComputePositions` object for Obs 2: `obs2 = ComputePositions()`. For the SCI field, enter the RA and Dec of your PSF reference star, and then enter the IDL X and Y of BS into the "Final IDL X" and "Final IDL Y" fields. These are accessible by copy-pasting from the output of Obs 1, or by reading from the `obs1.idl_coords_after_slew['bs']` tuple. Press `Compute Positions` to get the values to input into APT's Special Requirements field. If the ACQ and SCI targets are the same, these values will be the same as the Final IDL X/Y values, but if the ACQ and SCI targets are different, they will not be.


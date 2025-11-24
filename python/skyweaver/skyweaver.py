"""
Delay generation and metadata parsing for MeerKAT baseband recordings
"""
from __future__ import annotations

# stdlib imports
import logging
import textwrap
import ctypes
from typing import Any, Tuple
from dataclasses import dataclass
from typing_extensions import Self
import sys
from collections import defaultdict
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm
from scipy.interpolate import RegularGridInterpolator
# 3rd party imports
import h5py
import yaml
import numpy as np
import pandas as pd
import random
from rich.progress import track
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy import wcs
import matplotlib.pyplot as plt
from katpoint import Target, Antenna
from mosaic.beamforming import (
    DelayPolynomial,
    PsfSim,
    generate_nbeams_tiling,
    Tiling, BeamShape)

log = logging.getLogger("skyweaver")


class DelayModelHeader(ctypes.Structure):
    """Structure for packing delay model information
    """
    # pylint: disable=too-few-public-methods
    _pack_ = 1  # byte alignment
    _fields_ = [
        # Number of beams
        ("nbeams", ctypes.c_uint32),
        # Number of antennas
        ("nantennas", ctypes.c_uint32),
        # UNIX epoch of start of delay model validity
        ("start_epoch", ctypes.c_double),
        # UNIX epoch of end of delay model validity
        ("end_epoch", ctypes.c_double),
    ]

@dataclass
class DelayModel:
    """Wrapper class for a delay model
    """
    # Start of validity epoch
    start_epoch: Time
    # End of validity epoch
    end_epoch: Time
    # weight/delay array in (Nbeam, Nant, 3) format
    delays: np.ndarray

    @property
    def nbeams(self) -> int:
        """Return the number of beams in the model

        Returns:
            int: The number of beams
        """
        return self.delays.shape[0]

    @property
    def nantennas(self) -> int:
        """Return the number of antennas in the model

        Returns:
            int: The number of antennas
        """
        return self.delays.shape[1]

    def to_bytes(self) -> bytes:
        """Pack the delay model into bytes

        Returns:
            bytes: Byte packed header and delays
        """
        header = DelayModelHeader(
            self.nbeams,
            self.nantennas,
            self.start_epoch.unix,
            self.end_epoch.unix
        )
        body = self.delays.astype("float32").tobytes()
        return bytes(header) + body

    def validate_epoch(self, epoch: Time) -> bool:
        """Check that the delay model is valid for the given epoch

        Args:
            epoch (Time): A time value to check

        Returns:
            bool: True = valid, False = invalid
        """
        return (epoch >= self.start_epoch) and (epoch <= self.end_epoch)

    def get_phase(self, frequencies: Quantity, epoch: Time) -> np.ndarray:
        """Return an array of beam weights (geometric only)
           for the given time and frequencies

        Args:
            frequencies (Quantity): An array of frequencies
            epoch (Time): An epoch at which to evaluate the delay model

        Raises:
            ValueError: Raised when the passed epoch is invalid for this model

        Returns:
            np.ndarray: A complex64 type numpy array with dimensions
            (nbeams, nantennas, nfrequencies)

        Note:
            The scalar weights set for the antennas will be applied during
            phase determination, such that:
            phase = weight * e^(-i * 2pi * delay * frequency)
        """
        if not self.validate_epoch(epoch):
            raise ValueError(
                "Epoch outside of delay polynomial validity window")
        epoch_diff = (epoch - self.start_epoch).to(u.s).value
        weights = self.delays[:, :, 0]
        offsets = self.delays[:, :, 1]
        rates = self.delays[:, :, 2]
        true_delay = offsets + rates * epoch_diff
        phases = weights[:, :, np.newaxis] * \
            np.exp(-1j * np.pi * 2 * true_delay[:, :, np.newaxis] * frequencies.to(u.Hz).value)
        return phases


def get_phases(delays: list[DelayModel], frequencies: Quantity,
               epochs: Time) -> list[tuple[Time, np.ndarray]]:
    """Get multiple phasing solutions from a list of delay models

    Args:
        delays (list[DelayModel]): A list of delay models
        frequencies (Quantity): The frequencies for which to calculate delays
        epochs (Time): The epochs at which to determine phases

    Raises:
        ValueError: Raised if any given epoch is invalid for all delay models

    Returns:
        list[tuple[Time, np.ndarray]]: An list of tuples containing the epoch and
                                       corresponding phases. Phases are  complex64
                                       type numpy array with dimensions (nbeams,
                                       nantennas, nfrequencies)
    """
    delay_epochs = [i.start_epoch for i in delays]
    phases = []
    for epoch in epochs:
        idx = np.searchsorted(delay_epochs, epoch, side="right")
        if idx == 0 or not delays[idx - 1].validate_epoch(epoch):
            print(idx, delays[idx - 1].validate_epoch(epoch))
            raise ValueError(f"Epoch {epoch} is out of range for delay models")
        phases.append((epoch, delays[idx - 1].get_phase(frequencies, epoch)))
    return phases


class DelayEngine:
    """A class for generating delay/weight solutions for beamforming
    """

    def __init__(self, subarray: Subarray, phase_reference: Target) -> None:
        """Create an instance of DelayEngine

        Args:
            subarray (Subarray): The superset of antennas that will be beamformed
            phase_reference (Target): The phase centre of the pointing
        """
        self.subarray = subarray
        self.phase_reference = phase_reference
        self._targets: list[tuple[Target, Subarray]] = []
        # subarray sets are used as a performance optimisation
        # for making maps of scalar antenna weights
        self._subarray_sets: list = []

    @property
    def nbeams(self) -> int:
        """Return the total number of beams specified

        Returns:
            int: The number of beams
        """
        return len(self._targets)

    @property
    def targets(self) -> list[Target]:
        """Get a the list of targets

        Returns:
            list[Target]: List of katpoint Target objects
        """
        return [i[0] for i in self._targets]

    def _validate_subarray(self, subarray: Subarray) -> None:
        if subarray not in self.subarray:
            raise ValueError("Subarray is not a sub-set of full array")

    def add_tiling(
            self,
            tiling: Tiling,
            subarray: Subarray = None,
            prefix: str = None) -> None:
        """Add a tiling to the engine

        Args:
            tiling (Tiling): A mosaic Tiling object containing beam positions
            subarray (Subarray, optional): The subarray that should be used when beamforming.
                                           Defaults to the full array.
            prefix (str, optional): A prefix to be used when naming beams from this tiling.
                                    Defaults to the name of the tiling target.
        """
        if subarray is None:
            subarray = self.subarray
        else:
            self._validate_subarray(subarray)
        if prefix is None:
            prefix = tiling.beam_shape.bore_sight.name
        coordinates = tiling.get_equatorial_coordinates()
        self._subarray_sets.append(subarray)
        sub_array_idx = len(self._subarray_sets) - 1
        for ii, (ra, dec) in enumerate(coordinates):
            self._targets.append(
                (Target(f"{prefix}_{ii:04d},radec,{ra},{dec}"), sub_array_idx))
    


    def add_beam(self, target: Target, subarray: Subarray = None) -> None:
        """Add a single beam to the engine

        Args:
            target (Target): A katpoint target specifying the beam position
            subarray (Subarray, optional): The subarray that should be used when beamforming.
                                           Defaults to the full array.
        """
        if subarray is None:
            subarray = self.subarray
        else:
            self._validate_subarray(subarray)
        self._subarray_sets.append(subarray)
        sub_array_idx = len(self._subarray_sets) - 1
        self._targets.append((target, sub_array_idx))

       
    def _extract_weights(self) -> np.ndarray:
        # Here we extract the weights for each beam/antenna
        # as an optimisation we cache the antenna mask per
        # unique subarray set.
        weights = np.zeros((len(self._targets),
                            self.subarray.nantennas),
                           dtype="float32")
        masks = {}
        for beam_idx, (_, subarray_idx) in enumerate(self._targets):
            if subarray_idx not in masks:
                subarray = self._subarray_sets[subarray_idx]
                mask = np.zeros(self.subarray.nantennas, dtype="float32")
                for ant_idx, antenna in enumerate(
                        self.subarray.antenna_positions):
                    if antenna in subarray.antenna_positions:
                        mask[ant_idx] = 1.0
                masks[subarray_idx] = mask
            weights[beam_idx, :] = masks[subarray_idx][:]
        return weights.reshape((len(self._targets), self.subarray.nantennas, 1))

    def calculate_delays(self, start: Time, end: Time,
                         step: TimeDelta) -> list[float, float, np.ndarray]:
        """Calculate weight/delay solutions for all specified beams
           in intervals over a time span.

        Args:
            start (Time): The start epoch
            end (Time): The end epoch
            step (TimeDelta): The time step between solutions

        Returns:
            list[float, float, np.ndarray]: A list of delay solutions with validity
                                            start times and end times.

        Notes:
            Delays are returned in the format:
            [(start epoch, end epoch, delay model)]

            The delay model itself is a numpy ndarray with dimensions
            (nbeams, nantennas, 3) where the inner 3-length dimension
            contains (scalar antenna weight, delay offset, delay rate).
        """
        targets = [i[0] for i in self._targets]
        weights = self._extract_weights()
        delay_calc = DelayPolynomial(
            self.subarray.antenna_positions, self.phase_reference,
            targets, self.subarray.reference_antenna)
        models = []
        epochs = []
        current_time = start
        while current_time <= end:
            epochs.append(current_time)
            current_time += step
        for epoch in track(epochs, description="Generating delays..."):
            delays = delay_calc.get_delay_polynomials(
                epoch.to_datetime(),
                duration=step.to(u.s).value)
            # Delays are returned with shape (nbeam, nant, 2)
            # where the innermost dimension is (delays, rate)
            # Here we inject the beam weight into the delay array
            # to make the inner dimension (antenna weight, delay, rate).
            # The weight here refers to a scalar multiplier for the antenna
            # and is intended as a means to allow different subsets of antennas
            # to be used for beamforming.
            delays = np.concatenate((weights, delays), axis=-1)
            models.append(DelayModel(epoch, epoch + step, delays))
        return models


class Subarray:
    """Class to wrap a set of antennas
    """

    def __init__(self, antenna_positions: list[Antenna]) -> None:
        """Create a Subarray instance

        Args:
            antenna_positions (list[Antenna]): A list of katpoint Antenna objects
            
        Note:
            All antenns must share a common reference antenna
        """
        self.antenna_positions = antenna_positions
        self._check_references()

    def _check_references(self):
        ref = self.antenna_positions[0].array_reference_antenna()
        for ant in self.antenna_positions[1:]:
            if ant.array_reference_antenna() != ref:
                raise ValueError(f"Antenna {ant} does not share reference "
                                 f"position with antenna {self.antenna_positions[0]}")

    @property
    def reference_antenna(self) -> Antenna:
        """Get the common reference antenna

        Returns:
            Antenna: The array reference antenna
        """
        return self.antenna_positions[0].array_reference_antenna()

    @property
    def nantennas(self) -> int:
        """Get the number of antennas in the subarray

        Returns:
            int: The antenna count
        """
        return len(self.antenna_positions)

    @property
    def names(self) -> list[str]:
        """Retrun the names of the antennas in the subarray

        Returns:
            list[str]: Antenna names
        """
        return [i.name for i in self.antenna_positions]

    def __contains__(self, subarray: Subarray) -> bool:
        other_set: set = set(subarray.antenna_positions)
        self_set: set = set(self.antenna_positions)
        return other_set.issubset(self_set)

@dataclass
class CalibrationSolution:

    epoch: str
    antenna_pols : list[str]
    solution : np.ndarray

    def __init__(self, epoch, antenna_pols, solution):

        self.epoch = epoch
        self.antenna_pols = antenna_pols
        self.solution = solution

    def to_file(self,basename="gains"):

        with open(basename + f"_{self.epoch}.afp",'wb') as out:
            self.solution.tofile(out)

@dataclass
class SessionMetadata:
    """Class wrapping the contents of the HDF5 metadata files from FBFUSE-BVR
    """
    # Mapping of antenna names to XEPhem metadata
    antenna_positions: dict[str, Antenna]
    # Mapping of antenna names to ordinal indices in CBF stream
    antenna_feng_map: dict[str, int]
    # Observation centre frequency
    centre_frequency: Quantity
    # Observation bandwidth
    bandwidth: Quantity
    # Total number of frequency channels
    nchannels: int
    # The project identifier
    project_id: str
    # The schedule block identifier
    sb_id: str
    # The version of the CBF used for the observation
    cbf_version: str
    # The UNIX epoch at which the system was synchronised
    sync_epoch: Time
    # Time series of phase centres used during observation
    phase_centres: np.ndarray[str, str]
    # Time series of "suspect" flags indicating data validity
    # True values imply invalid data
    suspect_flags: np.ndarray[bool]
    # Complex gain solutions (if no phase-up is done)
    calibration_solutions: list[CalibrationSolution]

    def __str__(self) -> str:
        key_padding = 20

        def keyval_format(key, value):
            lines = textwrap.wrap(str(value))
            if not lines:
                lines = ["-"]
            output = []
            first = True
            for line in lines:
                if first:
                    output.append(f"{key+':':<{key_padding}}{line}")
                    first = False
                else:
                    output.append(f"{'':<{key_padding}}{line}")
            return "\n".join(output)
        windows = self.find_observing_windows()
        windows_formatted = []
        for ii, (start, end, target, _) in enumerate(windows):
            start_sample = int(
                (start.unix -
                 self.sync_epoch.unix) *
                2 *
                self.bandwidth.to(
                    u.Hz).value)
            end_sample = int((end.unix - self.sync_epoch.unix)
                             * 2 * self.bandwidth.to(u.Hz).value)
            windows_formatted.append(
                f"#{ii:<4}{target.name:<20}{start} until {end} (UTC)")
            windows_formatted.append(
                f"{'':<25}{start.unix} until {end.unix} (UNIX)")
            windows_formatted.append(
                f"{'':<25}{start_sample} until {end_sample} (SAMPLE CLOCK)")
        output = (
            f"{'Array configuration':-^50}",
            keyval_format("Nantennas", len(self.antenna_positions)),
            keyval_format(
                "Subarray", ",".join(sorted(self.antenna_positions.keys()))),
            keyval_format("Centre frequency", self.centre_frequency),
            keyval_format("Bandwidth", self.bandwidth),
            keyval_format("Nchannels", self.nchannels),
            keyval_format("Sync epoch (UNIX)", self.sync_epoch),
            keyval_format("Project ID", self.project_id),
            keyval_format("Schedule block ID", self.sb_id),
            keyval_format("CBF version", self.cbf_version),
            f"{'Pointings':-^50}", "\n".join(windows_formatted))
        return "\n".join(output)

    @classmethod
    def from_file(cls, fname) -> Self:
        """Parse an HDf5 metdata file

        Args:
            fname (str): Path to the metadata file to be read

        Returns:
            Self: SessionMetadata instance
        """
        with h5py.File(fname) as f:
            # Parse out antenna to F-engine mapping
            antenna_feng_map: dict[str, int] = {
                antenna.decode(): index for antenna, index in f["antenna_feng_map"][()]}
            # Parse out the antenna positions
            antenna_positions: dict = {}
            for antenna_descriptor in f["antenna_positions"]:
                kp_antenna: Antenna = Antenna(antenna_descriptor.decode())
                if kp_antenna.name in antenna_feng_map.keys():
                    antenna_positions[kp_antenna.name] = kp_antenna
            # Parse out general metadata
            metadata: dict = dict(f.attrs)

            calibration_solutions = []
            if 'calibration_solutions' in f:
                for epoch in f['calibration_solutions']:
                    antenna_pols = [str(ap) for ap in f['calibration_solutions'][epoch]]
                    gains = np.array([np.array(f['calibration_solutions'][epoch][ap]) for ap in antenna_pols])

                    # Convert to AFP order, to match dada file axes
                    gains = gains.reshape(-1,2,metadata['nchans']).transpose(0,2,1)
                    S = CalibrationSolution(epoch, antenna_pols, gains)
                    calibration_solutions.append(S)

            return cls(
                antenna_positions,
                antenna_feng_map,
                metadata["cfreq"] * u.Hz,
                metadata["bandwidth"] * u.Hz,
                metadata["nchans"],
                metadata["project_id"],
                metadata["sb_id"],
                metadata["cbf_version"],
                Time(metadata["sync_epoch"], format="unix"),
                f["phase_centres"][()].astype([ # pylint: disable=no-member
                    ("timestamp", "datetime64[us]"), ("value", "|S64")]),
                f["suspect_flags"][()].astype([ # pylint: disable=no-member
                    ("timestamp", "datetime64[us]"), ("value", "bool")]),
                calibration_solutions
            )


    def _drop_duplicate_values(self, ar: np.ndarray) -> np.ndarray:
        vals = ar["value"]
        idxs = np.concatenate(([0], np.where(vals[1:] != vals[:-1])[0] + 1))
        return ar[idxs]

    def _covert_to_windows(self, ar: np.ndarray,
                           sentinel_time: Time) -> list[tuple[Time, Time, Any]]:
        windows = []
        for ii, (start, value) in enumerate(ar):
            if ii + 1 == len(ar):
                end = sentinel_time
            else:
                end, _ = ar[ii + 1]
            windows.append((start, end, value))
        return windows

    def get_subarray(self, antenna_names: list[str] = None) -> Subarray:
        """Get a subarray of the antenna set

        Args:
            antenna_names (list[str], optional): List of antenna names to extract. 
            Defaults to the full array.

        Returns:
            Subarray: _description_
        """
        if antenna_names is None:
            return Subarray(list(self.antenna_positions.values()))
        antennas = []
        for name in antenna_names:
            if name not in self.antenna_positions:
                raise KeyError(
                    f"Antenna {name} is not available in the array")
            else:
                antennas.append(self.antenna_positions[name])
        return Subarray(antennas)

    def find_observing_windows(
        self, min_duration: Quantity = 60 * u.s, allow_suspect: bool = False
    ) -> list[tuple[Time, Time, Target, bool]]:
        """Find observing windows corresponding to phase centre and suspect flags

        Args:
            min_duration  (Quantity, optional): Only return observing windows 
                                                bigger than this length.
            allow_suspect (bool, optional): Also return observing windows with suspect data. 
                                            Defaults to False.


        Returns:
            list: A list of observing epochs specifying the start, end, phase centre 
                  and suspect flag.

        Notes:
            The timing of the suspect flags to phase centre changes is unpredictable 
            and so to be safe it is necessary to only return valid observing windows 
            above a given size. Typically invalid windows are less than 1 second long, 
            but we leave the default as 60 second under the expectation that observations 
            will be considerably longer than this.
        """
        # Rather than handle some special case for the end of timestamp
        # sequences, we here define a sentinel time far away in the future.
        # This is chosen to be a year after the start of the observation
        # for saftey.
        sentinel_time = self.phase_centres[0]["timestamp"] + \
            np.timedelta64(365, 'D')
        phase_centre_windows = self._covert_to_windows(
            self._drop_duplicate_values(self.phase_centres),
            sentinel_time)
        suspect_windows = self._covert_to_windows(
            self._drop_duplicate_values(self.suspect_flags),
            sentinel_time)
        overlaps: list = []
        for pc_start, pc_end, target_descriptor in phase_centre_windows:
            target = Target(target_descriptor.decode())
            for flag_start, flag_end, is_suspect in suspect_windows:
                if not allow_suspect and is_suspect:
                    continue
                if pc_start < flag_end and pc_end > flag_start:
                    overlap_start = max(pc_start, flag_start)
                    overlap_end = min(pc_end, flag_end)
                    duration = (overlap_end - overlap_start) / \
                        np.timedelta64(1, 's')
                    if duration < min_duration.to(u.s).value:
                        continue
                    overlaps.append((
                        Time(overlap_start),
                        Time(overlap_end),
                        target,
                        is_suspect
                    ))
        return overlaps

    def get_pointings(self) -> list[PointingMetadata]:
        """Fetch all pointings from the session

        Returns:
            list[PointingMetadata]: A list of pointing objects

        Notes:
            The pointing objects are wrappers around the windows that
            are returned by find_observing_windows().
        """
        windows = self.find_observing_windows(allow_suspect=False)
        pointings = []
        for start, end, target, _ in windows:
            pointings.append(PointingMetadata(target, start, end, self))
        return pointings



@dataclass
class PointingMetadata:
    """Class for tracking information about a specific pointing
    """
    # Boresight phasing position for the pointing
    phase_centre: Target
    # The start time of the pointing
    start_epoch: Time
    # The end time of the pointing
    end_epoch: Time
    # Reference to the session metadata
    session_metadata: SessionMetadata


@dataclass
class BeamSet:
    """Class for tracking information about a beam set

    A beam set is a collection of beams which are formed from
    a common subarray.
    """
    name: str
    anntenna_names: list[str]
    beams: list[Target]
    tilings: list[dict]


@dataclass
class BeamformerConfig:
    """Configuration information for the beamformer pipeline
    """
    # The total number of beams
    nbeams: int
    # The tscrunch factor
    tscrunch: int
    # The fscrunch factor
    # Note fscrunching is applied after incoherent dedispersion
    fscrunch: int
    # Which Stokes vector should be extracted (I=0, Q=1, U=2, V=3)
    stokes_mode: int
    # Subtract the incoherent beam from the coherent beams
    subtract_ib: bool
    # List of dispersion measures to coherently dedisperse to
    coherent_dms: list[float]
    # List of beam sets (sets of beams from a common subarray)
    beam_sets: list[BeamSet]
    # Filter parameters
    filter: dict[str]

    @classmethod
    def from_file(cls, config_file: str) -> Self:
        """Read a beamformer config from a YAML file

        Args:
            config_file (str): Path to the configuration file

        Returns:
            Self: An instance of BeamformerConfig
        """
        with open(config_file, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        bfc = data["beamformer_config"]
        filter_config = data.get("filter", {"enabled": False})
        beam_sets = []
        for bs in data["beam_sets"]:
            beam_sets.append(
                BeamSet(
                    bs["name"],
                    bs["antenna_set"],
                    bs["beams"],
                    bs["tilings"],
                ))

        return cls(
            bfc["total_nbeams"],
            bfc["tscrunch"],
            bfc["fscrunch"],
            bfc["stokes_mode"],
            bfc["subtract_ib"],
            bfc["coherent_dms"],
            beam_sets,
            filter_config
        )

def make_tiling(
        pointing: PointingMetadata,
        subarray: Subarray,
        tiling_desc: dict) -> Tuple[Tiling, BeamShape]:
    """Make a tiling using the complete mosaic tiling options

    Args:
        pointing (PointingMetadata): Information on the pointing being tiled
        subarray (Subarray): The subarray to use for PSF determination
        tiling_desc (dict): The mosaic tiling arguments

    Returns:
        Tiling: A mosaic Tiling object containing tiled beam positions

    Note:
        The tiling description dictionary has been kept close to the form
        it is has in the FBFUSE configuration authority at MeerKAT. Examples
        are below.

        Simple circular tiling:
        {
        "nbeams": 22,
        "overlap": 0.7
        }

        Hexagonal tiling with specified frequency and epoch offset:
        {
            'nbeams': 123,
            'reference_frequency': 1400000000.0, # Hz
            'edelay': 30.0, # second from start of scan
            'target': 'source0,radec,00:00:00,00:00:00',
            'shape': 'hexagon',
            'method': 'variable_overlap',
            'shape_parameters': [0.16656555, 0.0],
            'coordinate_type': 'galactic'
        }

        The following defaults are used:
        reference_frequency --> centre frequency of the current observation
        edelay --> half the duration of the observation
        target --> phase centre of the pointing
        shape --> circular
        method --> variable size (e.g. force the overlap but allow
                   the total tiled area to change)
        overlap --> half power point overlaps
        coordinate_type --> equatorial coordinates
    """
    # Need to parse all the tiling information
    ref_freq = tiling_desc.get("reference_frequency", None)
    if ref_freq is None:
        ref_freq = pointing.session_metadata.centre_frequency.to(u.Hz).value
    target = tiling_desc.get("target", None)
    if target is None:
        target = pointing.phase_centre
    else:
        target = Target(target)
    edelay = tiling_desc.get("edelay", None)
    if edelay is None:
        epoch = pointing.start_epoch + (
            pointing.end_epoch - pointing.start_epoch) / 2
    else:
        epoch = pointing.start_epoch + TimeDelta(edelay * u.s)
    shape = tiling_desc.get("shape", "circle")
    shape_params = tiling_desc.get("shape_parameters", None)
    coordinate_type = tiling_desc.get("coordinate_type", "equatorial")
    method = tiling_desc.get("method", "variable_size")
    nbeams = tiling_desc["nbeams"]
    overlap = tiling_desc.get("overlap", 0.5)

    antenna_strings: list[str] = [
        ant.format_katcp() for ant in subarray.antenna_positions
    ]
  
    psfsim: PsfSim = PsfSim(antenna_strings, ref_freq)
    psf_beam_shape: BeamShape = psfsim.get_beam_shape(target, epoch.unix, beam_number=400)
    #Build Mosaic command here. Remove T
    mosaic_epoch = epoch.iso.replace("T", " ")
    mosaic_epoch = mosaic_epoch.replace("-", ".")
    mosaic_antenna_string = ','.join([item.replace('m', '') for item in subarray.names])
    mosaic_command=f"python maketiling.py --freq {ref_freq} --source {target.body._ra} {target.body._dec} --datetime {mosaic_epoch} --subarray {mosaic_antenna_string} --verbose --tiling_method {method} --tiling_shape {shape} --ants antenna.csv --beamnum {nbeams} --overlap {overlap}"
    

    tiling: Tiling = generate_nbeams_tiling(
        psf_beam_shape, nbeams, overlap,
        method, shape,
        parameter=shape_params,
        coordinate_type=coordinate_type)
    return tiling, psf_beam_shape, mosaic_command

def pad_ra_dec(ra, dec):
    # Split RA and Dec into components
    ra_parts = ra.split(':')
    dec_parts = dec.split(':')
    
    # Pad RA hours and Dec degrees to two digits
    ra_parts[0] = ra_parts[0].zfill(2)  # Ensure two digits for RA hours
    dec_parts[0] = dec_parts[0].zfill(2)  # Ensure two digits for Dec degrees
    
    # Join the parts back together
    padded_ra = ':'.join(ra_parts)
    padded_dec = ':'.join(dec_parts)
    
    return padded_ra, padded_dec

class PsfInterpolator:
    """Handles interpolation of a PSF image onto arbitrary (x,y) coordinate grids."""
    def __init__(self, psf, boresight):
        self._psf = psf
        self._cdelt = psf.wcs_header['cdelt']
        image = self._psf.image

        x = self._calc_offsets(image.shape[0], self._cdelt[0]) + boresight[0] - self._cdelt[0] / 2
        y = self._calc_offsets(image.shape[1], self._cdelt[1]) - boresight[1] - self._cdelt[1] / 2

        self._interpolator = RegularGridInterpolator(
            (y, x), image, bounds_error=False, fill_value=-1
        )

    @staticmethod
    def _calc_offsets(npix, width):
        half = (npix * width) / 2
        return np.linspace(-half + width/2, half - width/2, npix)

    def __call__(self, pts):
        return self._interpolator(pts)
    

def initialize_beam_table(tilings):
    dtype = [
        ("tiling_idx", "uint32"),
        ("tiling_beam_idx", "uint32"),
        ("global_beam_idx", "uint32"),
        ("coordinates", [("x", "float32"), ("y", "float32")]),
        ("count", "uint32"),
        ("relative_count", "float32"),
        ("keep", "bool")
    ]
    total_beams = sum([t.beam_num for _, t in tilings])
    table = np.zeros(total_beams, dtype=dtype)
    n = 0

    for tiling_idx, (bs, tiling) in enumerate(tilings):
        for beam_idx, coord in enumerate(tiling.coordinates):
            table[n] = (
                tiling_idx,
                beam_idx,
                n,
                (coord[0], coord[1]),
                0,
                0.0,
                True
            )
            n += 1

    return table


def generate_grid(beam_table, resolution, margin_fraction=0.1):
    def generate_axis(coords):
        margin = abs((coords.max() - coords.min())) * margin_fraction 
        lower_bound = coords.min() - margin
        upper_bound = coords.max() + margin
        npix = int((upper_bound - lower_bound) / resolution)
        return np.linspace(lower_bound, upper_bound, npix)
    x = generate_axis(beam_table["coordinates"]["x"])
    y = generate_axis(beam_table["coordinates"]["y"])
    X, Y = np.meshgrid(x, y)
    grid_coords = np.stack((X, Y), axis=-1)
    return x, y, grid_coords


def compute_gain_and_mask(tilings, all_beams, grid_coords, gain_cap, gain_tol):
    gain = np.zeros(grid_coords.shape[:2], dtype="float32")
    mask = np.full(grid_coords.shape[:2], -1, dtype="int32")

    for tiling_idx, (bs, tiling) in enumerate(tilings):
        array_gain = len(bs.anntenna_names)
        tiling_data = all_beams[all_beams["tiling_idx"] == tiling_idx]

        for row in tiling_data:
            interp = PsfInterpolator(tiling.beam_shape.psf, row["coordinates"])
            gmap = interp(grid_coords)
            gmap[gmap < gain_cap] = -1.0

            val = gmap * array_gain
            update = val >= (gain_tol * gain)
            mask[update] = row["global_beam_idx"]
            gain = np.maximum(gain, val)

    return gain, mask


def update_beam_counts(tilings, all_beams, mask, low_influence_tol):
    for tiling_idx, (bs, tiling) in enumerate(tilings):
        tiling_data = all_beams[all_beams["tiling_idx"] == tiling_idx]
        lo, hi = tiling_data["global_beam_idx"].min(), tiling_data["global_beam_idx"].max()

        valid = (mask >= lo) & (mask <= hi)
        unique, counts = np.unique(mask[valid], return_counts=True)

        m = np.median(counts)
        all_beams["count"][unique] = counts
        all_beams["relative_count"][unique] = counts / m

    all_beams["keep"] = all_beams["relative_count"] >= low_influence_tol

def plot_tiling_influence(tiling_idx, tiling_data, pdf=None):
    fig, ax = plt.subplots(figsize=(6, 5))

    counts = tiling_data["count"]
    bad = tiling_data[tiling_data["keep"] == False]

    norm = SymLogNorm(
        linthresh=0.9,
        linscale=1.0,
        vmin=counts.min(),
        vmax=counts.max()
    )

    sc = ax.scatter(
        tiling_data["coordinates"]["y"],
        tiling_data["coordinates"]["x"],
        c=counts,
        cmap="Wistia",
        norm=norm
    )

    ax.scatter(
        bad["coordinates"]["y"],
        bad["coordinates"]["x"],
        facecolor="none",
        edgecolor="k",
        lw=1.5,
        label="Excised Beams"
    )

    plt.colorbar(sc, label="Beam Pixel Count")

    ax.set_title(f"Tiling {tiling_idx}: Influence Counts")
    ax.set_xlabel("Y Offset (deg)")
    ax.set_ylabel("X Offset (deg)")
    ax.legend()
    plt.tight_layout()
    if pdf:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

def plot_gain_and_mask(gain, mask, x, y, pdf=None, suffix=""):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    if suffix:
        suffix = f" ({suffix})"

    ax[0].set_title("Peak Gain Map" + suffix)
    im0 = ax[0].imshow(gain, extent=[x.min(), x.max(), y.min(), y.max()], 
                       origin="lower", cmap="hot")
    plt.colorbar(im0, ax=ax[0], label="Gain")

    ax[1].set_title("Beam Dominance Mask" + suffix)
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    im1 = ax[1].imshow(mask, extent=[x.min(), x.max(), y.min(), y.max()],
                       cmap="gist_ncar", origin="lower")
    plt.colorbar(im1, cax=cax, label="Beam Index")

    for axis in ax:
        axis.set_xlabel("Offset (deg)")
        axis.set_ylabel("Offset (deg)")

    plt.tight_layout()
    if pdf:
        pdf.savefig()
        plt.close()
    else:
        plt.show()

def filter_tilings(beamformer_config, bs_tilings, outfile):
    log.info("Filtering generated beam sets")
    bc = beamformer_config
    tilings = []
    for bs, subtilings in zip(bc.beam_sets, bs_tilings):
        for tiling in subtilings:
            tilings.append((bs, tiling))
    all_beams = initialize_beam_table(tilings)
    x, y, grid_coords = generate_grid(all_beams, bc.filter["resolution"], margin_fraction=bc.filter["margin"])
    log.debug("Computing gain and beam dominance maps")
    gain, mask = compute_gain_and_mask(tilings, all_beams, grid_coords, bc.filter["gain_cap"], bc.filter["gain_tol"])
    log.debug("Updating beam influence counts")
    update_beam_counts(tilings, all_beams, mask, bc.filter["low_influence_tol"])
    log.debug("Generating output plots")
    with PdfPages(f"{outfile}_filter.pdf") as pdf:
        # High-level plot
        plot_gain_and_mask(gain, mask, x, y, pdf, "pre")
        
        # Per-tiling influence plots
        for tiling_idx, _ in enumerate(tilings):
            tiling_data = all_beams[all_beams["tiling_idx"] == tiling_idx]
            plot_tiling_influence(tiling_idx, tiling_data, pdf)
    
        # Modify tilings in-place
        for tiling_idx, (bs, tiling) in enumerate(tilings):
            keepers = all_beams[(all_beams["tiling_idx"] == tiling_idx) & all_beams["keep"]]
            tiling.coordinates = keepers["coordinates"].tolist()
            tiling.num_beams = len(tiling.coordinates)
        
        all_beams_post = initialize_beam_table(tilings)
        gain_post, mask_post = compute_gain_and_mask(tilings, all_beams_post, grid_coords, 
                                                     bc.filter["gain_cap"], bc.filter["gain_tol"])

        # High-level plots
        plot_gain_and_mask(gain_post, mask_post, x, y, pdf, "post")
    
    return bs_tilings

def create_delays(
        session_metadata: SessionMetadata,
        beamformer_config: BeamformerConfig,
        pointing: PointingMetadata,
        start_epoch: Time = None,
        end_epoch: Time = None,
        step: TimeDelta = 4 * u.s,
        outfile: str = None) -> list[DelayModel]:
    """Create a set of delay models

    Args:
        session_metadata (SessionMetadata): A session metadata object
        beamformer_config (BeamformerConfig): A beamformer configuraiton object
        pointing (PointingMetadata): A pointing metadata object
        start_epoch (Time, optional): The start of the window to produce delays for.
                                      Defaults to the start of the pointing.
        end_epoch (Time, optional): The end of the window to produce delays until.
                                    Defaults to the end of the pointing.
        step (TimeDelta, optional): The step size between consequtive solutions.
                                    Defaults to 4*u.s.
        outfile (str, optional): The path to write the delay models to. Defaults is None.

    Returns:
        list[DelayModel]: A list of delay models
    """
    if start_epoch is None:
        start_epoch = pointing.start_epoch
    if end_epoch is None:
        end_epoch = pointing.end_epoch
    om = session_metadata
    bc = beamformer_config
    log.info("Creating delays for target %s", pointing.phase_centre.name)
    log.info("Start epoch: %s (UNIX %f)", start_epoch.isot, start_epoch.unix)
    log.info("End epoch: %s (UNIX %f)", end_epoch.isot, end_epoch.unix)
    log.info("Step size: %s", step.to(u.s))
    full_subarray = om.get_subarray()
    de = DelayEngine(full_subarray, pointing.phase_centre)
    # Initialize beam_set_lookup dictionary to track unique beam sets
    beam_set_lookup = {}
    beam_set_id = 0
    plot_beams = []
    neighbouring_beams = []
    #Initialise the known beams. PSF size is calculated for 50% overlap
    nbeams_requested = 1
    overlap = 0.5
    #Iterate beams first through all beam sets
    for bs in bc.beam_sets:
        sorted_antennas = sorted(bs.anntenna_names)
        subarray_subset = om.get_subarray(sorted_antennas)               
        antenna_string = ','.join(sorted_antennas)
        if bs.beams is not None:
            for target_desc in bs.beams:
                target = Target(target_desc)
                #Add the beam to the delay engine
                de.add_beam(target, subarray_subset)
                beam_key = (antenna_string, overlap, nbeams_requested)
                if beam_key not in beam_set_lookup:
                    beam_set_lookup[beam_key] = beam_set_id
                    beam_set_id += 1
                ra = str(target.body._ra)
                dec = str(target.body._dec)
                # Pad RA and Dec to two digits
                ra, dec = pad_ra_dec(ra, dec)
                name = target.name
                current_beam_set_id = beam_set_lookup[beam_key]
                #Get the PSF Beam shape at 50 % overlap for beams defined in yaml.
                tiling_desc = {
                    "nbeams": nbeams_requested,
                    "overlap": overlap,
                    "target": target_desc,
                }
                _, psf_beam_shape,_ = make_tiling(pointing, subarray_subset, tiling_desc)
                #Add the beam to the plot beams list
                plot_beams.append((name, ra, dec, round(psf_beam_shape.axisH, 5), round(psf_beam_shape.axisV, 5), round(psf_beam_shape.angle, 5), current_beam_set_id, 0.5, len(sorted_antennas), 'known'))
                neighbouring_beams.append((name, ra, dec, round(psf_beam_shape.axisH, 5), round(psf_beam_shape.axisV, 5), round(psf_beam_shape.angle, 5), current_beam_set_id, 0.5, len(sorted_antennas), 'known'))
    
    #Iterate through all tilings
    bs_tilings = []
    bs_list = []
    for bs in bc.beam_sets:
        tilings = []
        sorted_antennas = sorted(bs.anntenna_names)
        subarray_subset = om.get_subarray(sorted_antennas)
        antenna_string = ','.join(sorted_antennas)
        if bs.tilings is not None:
            output_prefix = f"{outfile}_{bs.name}" if outfile is not None else None
            for tiling_desc in bs.tilings:
                tiling, psf_beamshape, mosaic_command = make_tiling(pointing, subarray_subset, tiling_desc)
                #Add the tiling to the delay engine
                de.add_tiling(tiling, subarray_subset)
                tilings.append(tiling)
        bs_tilings.append(tilings)
    
    if bc.filter.get("enable", False):
        bs_tilings = filter_tilings(bc, bs_tilings, outfile)
     
    for bs, tilings in zip(bc.beam_sets, bs_tilings):
        sorted_antennas = sorted(bs.anntenna_names)
        subarray_subset = om.get_subarray(sorted_antennas)
        antenna_string = ','.join(sorted_antennas)
        if bs.tilings is not None:
            output_prefix = f"{outfile}_{bs.name}" if outfile is not None else None
            for tiling in tilings:
                _, psf_beamshape, mosaic_command = make_tiling(pointing, subarray_subset, tiling_desc)
                cb_beamshape = tiling.meta["axis"][:3] #axisH, axisV, angle
                overlap = tiling.meta["axis"][-1]
                nbeams_requested = tiling_desc['nbeams']
                beam_key = (antenna_string, overlap, nbeams_requested)
                if beam_key not in beam_set_lookup:
                    beam_set_lookup[beam_key] = beam_set_id
                    beam_set_id += 1
                current_beam_set_id = beam_set_lookup[beam_key]
                coords = tiling.get_equatorial_coordinates()
                coords = SkyCoord(coords, unit=u.deg)
                mosaic_command+=f" --tiling_plot {output_prefix}_bid_{beam_set_id}.png --tiling_coordinate {output_prefix}_bid_{beam_set_id}.csv"
                log.info("Mosaic command written to %s", f"{outfile}.mosaic")
                log.info(f"Writing PSF of BeamSet {bs.name} to {output_prefix}_bid_{beam_set_id}.fits")
                psf_beamshape.psf.write_fits(f"{output_prefix}_bid_{beam_set_id}.fits")
                log.info(f"Writing PSF Plot of BeamSet {bs.name} to {output_prefix}_bid_{beam_set_id}.png")
                psf_beamshape.plot_psf(f"{output_prefix}_bid_{beam_set_id}.png")
                with open(f"{outfile}.mosaic", "a") as f:
                    f.write(mosaic_command + "\n")
                for index, coord in enumerate(coords):
                    ra_hms = coord.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
                    dec_dms = coord.dec.to_string(unit=u.degree, sep=':', precision=1, alwayssign=True, pad=True)
                    plot_beams.append((f"{bs.name}_{index:03d}", ra_hms, dec_dms, round(cb_beamshape[0], 5), round(cb_beamshape[1], 5), round(cb_beamshape[2], 5), current_beam_set_id, overlap, len(sorted_antennas), 'tiling'))
                    neighbouring_beams.append((f"{bs.name}_{index:03d}", ra_hms, dec_dms, round(psf_beamshape.axisH, 5), round(psf_beamshape.axisV, 5), round(psf_beamshape.angle, 5), current_beam_set_id, 0.5, len(sorted_antennas), 'tiling'))
    
    # What is this for???    
    target = tiling_desc.get("target", None)
    if target is None:
        print("No target specified")
        target = pointing.phase_centre
    else:
        target = Target(target)
    column_list = ['name', 'ra', 'dec', 'x', 'y', 'angle', 'beam_set_id', 'overlap', 'nantennas', 'type']
    plot_beams_df = pd.DataFrame(plot_beams, columns=column_list)
    #Plot beams has beam shape for the overlap requested, wheras neighbouring beams has the PSF shape (50% overlap)
    neighbouring_beams_df = pd.DataFrame(neighbouring_beams, columns=column_list)
    log.info("Calculating solutions for %d antennas and %d beams", full_subarray.nantennas, de.nbeams)
    delays = de.calculate_delays(start_epoch, end_epoch, step)
    log.info("Beams and tilings written to %s.targets", outfile)
    plot_beams_df.to_csv(outfile + ".targets", index=False)
    boresight_coords = psf_beamshape.bore_sight.equatorial
    plot_multiple_tilings(pointing, neighbouring_beams_df, plot_beams_df, boresight_coords, outfile)
    return delays, de.targets, bs_tilings

def plot_multiple_tilings(pointing, neighbouring_beams_df, plot_beams_df, boresight_coords, outfile, HD=True, beam_size_scaling=1.0, annotate_beam_names=False):
    # Initialize WCS projection
    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [0, 0]
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs_properties.wcs.crval = boresight_coords
    center = boresight_coords
    
    # Use neutral scaling for overlap detection
    wcs_properties.wcs.cdelt = [1, 1]  # Neutral pixel scale

    detector = BeamOverlapDetector(neighbouring_beams_df, wcs_properties)
    overlaps_df = detector.find_overlapping_beams()
    overlaps_df.to_csv(outfile + "_overlapping_beams.csv", index=False)
    
    # Update WCS projection for plotting
    step = 1 / 10000000000.0
    wcs_properties.wcs.cdelt = [-step, step]
    resolution = step

    thisDPI = 300
    if HD:
        width = 3200.
        extra_source_text_size = 8
    else:
        width = 800.
        extra_source_text_size = 3

    fig = plt.figure(figsize=(width/thisDPI, width/thisDPI), dpi=thisDPI)
    axis = fig.add_subplot(111, aspect='equal', projection=wcs_properties)
    
    # Define color palette
    color_palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']

    # Extract relevant data from DataFrame in one step for efficiency and clarity
    beam_ra, beam_dec, beam_name, beam_set_id, beam_x, beam_y, beam_angle, beam_type, nantennas, overlap = (
        plot_beams_df['ra'].astype(str),
        plot_beams_df['dec'].astype(str),
        plot_beams_df['name'].astype(str),
        plot_beams_df['beam_set_id'].astype(int),
        plot_beams_df['x'].astype(float),
        plot_beams_df['y'].astype(float),
        plot_beams_df['angle'].astype(float),
        plot_beams_df['type'].astype(str),
        plot_beams_df['nantennas'].astype(int),
        plot_beams_df['overlap'].astype(float)
    )
   
    # Get equatorial coordinates
    equatorialCoordinates = SkyCoord(beam_ra, beam_dec, frame='fk5', unit=(u.hourangle, u.deg))
    beam_coordinate = np.array(wcs_properties.wcs_world2pix(np.array([equatorialCoordinates.ra.deg, equatorialCoordinates.dec.deg]).T, 0))

    # Plot the boresight for reference
    axis.plot(0, 0, marker='+', markersize=15, color='black')

    # Store labels for legend creation based on beam_set_id
    labels = {}

    # Loop through all beams to create ellipses and prepare legend text
    for idx in range(len(beam_coordinate)):
        coord = beam_coordinate[idx]

        # Create the ellipse for each beam
        ellipse = Ellipse(
            xy=coord,
            width=2.0 * beam_x.iloc[idx] * beam_size_scaling / resolution,
            height=2.0 * beam_y.iloc[idx] * beam_size_scaling / resolution,
            angle=beam_angle.iloc[idx]
        )
        ellipse.fill = False
        color = color_palette[beam_set_id.iloc[idx] % len(color_palette)]
        ellipse.set_edgecolor(color)
        axis.add_artist(ellipse)

        # Prepare legend text based on the beam type criteria
        if beam_type.iloc[idx] == 'known':
            legend_text = f"Beam set{beam_set_id.iloc[idx]}: {beam_type.iloc[idx]}, {nantennas.iloc[idx]} antennas"
            axis.annotate(
                beam_name.iloc[idx],  
                xy=coord,            
                xytext=(coord[0] + 5, coord[1] + 5),  # Text position, slightly offset
                fontsize=8 
            )
        else:
            legend_text = f"Beam set{beam_set_id.iloc[idx]}: {beam_name.iloc[idx].split('_')[0]}, {nantennas.iloc[idx]} antennas, {overlap.iloc[idx]} overlap"
            # Add an annotation with the beam name. Be careful with the font size, as it can be too large for the plot
            if annotate_beam_names:
                axis.annotate(beam_name.iloc[idx], xy=coord, xytext=(coord[0], coord[1]), fontsize=extra_source_text_size)
        # If the legend entry for this beam_set_id does not exist yet, create it
        if beam_set_id.iloc[idx] not in labels:
            labels[beam_set_id.iloc[idx]] = (ellipse, legend_text)

    # Create a legend for all beam categories (all types)
    handles = [labels[k][0] for k in labels]
    legend_texts = [labels[k][1] for k in labels]
    axis.legend(handles, legend_texts, loc='upper right')
    margin = 1.1 * max(np.sqrt(np.sum(np.square(beam_coordinate), axis=1)))
    axis.set_xlim(center[0] - margin, center[0] + margin)
    axis.set_ylim(center[1] - margin, center[1] + margin)
    axis.set_xlabel('RA', fontsize=30)
    axis.set_ylabel('Dec', fontsize=30)

    # Save the original, unzoomed plot
    # output_filename_unzoomed = f"{outfile}_original.png"
    # log.info(f"Saving original unzoomed plot to {output_filename_unzoomed}")
    # plt.tight_layout()
    # plt.savefig(output_filename_unzoomed)

    # Now zoom into the tiling beam sets and save separate plots
    tiling_df = plot_beams_df[plot_beams_df['type'] == 'tiling']

    for beam_set, group_df in tiling_df.groupby('beam_set_id'):
        equatorialCoordinates_tiling = SkyCoord(group_df['ra'].astype(str), group_df['dec'].astype(str), frame='fk5', unit=(u.hourangle, u.deg))
        beam_coordinate_tiling = np.array(wcs_properties.wcs_world2pix(np.array([equatorialCoordinates_tiling.ra.deg, equatorialCoordinates_tiling.dec.deg]).T, 0))
        # Recalculate margin for this specific beam set
        margin = 1.3 * max(np.sqrt(np.sum(np.square(beam_coordinate_tiling), axis=1)))

        # Set the new axis limits for zooming into this region
        axis.set_xlim(center[0] - margin, center[0] + margin)
        axis.set_ylim(center[1] - margin, center[1] + margin)

        # Save the zoomed plot for the current beam_set_id
        output_filename_zoomed = f"{outfile}_tiling_beamset_{beam_set}.png"
        log.info(f"Saving zoomed plot for beam_set_id {beam_set} to {output_filename_zoomed}")
        plt.title(f"Boresight: {pointing.phase_centre.name}, UTC Start {pointing.start_epoch.isot}, Zoomed into Beam set {beam_set}")
        plt.tight_layout()
        plt.savefig(output_filename_zoomed)


class BeamOverlapDetector:
    def __init__(self, neighbour_df, wcs_properties):
        self.neighbour_df = neighbour_df
        self.wcs_properties = wcs_properties
        self._convert_coordinates()
        self._convert_to_pixel_coordinates()

    def _convert_coordinates(self):
        # Convert RA and Dec from strings to degrees
        coords = SkyCoord(ra=self.neighbour_df['ra'].values, dec=self.neighbour_df['dec'].values, unit=(u.hourangle, u.deg))
        self.neighbour_df['ra_deg'] = coords.ra.deg
        self.neighbour_df['dec_deg'] = coords.dec.deg
        
    def _convert_to_pixel_coordinates(self):
        # Convert RA and Dec to pixel coordinates using WCS
        sky_coords = SkyCoord(ra=self.neighbour_df['ra_deg'].values, dec=self.neighbour_df['dec_deg'].values, unit='deg')
        pixel_coords = np.array(self.wcs_properties.wcs_world2pix(np.array([sky_coords.ra.deg, sky_coords.dec.deg]).T, 0))
        self.neighbour_df['x_pix'], self.neighbour_df['y_pix'] = pixel_coords[:, 0], pixel_coords[:, 1]
        

    def ellipse_parametric(self, t, a, b, x0, y0, theta):
        cos_t = np.cos(t)
        sin_t = np.sin(t)
        x = x0 + a * cos_t * np.cos(theta) - b * sin_t * np.sin(theta)
        y = y0 + a * cos_t * np.sin(theta) + b * sin_t * np.cos(theta)
        return x, y

    def point_in_ellipse(self, x, y, ellipse):
        x0, y0, a, b, theta = ellipse["x0"], ellipse["y0"], ellipse["a"], ellipse["b"], ellipse["theta"]
        cos_theta = np.cos(-theta)
        sin_theta = np.sin(-theta)
        xr = cos_theta * (x - x0) - sin_theta * (y - y0)
        yr = sin_theta * (x - x0) + cos_theta * (y - y0)
        return (xr**2 / a**2) + (yr**2 / b**2) <= 1
    
    def check_containment(self, ellipse1, ellipse2):
        # Check if the center of ellipse1 is inside ellipse2
        x0, y0 = ellipse1["x0"], ellipse1["y0"]
        return self.point_in_ellipse(x0, y0, ellipse2)

    def discrete_overlap(self, ellipse1, ellipse2, num_points=100):
        # First, check if one ellipse contains the center of the other
        if self.check_containment(ellipse1, ellipse2) or self.check_containment(ellipse2, ellipse1):
            return True

        # Check points on the perimeters of both ellipses
        t_values = np.linspace(0, 2 * np.pi, num_points)
        #psf_x, psf_y -> semi-major axis (a), semi-minor axis (b) of the ellipse
        x1, y1 = self.ellipse_parametric(t_values, ellipse1["a"], ellipse1["b"], ellipse1["x0"], ellipse1["y0"], ellipse1["theta"])
        
        # Check if any points on the perimeter of ellipse1 lie inside ellipse2
        for x, y in zip(x1, y1):
            if self.point_in_ellipse(x, y, ellipse2):
                return True    
        return False
    
    def find_nearest_neighbours(self, beam, beams, n=6):
        # Create a SkyCoord object for the target beam
        target_coord = SkyCoord(ra=beam['ra_deg'] * u.deg, dec=beam['dec_deg'] * u.deg)

        # Create SkyCoord objects for all other beams
        all_coords = SkyCoord(ra=[b['ra_deg'] for b in beams] * u.deg,
                              dec=[b['dec_deg'] for b in beams] * u.deg)

        # Calculate separations
        separations = target_coord.separation(all_coords)
        
        # Sort by separation and get the nearest neighbours
        nearest_indices = np.argsort(separations)[1:n+1]  # Skip the first one (itself)

        # Return the list of nearest beam names in descending order of separation
        nearest_neighbours = [beams[i]['name'] for i in nearest_indices]
        return nearest_neighbours

    def find_overlapping_beams(self):
        beams = []
        for _, row in self.neighbour_df.iterrows():
            beams.append({
                "name": row['name'],
                "x0": row['x_pix'],
                "y0": row['y_pix'],
                "a": row['x'],
                "b": row['y'],
                "theta": np.radians(row['angle']),
                "beam_set_id": row['beam_set_id'],
                "ra_deg": row['ra_deg'],
                "dec_deg": row['dec_deg']
            })
        beam_overlap_dict = {}
        for i, beam1 in enumerate(beams):
            # Initialize overlap list for beam1 if not already in the dict
            if beam1['name'] not in beam_overlap_dict:
                beam_overlap_dict[beam1['name']] = {
                    "name": beam1['name'],
                    "x_pix": beam1['x0'],
                    "y_pix": beam1['y0'],
                    "a": beam1['a'],
                    "b": beam1['b'],
                    "angle": np.degrees(beam1['theta']),
                    "beam_set_id": beam1['beam_set_id'],
                    "overlapping_beams": [],  # Empty list initially
                    "neighbouring_beams": []  # Empty list for nearest neighbours
                }
            
            # Inner loop to check overlaps
            for j in range(i + 1, len(beams)):  # Only check pairs once
                beam2 = beams[j]

                # Avoid checking beams with the same beam_set_id. They have same PSF shape
                if beam1['beam_set_id'] == beam2['beam_set_id']:
                    continue

                # Initialize overlap list for beam2 if not already in the dict
                if beam2['name'] not in beam_overlap_dict:
                    beam_overlap_dict[beam2['name']] = {
                        "name": beam2['name'],
                        "x_pix": beam2['x0'],
                        "y_pix": beam2['y0'],
                        "a": beam2['a'],
                        "b": beam2['b'],
                        "angle": np.degrees(beam2['theta']),
                        "beam_set_id": beam2['beam_set_id'],
                        "overlapping_beams": [],  # Empty list initially
                        "neighbouring_beams": []  # Empty list for nearest neighbours
                    }
                
                # Check if beam1 and beam2 overlap
                if self.discrete_overlap(beam1, beam2):  
                    # Update overlap lists for both beams
                    beam_overlap_dict[beam1['name']]['overlapping_beams'].append(beam2['name'])
                    beam_overlap_dict[beam2['name']]['overlapping_beams'].append(beam1['name'])
            
            # Find the 6 nearest neighbours for beam1 only if beam set IDs match
            matching_beams = [b for b in beams if b['beam_set_id'] == beam1['beam_set_id']]
            nearest_neighbours = self.find_nearest_neighbours(beam1, matching_beams)
            beam_overlap_dict[beam1['name']]['neighbouring_beams'] = nearest_neighbours

        # Convert the dictionary to a DataFrame
        overlap_results = pd.DataFrame(beam_overlap_dict.values())
        overlap_results = overlap_results.sort_values(by=['beam_set_id', 'name'])

        return overlap_results


    




    
def main():
    """
    What does this thing actually do?

    1. Read the observation metadata file - Done
    2. Read some kind of beamformer configuration file
    4. From the beamformer configuration: - Done
        - generate tilings
        - generate beam pointings
        - generate antenna to beam maps
    5. Using the antenna database and the metadata: - Done
        - Generate all delays required for the full observation
          and write them to a (temporary?) file
        - Calculate the PSFs for each beam as a function of time
        - Generate a sky coverage map for the run
    6. Write a configuration header for the beamformer specifying the configuration
    7. Compile the beamformer for the specific configuration
    8. Excute the beamformer passing the weights file a pipeline configuration and the input data
    9. ...
    10. Profit


    Need to define:
    - Configuration file format (preferrably YAML)
    - Antenna database file format (preferrably YAML)
    - Delay format (preferrably some packed binary format)
    """

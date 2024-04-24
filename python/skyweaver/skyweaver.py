"""
Delay generation and metadata parsing for MeerKAT baseband recordings
"""
from __future__ import annotations

# stdlib imports
import logging
import textwrap
import ctypes
from typing import Self, Any
from dataclasses import dataclass

# 3rd party imports
import h5py
import yaml
import numpy as np
from rich.progress import track
from astropy.time import Time, TimeDelta
from astropy import units as u
from astropy.units import Quantity
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
    # ITRF reference position of the observatory
    itrf_reference: tuple[float, float, float]
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
            return cls(
                antenna_positions,
                antenna_feng_map,
                metadata["cfreq"] * u.Hz,
                metadata["bandwidth"] * u.Hz,
                tuple(metadata["itrf_reference"]),
                metadata["nchans"],
                metadata["project_id"],
                metadata["sb_id"],
                metadata["cbf_version"],
                Time(metadata["sync_epoch"], format="unix"),
                f["phase_centres"][()].astype([ # pylint: disable=no-member
                    ("timestamp", "datetime64[us]"), ("value", "|S64")]),
                f["suspect_flags"][()].astype([ # pylint: disable=no-member
                    ("timestamp", "datetime64[us]"), ("value", "bool")])
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
                    if overlap_end == sentinel_time:
                        # We convert sentinel time to None to remove any ambiguity
                        # a value of None implies "to the end of the session"
                        overlap_end = None
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
        beam_sets = []
        for bs in data["beam_sets"]:
            beam_sets.append(
                BeamSet(
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
            beam_sets
        )


def make_tiling(
        pointing: PointingMetadata,
        subarray: Subarray,
        tiling_desc: dict) -> Tiling:
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
    beam_shape: BeamShape = psfsim.get_beam_shape(target, epoch.unix)
    tiling: Tiling = generate_nbeams_tiling(
        beam_shape, nbeams, overlap,
        method, shape,
        parameter=shape_params,
        coordinate_type=coordinate_type)
    return tiling


def create_delays(
        session_metadata: SessionMetadata,
        beamformer_config: BeamformerConfig,
        pointing: PointingMetadata,
        start_epoch: Time = None,
        end_epoch: Time = None,
        step: TimeDelta = 4 * u.s) -> list[DelayModel]:
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
    bs_tilings = []
    for bs in bc.beam_sets:
        tilings = []
        subarray_subset = om.get_subarray(bs.anntenna_names)
        # add beams first as these are likely more important than tilings
        if bs.beams is not None:
            for target_desc in bs.beams:
                de.add_beam(Target(target_desc), subarray_subset)
        if bs.tilings is not None:
            for tiling_desc in bs.tilings:
                tiling = make_tiling(pointing, subarray_subset, tiling_desc)
                de.add_tiling(tiling, subarray_subset)
                tilings.append(tiling)
        bs_tilings.append(tilings)
    log.info(
        "Calculating solutions for %d antennas and %d beams", 
        full_subarray.nantennas, de.nbeams)
    delays = de.calculate_delays(start_epoch, end_epoch, step)
    return delays, de.targets, bs_tilings


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
    
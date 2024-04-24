"""
Command line parsing for the skyweaver application
"""
import argparse
import logging
import sys
import pprint

import secrets

import astropy.units as u
import skyweaver

log = logging.getLogger("skyweaver")

"""
What kind of cli do we want:
app: skyweaver
alias: sw

sw metadata show <metafile>
sw metadata show --key=<key> <metafile>
sw delays create --outfile <outfile> <metafile> <beamformer_config>
sw delays show <delayfile>
sw beamformer run --log-file=<logfile> --log-level=<LEVEL> <delayfile> <metafile> <dataconfig>
sw beamformer monitor <>

sw bfconfig create --template <blah> --overides=??? --metafile <> --window=4

We need to be able to do checkpointed recovery somehow aswell. This can probably be 
done per input file.
"""
class ColouredLogFormatter(logging.Formatter):
    """Colourise logging messages
    """
    grey = '\x1b[38;21m'
    blue = '\x1b[38;5;39m'
    yellow = '\x1b[38;5;226m'
    red = '\x1b[38;5;196m'
    bold_red = '\x1b[31;1m'
    reset = '\x1b[0m'

    def __init__(self, fmt: str) -> None:
        """Create a formatter

        Args:
            fmt (str): Formatting string for log messages
        """
        super().__init__()
        self.fmt = fmt
        self._formats = {
            logging.DEBUG: self.grey + self.fmt + self.reset,
            logging.INFO: self.blue + self.fmt + self.reset,
            logging.WARNING: self.yellow + self.fmt + self.reset,
            logging.ERROR: self.red + self.fmt + self.reset,
            logging.CRITICAL: self.bold_red + self.fmt + self.reset
        }

    def format(self, record: logging.LogRecord) -> logging.LogRecord:
        """Format a log record

        Args:
            record (logging.LogRecord): A log record to format

        Returns:
            logging.LogRecord: A formatted log record
        """
        log_fmt = self._formats.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

def init_logger(log_level: str, log_file: str = None) -> None:
    """Initialise skyweaver logging.

    Args:
        log_level (str): The severity threshold for logging (e.g. DEBUG, INFO, etc.)
        log_file (str, optional): A file to log to. Defaults to no file logging.
    """
    root_logger = logging.getLogger()
    if root_logger.hasHandlers():
        root_logger.handlers.clear()
    logger = logging.getLogger("skyweaver")
    # Set the default filter level for logging
    logger.setLevel(log_level.upper())
    # Create coloured logging formatter
    formatter = ColouredLogFormatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # Always setup a stdout handler
    stdout_handler = logging.StreamHandler(sys.stdout)
    # This handler needs to also have a level to deal with
    # mosaic filtering, see below
    stdout_handler.setLevel(log_level.upper())
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)
    if log_file is not None:
        # Set up a handler to route to file
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level.upper())
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    # It is useful to see logs from mosaic, but they are a bit of a mess
    # so here we add a filter to catch these and convert INFO to DEBUG
    mosaic_logger = logging.getLogger("mosaic")
    mosaic_logger.setLevel(log_level.upper())

    def mosaic_filter(record):
        if record.levelno == logging.INFO:
            record.levelno = logging.DEBUG
            record.levelname = logging.getLevelName(logging.DEBUG)
        return True

    # For stupid reasons filters don't propagate to child loggers
    # so we just do it manually here to avoid the issues that come
    # with setting the filter on the handler
    # pylint: disable=E1103
    for key in logging.root.manager.loggerDict.keys():
    # pylint: enable=E1103
        if key.startswith("mosaic"):
            logging.getLogger(key).addFilter(mosaic_filter)
    mosaic_logger.addHandler(stdout_handler)
    if log_file is not None:
        mosaic_logger.addHandler(file_handler)

def metadata_show(metafile: str, verbose: bool = False) -> None:
    """Display the contents of a metadata file.

    Args:
        verbose (bool, optional): Show a full breakdown of the file. Defaults to a summary.
    """
    om = skyweaver.SessionMetadata.from_file(metafile)
    if verbose:
        pprint.pprint(om)
    else:
        print(om)

def delays_create(
    metafile: str,
    bfconfig: str,
    pointing_idx: int = None,
    step: float = 4.0,
    outfile: str = None):
    """Create a delay file

    Args:
        metafile (str): Path to the session metadata file
        bfconfig (str): Path to the beamformer config file
        pointing_idx (int, optional): The pointing to generate delays for. Defaults to None.
        step (float, optional): The time step between delay solutions. Defaults to 4.0 seconds.
        outfile (str, optional): The file to write delay models to. 
                                 Defaults to a standard output filename.
    """
    sm = skyweaver.SessionMetadata.from_file(metafile)
    bc = skyweaver.BeamformerConfig.from_file(bfconfig)
    pointings = sm.get_pointings()
    npointings = len(pointings)
    if pointing_idx is None and npointings > 1:
        raise ValueError("Multiple pointings in session but no --pointing-idx specified")
    elif pointing_idx is None:
        pointing_idx = 0
    if pointing_idx >= npointings:
        raise ValueError("Pointing idx {} requested but only {} pointings in session")
    step = step * u.s
    pointing = pointings[pointing_idx]
    delays, _, _ = skyweaver.create_delays(sm, bc, pointing, step=step)
    if outfile is None:
        fname = "swdelays_{}_{}_to_{}_{}.bin".format(
            pointing.phase_centre.name,
            int(pointing.start_epoch.unix),
            int(pointing.start_epoch.unix),
            secrets.token_hex(3)
        )
    else:
        fname = outfile
    log.info("Writing delay model to file %s", fname)
    with open(fname, "wb") as fo:
        for delay_model in delays:
            fo.write(delay_model.to_bytes())

def parse_default_args(args):
    """Execute functions for common arguments
    """
    init_logger(args.log_level, args.log_file)
    
def add_defaults_args(parser):
    """Inject common arguments onto a parser
    """
    parser.add_argument("--log-level", help="Set the log level", type=str)
    parser.add_argument("--log-file", help="Specify and output logging file", type=str)

def subparser_create_wrapper(parent, *args, **kwargs):
    """A wrapper to add default args on subparser creation
    """
    parser = parent.add_parser(*args, **kwargs)
    add_defaults_args(parser)
    return parser

def cli():
    """Parse the skyweaver command line
    """
    parser = argparse.ArgumentParser(prog="skyweaver", add_help=True)
    l1subparsers = parser.add_subparsers(help="sub-command help")

    # sw metadata
    metadata = l1subparsers.add_parser("metadata", help="Tools for observation metadata files")
    metadata_subparsers = metadata.add_subparsers(help="sub-command help")

    # sw metadata show
    metadata_show_parser = subparser_create_wrapper(
        metadata_subparsers, "show", help="Show metadata")
    metadata_show_parser.add_argument("--verbose", action="store_true")
    metadata_show_parser.add_argument("metafile", metavar="FILE")
    metadata_show_parser.set_defaults(
        func=lambda args: metadata_show(args.metafile, args.verbose))

    # sw delays
    delays = l1subparsers.add_parser("delays", help="Tools for delay files")
    delays_subparsers = delays.add_subparsers(help="sub-command help")

    # sw delays create
    delays_create_parser = subparser_create_wrapper(
        delays_subparsers, "create", help="Create a delay file")
    delays_create_parser.add_argument("--step", dest="step", default=4, type=float,
                                    help="Step size in seconds between delay solutions")
    delays_create_parser.add_argument("--pointing-idx", dest="pointing_idx", default=None, type=int,
                                    help="The sub observation to create delays for")
    delays_create_parser.add_argument(
        "--outfile", dest="outfile", default=None, type=str,
        help="Output file to write delays to. Defaults to auto generated filename.")
    delays_create_parser.add_argument("metafile", metavar="METAFILE", 
                                      help="An HDF5 FBFUSE-BVR metadata file")
    delays_create_parser.add_argument("bfconfig", metavar="BFCONFIG", 
                                      help="A YAML beamformer configuration file")
    delays_create_parser.set_defaults(
        func=lambda args: delays_create(
            args.metafile,
            args.bfconfig,
            args.pointing_idx,
            args.step,
            args.outfile))

    # parse and execute
    args = parser.parse_args()
    parse_default_args(args)
    args.func(args)

if __name__ == "__main__":
    cli()
    
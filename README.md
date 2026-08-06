# skyweaver
Implementation of an offline FBFUSE beamformer for the MeerKAT telescope

# Installation

It is easiest to use the software inside a Docker container.  Two dockerfiles is included part of this repository: To compile the c++ and python parts respectively. 

# Usage

## Step 1:  Get delays for the beamformer

### Start with the skyweaver python CLI 

```bash
alias sw="python /path/to/python/skyweaver/cli.py"
sw --help
```
This will print:
```console
usage: skyweaver [-h] {metadata,delays} ...

positional arguments:
  {metadata,delays}  sub-command help
    metadata         Tools for observation metadata files
    delays           Tools for delay files

optional arguments:
  -h, --help         show this help message and exit
```

### Get the metadata for the corresponding observation 

This is done outside this repository

### Obtain metadata information

```bash
sw metadata show <metadata_file>
```

This will produce an output like the following:
```console
sw metadata show bvrmetadata_2024-02-16T10\:50\:46_72275.hdf5
---------------Array configuration----------------
Nantennas:          57
Subarray:           m000,m002,m003,m004,m005,m007,m008,m009,m010,m011,m012,m014,m015,m016,
                    m017,m018,m019,m020,m021,m022,m023,m024,m025,m026,m027,m029,m030,m031,
                    m032,m033,m034,m035,m036,m037,m038,m039,m040,m041,m042,m043,m044,m045,
                    m046,m048,m049,m050,m051,m053,m054,m056,m057,m058,m059,m060,m061,m062,
                    m063
Centre frequency:   1284000000.0 Hz
Bandwidth:          856000000.0 Hz
Nchannels:          4096
Sync epoch (UNIX):  1708039531.0
Project ID:         -
Schedule block ID:  -
CBF version:        cbf_dev
--------------------Pointings---------------------
#0   J1644-4559          2024-02-16T11:16:08.957000000 until 2024-02-16T11:21:04.865000000 (UTC)
                         1708082168.957 until 1708082464.865 (UNIX)
                         72996182384029 until 73502776880016 (SAMPLE CLOCK)
#1   J1644-4559_Offset1  2024-02-16T11:21:22.092000000 until 2024-02-16T11:26:15.682000000 (UTC)
                         1708082482.092 until 1708082775.682 (UNIX)
                         73532269504013 until 74034895583866 (SAMPLE CLOCK)
#2   J1644-4559_Offset2  2024-02-16T11:26:34.536000000 until 2024-02-16T11:31:25.723000000 (UTC)
                         1708082794.536 until 1708083085.723 (UNIX)
                         74067173632022 until 74565685776084 (SAMPLE CLOCK)
#3   M28                 2024-02-16T11:31:52.503000000 until 2024-02-16T12:01:51.651000000 (UTC)
                         1708083112.503 until 1708084911.651 (UNIX)
                         74611533136035 until 77691674512039 (SAMPLE CLOCK)
#4   J0437-4715          2024-02-16T12:03:06.117000000 until 2024-02-16T12:08:04.364000000 (UTC)
                         1708084986.117 until 1708085284.364 (UNIX)
                         77819160304176 until 78329759168140 (SAMPLE CLOCK)
```

### Create a config file in .yml format
Here is an example - there are comments to explain each parameter. 

```.yml 
created_by: Vivek
beamformer_config:
  # The total number of beams to be produced (must be a multiple of 32). This needs to be <= the number that SKYWEAVER is compiled for.
  total_nbeams: 800
  # The number of time samples that will be accumulated after detection, inside the beamformer
  tscrunch: 4
  # The number of frequency channels that will be accumulated after detection, inside the beamformer
  # Will be coerced to 1 if coherent dedispersion is specified.
  fscrunch: 1
  # The Stokes product to be calculated in the beamformer (I=0, Q=1, U=2, V=3)
  stokes_mode: 0
  # Enable CB-IB subtraction in the beamformer
  subtract_ib: True

  # Dispersion measure for coherent / incoherent dedispersion in pc cm^-3
  # A dispersion plan definition string "
  #           "(<coherent_dm>:<start_incoherent_dm>:<end_incoherent_dm>:<dm_step>:<tscrunch>) or "
  #           "(<coherent_dm>:<tscrunch>) "
  #           "or (<coherent_dm>)")
# Each DD plan is a "Stream" with zero indexed stream-ids 

ddplan:
  - "478.6:478.6:478.6:1:1" #stream-id=0
  - "0.00:478.6:478.6:1:1" #stream-id=1

# every beamset can contain arbitrary set of antennas, corresponding targeted beams, and tiled beams
# total number of beams across all beamsets should be <= the number of beams that SKYWEAVER is compiled for.
beam_sets:

  - antenna_set: ['m000','m002','m003','m004','m005','m007','m008','m009','m010','m011',
                  'm012','m014','m015','m016','m017','m018','m019','m020','m021','m022',
                  'm023','m024','m025','m026','m027','m029','m030','m031','m032','m033',
                  'm034','m035','m036','m037','m038','m039','m040','m041','m042','m043',
                  'm044','m045','m046','m048','m049','m050','m051','m053','m054','m056',
                  'm057','m058','m059','m060','m061','m062','m063']
    beams: []
    tilings:
      - nbeams: 32
        reference_frequency: null
        target: "J1644-4559,radec,16:44:49.273,-45:59:09.71"
        overlap: 0.9
```



### Create delay file for the corresponding pointing

```bash
sw delays create --pointing-idx 0 --outfile J1644-4559_pointing_0.delays --step 4 bvrmetadata_2024-02-16T10\:50\:46_72275.hdf5 J1644-4559_boresight.yaml
```

This produces a `.delays` file used for beamforming, and a `.targets` file that contains beam metadata. There are also other files produced here for reproducibility and for visualisation. 

## Step 2: Initialise input and compile skyweaver

### Create a list of dada files that correspond to the pointing

```console
ls /b/u/vivek/00_DADA_FILES/J1644-4559/2024-02-16-11\:16\:08/L/48/*dada -1 > /bscratch/vivek/skyweaver_tests/J1644-4559_boresight_dadafiles.list
```

### Compile skyweaver

This is done inside the dockerfile too. Either edit that to produce a docker image that has the software precompiled, or compile separately. 

```bash
    cmake -S . -B $cmake_tmp_dir -DENABLE_TESTING=0 -DCMAKE_INSTALL_PREFIX=$install_dir -DARCH=native -DPSRDADA_INCLUDE_DIR=/usr/local/include/psrdada -DPSRDADACPP_INCLUDE_DIR=/usr/local/include/psrdada_cpp -DSKYWEAVER_NANTENNAS=64 -DSKYWEAVER_NBEAMS=${nbeams} -DSKYWEAVER_NCHANS=64 -DSKYWEAVER_IB_SUBTRACTION=1 -DCMAKE_BUILD_TYPE=RELEASE -DSKYWEAVER_CB_TSCRUNCH=${tscrunch} -DSKYWEAVER_IB_TSCRUNCH=${tscrunch};
    cd $cmake_tmp_dir
    make -j 16
```
This compilation produces two binaries: `skyweavercpp` and `skycleaver`
## Step 3: Run the beamformer

```bash

/path/to/skyweavercpp --input-file J1644-4559_boresight_dadafiles.list --delay-file J1644-4559_pointing_0.delays --output-dir=/bscratch/vivek/skyweaver_out --gulp-size=32768 --log-level=warning --output-level=12 --stokes-mode I
```

Change the output level to 7 for bright pulsars like J1644-4559. 

This will produce `.tdb` files for the corresponding bridge. Run Step 3 for ALL 64 bridges with their corresponding dada file lists. These are DADA format files with the dimensions of TIME, INCOHERENT DM and BEAM as the order. For stokes I mode, The datatype is `int8_t`. For IQUV it is `char4`. 

## Steo 4: Cleave all bridges to form Filterbanks

Here we cleave the 64 TDB[I/Q/U/V/IV/QU/IQUV] files to produce `NDM*NBEAMS*NSTOKES` number of T(F=64) files. 

to run this, do

```bash

/path/to/skycleaver -r /bscratch/vivek/skyweaver_out --output-dir /bscratch/vivek/skycleaver_out --nsamples-per-block 65536 --nthreads 32 --stream-id 0 --targets_file/bscratch/vivek/skyweaver_out/swdlays_J1644-4559.targets --out-stokes I --required_beams 0 

```

This will produce a standard sigproc format `.fil` file that can be used for traditional processing. 




# ncspiR

An R script to compute the SPI index on netCDF files.

Example usage from a Unix terminal, using the example [input file](input/example_input_std_short_small.nc) included in this repo:
```
> ./ncspiR.R --progress input/example_input_std_short_small.nc output/example_output_std_short_small.nc
INFO [2019-05-16 10:43:10]  ### Starting ncspiR.R version 0.2.1 from Adriano Fantini (afantini@ictp.it) ###
INFO [2019-05-16 10:43:10] Loading packages
INFO [2019-05-16 10:43:10] Reading metadata
INFO [2019-05-16 10:43:11] Reading data
INFO [2019-05-16 10:43:13] Starting computation using 12 threads
   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed = 01m 49s
INFO [2019-05-16 10:45:02] Ended computation
INFO [2019-05-16 10:45:03] Initialising output file
INFO [2019-05-16 10:45:04] Filling output file
INFO [2019-05-16 10:45:04]  ### Goodbye! ###
```
The output should be identical to [this file](output/example_output_std_short_small.nc).

Depending on the files you want to process and on your particular system, running with less cores than the default (all cores available) or even in serial (`--n 1` option) might speed up computation.

Many options are available, use `--help` to list them:
```
> ./ncspiR.R --help
Usage: ./ncspiR.R [options] INPUT OUTPUT


Options:
        -t TIMESCALE, --timescale=TIMESCALE
                SPI timescale in months [default: 12]

        -n NTHREADS, --nthreads=NTHREADS
                Number of threads to be used. Set to NA for automatic detection (all available threads in the node will be used) [default: NA]

        -v VARIN, --varin=VARIN
                Variable name to use from input file [default: pr]

        -o VAROUT, --varout=VAROUT
                Variable name to use in output file [default: SPI]

        --nafraction=NAFRACTION
                Max fraction of accepted NAs in a timeseries (0.05 = 5%). If more than this, that cell point becomes NA. [default: 0.05]

        --maxspi=MAXSPI
                SPI output values greater than this (in absolute value) will become NA. [default: 1e+10]

        --refstart=REFSTART
                ref.start parameter to pass to SPEI::spi, a character in the YEAR-MON format (e.g. 1976-01). [default: NULL]

        --refend=REFEND
                ref.end parameter to pass to SPEI::spi, a character in the YEAR-MON format (e.g. 2005-12). [default: NULL]

        -l LOGFILE, --logfile=LOGFILE
                Optional file to write logs to. [default: NULL]

        --progress
                Show a progress bar - this slightly decreases performance

        --assume_monthly
                Assume the input file has the correct monthly periodicity, and only use the time of the first timestep to define times

        --compress
                Activate netCDF compression (with deflate level 1) for the SPI variable

        --debug
                Print additional debug output. This flag is also useful if you want to check that the options were correctly understood

        --dryrun
                Perform a dry run, do not compute nor write anything to file

        -h, --help
                Show this help message and exit


#================= DESCRIPTION =================#
ncspiR.R version 0.2.1 from Adriano Fantini (afantini@ictp.it)

Script to calculate the SPI index from a monthly netCDF file containing precipitation data.
Input file MUST be monthly. Does everything in memory, so make sure your dataset fits in memory!
 Very few checks are performed, so also make sure you know what you are doing.
 Input files must follow the CF Conventions >= 1.5 (http://cfconventions.org/).
 This program is parallel by default, but is not capable of crossing node boundaries (cannot currently run on multiple nodes).

REQUIRED PACKAGES: pacman, optparse, glue, SPEI, magrittr, futile.logger, lubridate, ncdf4, PCICt, pbapply
MISSING  PACKAGES: none, well done!

Note that ncspiR.R will *attempt* to install any missing packages by itself.

#================= GET IN TOUCH ================#
For feature requests and bug reports, please get in touch on github:
https://github.com/adrfantini/ncspiR/issues
```

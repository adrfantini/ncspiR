#!/usr/bin/env Rscript

program_name = 'ncspiR.R'
description = 'Script to calculate the SPI index from a monthly netCDF file containing precipitation data.'
further_description = 'Input file MUST be monthly. Does everything in memory, so make sure your dataset fits in memory! \nVery few checks are performed, so also make sure you know what you are doing.\n Input files must follow the CF Conventions >= 1.5 (http://cfconventions.org/).'
author = 'Adriano Fantini'
version = '0.1'
contact = 'afantini@ictp.it'
gh_url = 'https://github.com/adrfantini/ncspiR'
required_pkgs = c(
    'pacman',
    'optparse',
    'glue',
    'stars',
    'SPEI',
    'magrittr',
    'futile.logger',
    'lubridate',
    'raster'
)

#============= INITIALIZATION =============

fatalerror = function() stop('fatal error', call.=FALSE)
if (!'pacman' %in% rownames(installed.packages())) {
    cat('Please first install package pacman from the R command line with:
        install.packages(pacman)
    pacman will take care of the installation of other necessary packages
    ')
    fatalerror()
}
suppressPackageStartupMessages(library(pacman))
options("pac_update" = FALSE) # Do not try to update relevant packages on package install with p_load
missing = required_pkgs[!p_isinstalled(required_pkgs)]
if (length(missing) > 0L) {
    missing = paste0(missing, collapse = ', ')
} else {
    missing = 'none, well done!'
}
suppressPackageStartupMessages(p_load(optparse))
suppressPackageStartupMessages(p_load(glue))
suppressPackageStartupMessages(p_load(futile.logger))
flog.fatal = function(...) {futile.logger::flog.fatal(...); fatalerror()}

option_list = list(make_option(c("-t", "--timescale"),
                                type="integer",
                                default=12,
                                help="SPI timescale in months [default: %default]"),
                    make_option(c("-n", "--nthreads"),
                                type="integer",
                                default=NA,
                                help="Number of threads to be used. Set to NA for automatic detection (all available threads will be used) [default: %default]"),
                    make_option(c("-v", "--varin"),
                                type="character",
                                default='pr',
                                help="Variable name to use from input file [default: %default]"),
                    make_option(c("-o", "--varout"),
                                type="character",
                                default='SPI',
                                help="Variable name to use in output file [default: %default]"),
                    make_option("--nathreshold",
                                type="double",
                                default=5/100,
                                help="Max fraction of accepted NAs in a timeseries (0.05 = 5%). If more than this, that cell point becomes NA. [default: %default]"),
                    make_option("--maxspi",
                                type="double",
                                default=1e10,
                                help="SPI output values greater than this (in absolute value) will become NA. [default: %default]"),
                    make_option("--refstart",
                                type="character",
                                default=NULL,
                                help="ref.start parameter to pass to SPEI::spi, a character in the YEAR-MON format (e.g. 1976-01). [default: %default]"),
                    make_option("--refend",
                                type="character",
                                default=NULL,
                                help="ref.end parameter to pass to SPEI::spi, a character in the YEAR-MON format (e.g. 2005-12). [default: %default]"),
                    make_option(c('-l', "--logfile"),
                                type="character",
                                default=NULL,
                                help="Optional file to write logs to. [default: %default]"),
                    make_option("--debug",
                                action="store_true",
                                help="Print additional debug output. This flag is also useful if you want to check that the options were correctly understood [default: %default]")
                    )
parser = OptionParser(
    usage = "%prog [options] INPUT OUTPUT",
    option_list=option_list,
    epilogue=glue("

    #================= DESCRIPTION =================#
    {program_name} version {version} from {author} ({contact})

    {description}
    {further_description}

    REQUIRED PACKAGES: {paste(required_pkgs, collapse=', ')}
    MISSING: {missing}

    Note that {program_name} will *attempt* to install any missing packages by itself.

    #================= GET IN TOUCH ================#
    For feature requests and bug reports, please get in touch on github:
    {gh_url}/issues

    ")
)

#============= INPUT ARGUMENTS =============

arguments = parse_args(parser, positional_arguments = 2) # , args = c('input/example_input_std_short.nc', 'output/example_output_std_short.nc'))
opt = arguments$options

debug = isTRUE(opt$debug)
if (debug) invisible(flog.threshold(DEBUG))

logfile = opt$logfile
if (!is.null(logfile)) {
    if (!file.exists(logfile)) ftry(file.create(logfile))
    if (file.access(logfile, mode=2) == -1) {
        flog.fatal('cannot write to log file, check permissions and path')
    }
    flog.appender(appender.tee(logfile))
}

flog.debug("Parsing input arguments")

ts = opt$timescale
if (!is.integer(ts)) flog.fatal('timescale must be integer, got "%s".', ts)
if (ts < 1) flog.fatal('timescale must be greater than 1, got %d.', ts)
flog.debug("Timescale: %d months", ts)

fn_in = arguments$args[1]
flog.debug("Input file: %s", fn_in)
if (!file.exists(fn_in)) flog.fatal('input file does not exist')
if (file.access(fn_in, mode=4) == -1) flog.fatal('input file cannot be read, check permissions and path')

fn_out = arguments$args[2]
flog.debug("Output file: %s", fn_out)
if (file.exists(fn_out)) flog.fatal('input file already exists')
if (file.access(dirname(fn_in), mode=2) == -1) flog.fatal('output file cannot be created, check permissions and path')

var_in = opt$varin
flog.debug("Input variable name: %s", var_in)
var_out = opt$varout
flog.debug("Output variable name: %s", var_out)

nthreads = opt$nthreads
if (is.na(nthreads)) nthreads = parallel::detectCores()
if (!is.integer(nthreads)) flog.fatal('Number of threads must be integer, got "%s".', nthreads)
if (nthreads < 1) flog.fatal('Number of threads must be greater than 1, got %d.', nthreads)
flog.debug("Number of threads: %d", nthreads)

na_thr = opt$nathreshold
if (!is.numeric(na_thr)) flog.fatal('NA threshold must be integer, got "%s".', na_thr)
if (na_thr < 0) flog.fatal('NA threshold must be greater than 0, got %g.', na_thr)
if (na_thr >= 1) {
    flog.warn('NA threshold must be smaller than 1, got %g; rescaling by a factor 100', na_thr)
    na_thr = na_thr/100
}
flog.debug("NA threshold: %g (%g%%)", na_thr, na_thr*100)

spi_na_thr = opt$maxspi
if (!is.numeric(na_thr)) flog.fatal('SPI max threshold must be integer, got "%s".', spi_na_thr)
flog.debug("SPI max threshold: %g", spi_na_thr)

ref_start = opt$refstart
ref_end = opt$refend
notmonth = function(x) {
    x>12 | x<1
}
if (!is.null(ref_start)) {
    ref_start = as.numeric(strsplit(ref_start, '-')[[1]])
    if (length(ref_start) != 2L || notmonth(ref_start[2])) flog.fatal('Did not understand the reference start value, got %s', paste(ref_start, collapse='-'))
    flog.debug("Reference start: %d-%d", ref_start[1], ref_start[2])
}
if (!is.null(ref_end  )) {
    ref_end   = as.numeric(strsplit(ref_end  , '-')[[1]])
    if (length(ref_end) != 2L || notmonth(ref_end[2])) flog.fatal('Did not understand the reference end value, got %s', paste(ref_end, collapse='-'))
    flog.debug("Reference end: %d-%d", ref_end[1], ref_end[2])
}

#============= PACKAGE LOAD =============

flog.info("Loading packages")
suppressPackageStartupMessages(p_load(required_pkgs, character.only = TRUE))
flog.debug("Loaded packages")

#============= READ INPUT DATA =============

flog.info("Reading metadata")
nc_in = suppressWarnings(fn_in %>% read_stars(sub = var_in, proxy = TRUE))

if (!st_dimensions(nc_in)$time$refsys %in% c('PCICt', 'POSIXct', 'POSIXlt')) {
    flog.fatal('Input calendar not understood. Are you sure your input file follows the CF Conventions?')
}
times = suppressWarnings(st_get_dimension_values(nc_in, 'time') %>% as.POSIXct %>% round('day'))
day(times) = 15
times = as.Date(times)
#check periodicity
expected_times = times[1] + months(1:length(times)-1)
if ( !isTRUE(all( times == expected_times) ) ) flog.fatal('Incorrect monthly periodicity in the input file. Missing any timesptes?')

start_ym = c(year(times[1]), month(times[1]))
end_ym = c(year(times[length(times)]), month(times[length(times)]))
flog.debug("First month in the file: %s", paste(start_ym, collapse='-'))
flog.debug("Last  month in the file: %s", paste(end_ym, collapse='-'))

flog.info("Reading data")
nc_in = suppressWarnings(fn_in %>% read_stars(sub = var_in, proxy = FALSE))

#============= COMPUTE =============

calc_SPI = function(d, ts, thr, ref.s=NULL, ref.e=NULL, first_ym, last_ym) {
    d <- as.numeric(d)
    if (sum(is.na(d))/length(d) > thr) {
        rep(NA, length(d))
    } else {
        if (is.null(ref.s)) {
            ref.s = first_ym
        }
        if (is.null(ref.e)) {
            ref.e = last_ym
        }
        d = ts(d, frequency=12, start=first_ym)
        as.vector(fitted(SPEI::spi(d, ts, na.rm=TRUE, ref.start=ref.s, ref.end=ref.e)))
    }
}
if (nthreads > 1 ) {
    p_load(parallel)
    cluster = makeCluster(nthreads)
} else {
    cluster = NULL
}

flog.info("Starting computation using %d threads", nthreads)
spi_res <- st_apply(nc_in,
    which((st_dimensions(nc_in) %>% names) != 'time'),
    calc_SPI, CLUSTER = cluster, PROGRESS = TRUE,
    ts = ts, thr = na_thr, ref.s = ref_start, ref.e = ref_end, first_ym = start_ym, last_ym = end_ym
)
if (nthreads > 1 ) stopCluster(cluster)
names(spi_res) = var_out
flog.info("Ended computation")

if (debug) {
    n_over_thr = sum(abs(spi_res[[var_out]]) > spi_na_thr, na.rm=TRUE)
    flog.debug('Setting %d values above the SPI threshold (%g) to NA', n_over_thr, spi_na_thr)
}
spi_res[abs(spi_res) > spi_na_thr] = NA

#============= WRITE OUTPUT =============

flog.info('Writing to file %s', fn_out)
spell_res %>%
    as('Raster') %>%
    setZ(times) %>%
    writeRaster(fn_out, varname = var_out, varunit = '1', longname = 'SPI index', zname = 'time')

flog.info('All done')

# rnc_in = open.nc(fn_in)
# rnc_out = create.nc(fn_out)
# modifyNcdfCopyMetadata(rnc_in, rnc_out)

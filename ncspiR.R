#!/usr/bin/env Rscript

program_name = 'ncspiR.R'
description = 'Script to calculate the SPI index from a monthly netCDF file containing precipitation data.'
further_description = 'Input file MUST be monthly. Does everything in memory, so make sure your dataset fits in memory!\n Very few checks are performed, so also make sure you know what you are doing.\n Input files must follow the CF Conventions >= 1.5 (http://cfconventions.org/).\n This program is parallel by default, but is not capable of crossing node boundaries (cannot currently run on multiple nodes).'
author = 'Adriano Fantini'
version = '0.2'
contact = 'afantini@ictp.it'
gh_url = 'https://github.com/adrfantini/ncspiR'
required_pkgs = c(
    'pacman',
    'optparse',
    'glue',
    'SPEI',
    'magrittr',
    'futile.logger',
    'lubridate',
    'ncdf4',
    'PCICt',
    'pbapply'
)
version_pkgs = c( # Set requirements for package versions
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
# options(error = function() { flog.fatal(geterrmessage()) ; quit(runLast=FALSE)}) # Override R's default error handling

option_list = list(make_option(c("-t", "--timescale"),
                                type="integer",
                                default=12,
                                help="SPI timescale in months [default: %default]"),
                    make_option(c("-n", "--nthreads"),
                                type="integer",
                                default=NA,
                                help="Number of threads to be used. Set to NA for automatic detection (all available threads in the node will be used) [default: %default]"),
                    make_option(c("-v", "--varin"),
                                type="character",
                                default='pr',
                                help="Variable name to use from input file [default: %default]"),
                    make_option(c("-o", "--varout"),
                                type="character",
                                default='SPI',
                                help="Variable name to use in output file [default: %default]"),
                    make_option("--nafraction",
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
                    make_option("--progress",
                                action="store_true",
                                help="Show a progress bar - this slightly decreases performance"),
                    make_option("--compress",
                                action="store_true",
                                help="Activate netCDF compression (with deflate level 1) for the SPI variable"),
                    make_option("--debug",
                                action="store_true",
                                help="Print additional debug output. This flag is also useful if you want to check that the options were correctly understood")
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
        flog.fatal('Cannot write to log file, check permissions and path')
    }
    invisible(flog.appender(appender.tee(logfile)))
    if (file.exists(logfile)) flog.debug('Appending output to logfile %s', logfile)
}

flog.info(glue(' ### Starting {program_name} version {version} from {author} ({contact}) ### '))

progress = isTRUE(opt$progress)
if (progress) {
    flog.debug('Progress bar will be shown')
} else {
    flog.debug('Progress bar will not be shown')
}

compress = isTRUE(opt$compress)
deflate_level = 1
if (compress) flog.debug('Deflate compression activated')

flog.debug("Parsing input arguments")

ts = opt$timescale
if (!is.integer(ts)) flog.fatal('timescale must be integer, got "%s".', ts)
if (ts < 1) flog.fatal('timescale must be greater than 1, got %d.', ts)
flog.debug("Timescale: %d months", ts)

fn_in = arguments$args[1]
flog.debug("Input file: %s", fn_in)
if (!file.exists(fn_in)) flog.fatal('Input file does not exist')
if (file.access(fn_in, mode=4) == -1) flog.fatal('Input file cannot be read, check permissions and path')

fn_out = arguments$args[2]
flog.debug("Output file: %s", fn_out)
if (file.exists(fn_out)) flog.fatal('Output file already exists')
if (file.access(dirname(fn_out), mode=2) == -1) flog.fatal('Output file cannot be created, check permissions and path')

var_in = opt$varin
flog.debug("Input variable name: %s", var_in)
var_out = opt$varout
flog.debug("Output variable name: %s", var_out)

nthreads = opt$nthreads
if (is.na(nthreads)) nthreads = parallel::detectCores()
if (!is.integer(nthreads)) flog.fatal('Number of threads must be integer, got "%s".', nthreads)
if (nthreads < 1) flog.fatal('Number of threads must be greater than 1, got %d.', nthreads)
flog.debug("Number of threads: %d", nthreads)

na_thr = opt$nafraction
if (!is.numeric(na_thr)) flog.fatal('NA threshold must be integer, got "%s".', na_thr)
if (na_thr < 0) flog.fatal('NA threshold must be greater than 0, got %g.', na_thr)
if (na_thr >= 1) {
    flog.warn('NA threshold must be smaller than 1, got %g; rescaling by a factor 100', na_thr)
    na_thr = na_thr/100
}
flog.debug("NA threshold: %g (%g%%)", na_thr, na_thr*100)

spi_max_thr = opt$maxspi
if (!is.numeric(spi_max_thr)) flog.fatal('SPI max threshold must be integer, got "%s".', spi_max_thr)
flog.debug("SPI max threshold: %g", spi_max_thr)

ref_start = opt$refstart
ref_end = opt$refend
notmonth = function(x) {
    x>12 | x<1
}
if (is.null(ref_end) != is.null(ref_start)) {
    flog.fatal('--refstart and --refend must either be both present, or both absent')
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

    # Check ref_end is successive to ref_start
    ref_limits_invalid = (ref_start[1] > ref_end[1]) | (ref_start[1] == ref_end[1] & ref_start[2] >= ref_end[2])
    if (ref_limits_invalid) flog.fatal('--refend (%s) must be successive to --refstart (%s)', paste(ref_end, collapse='-'), paste(ref_start, collapse='-'))
}

#============= PACKAGE LOAD =============

flog.info("Loading packages")
suppressPackageStartupMessages(p_load(required_pkgs, character.only = TRUE))
version_pkgs %>% names %>% sapply(function(x) {
    ver = p_ver(x)
    if ( ver < version_pkgs[x] ) flog.fatal('Version %s of package %s required, got %s. Maybe your (probably CRAN) version is behind the github one?', version_pkgs[x], x, as.character(ver))
}) %>% invisible
flog.debug("Loaded packages")

#============= READ INPUT DATA =============

flog.info("Reading metadata")
nc_in = nc_open(fn_in)

# Read time
time_var = 'time'
times = nc_in %>% ncvar_get(time_var)
time_units = nc_in %>% ncatt_get(time_var, 'units')
if (!time_units$hasatt) flog.fatal('Cannot find time units!')
flog.debug('Time units: %s', time_units$value)
time_units = strsplit(time_units$value, " ")[[1]]
time_cal = nc_in %>% ncatt_get(time_var, 'calendar')
if (!time_cal$hasatt) {
    flog.info('No input calendar, assuming standard')
    time_cal = 'gregorian'
} else {
    time_cal = time_cal$value %>% tolower
    flog.debug('Input calendar: %s', time_cal)
}

pcict_calendars = c('365', '365_day', 'noleap', '360', '360_day', 'gregorian', 'proleptic_gregorian', 'standard')
if (time_cal %in% pcict_calendars) {
    # Use PCIct to deal with times
    time_offset_unit <- time_units[1]
    time_tz <- time_units[5]
    time_start <- strsplit(time_units[4], ":")[[1]]
    if (length(time_start) != 3 || time_start[1] > 24 || time_start[2] > 60 || time_start[3] > 60 || any(time_start < 0)) {
        flog.warn("%s is not a valid start time. Assuming 00:00:00", time_start)
        time_units[4] <- "00:00:00"
    }
    if (! time_tz %in% OlsonNames()) {
        flog.debug("%s is not a valid timezone. Assuming UTC", time_tz)
        time_tz <- "UTC"
    }
    time_start <- ymd_hms(paste(time_units[3], time_units[4]), tz=time_tz)
    # Find the correct lubridate time function based on the unit
    time_f <- switch(tolower(time_offset_unit),
        seconds=seconds, second=seconds, sec=seconds,
        minutes=minutes, minute=minutes, min=minutes,
        hours=hours,     hour=hours,     h=hours,
        days=days,       day=days,       d=days,
        months=months,   month=months,   m=months,
        years=years,     year=years,     yr=years,
        NA
    )
    times = as.PCICt(time_start, cal=time_cal) + as.numeric(time_f(floor(times))) # TODO this floor here is a hack, to work around files which have non-integer times. Works in most cases.
} else {
    flog.fatal('This program does not support calendar ""%s"', time_cal)
}

# Simplify dates, since we only care about months
times = suppressWarnings(times %>% as.POSIXct %>% round('day'))
day(times) = 15
times = as.Date(times)

# Check times
expected_times = times[1] + months(1:length(times)-1)
if ( !isTRUE(all( times == expected_times) ) ) flog.fatal('Incorrect monthly periodicity in the input file. Missing any timesteps?')

# Check reference times are acceptable
nc_start = c(year(times[1]), month(times[1]))
nc_end = c(year(times[length(times)]), month(times[length(times)]))
flog.debug("First month in the file: %s", paste(nc_start, collapse='-'))
flog.debug("Last  month in the file: %s", paste(nc_end, collapse='-'))

# Check reference start and end (again)
if (!is.null(ref_start)) {
    ref_start_invalid = (nc_start[1] > ref_start[1]) | (nc_start[1] == ref_start[1] & nc_start[2] > ref_start[2])
    if (ref_start_invalid) flog.fatal('--refstart (%s) antecedent to the first timestep in the file (%s)', paste(ref_start, collapse='-'), paste(nc_start, collapse='-'))
}
if (!is.null(ref_end)) {
    ref_end_invalid = (nc_end[1] < ref_end[1]) | (nc_end[1] == ref_end[1] & nc_end[2] < ref_end[2])
    if (ref_end_invalid) flog.fatal('--refend (%s) successive to the last timestep in the file (%s)', paste(ref_end, collapse='-'), paste(nc_end, collapse='-'))
}

ncdim_get = function(nc, varid) {
    nc$var[[varid]]$dim %>% sapply(`[[`, 'name')
}

# Inspect input variable to read
nc_var_in = nc_in$var[[var_in]]
if (is.null(nc_var_in)) flog.fatal('Input file %s does not contain variable %s', fn_in, var_in)
nc_var_in_dims = nc_in %>% ncdim_get(var_in)

# Read data
flog.info("Reading data")
spi_in = nc_in %>% ncvar_get(var_in)

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
# For testing purposes
# calc_SPI = function(d, ts, thr, ref.s=NULL, ref.e=NULL, first_ym, last_ym) {
#     rep(0, length(d))
# }

if (nthreads > 1 ) {
    p_load(parallel)
    cluster = makeCluster(nthreads)
} else {
    cluster = NULL
}

flog.info("Starting computation using %d threads", nthreads)
if (progress) {
    pboptions('type' = 'timer')
} else {
    pboptions('type' = 'none')
}
spi_out <- pbapply(spi_in,
    which(nc_var_in_dims != 'time'),
    calc_SPI, cl = cluster,
    ts = ts, thr = na_thr, ref.s = ref_start, ref.e = ref_end, first_ym = nc_start, last_ym = nc_end
)
if (nthreads > 1 ) stopCluster(cluster)
flog.info("Ended computation")

if (debug) {
    n_over_thr = sum(abs(spi_out) > spi_max_thr, na.rm=TRUE)
    flog.debug('Setting %d values above the SPI threshold (%g) to NA', n_over_thr, spi_max_thr)
}
spi_out[abs(spi_out) > spi_max_thr] = NA

# For some reason it is necessary to permutate the output. Further checks are necessary to verify this is the correct permutation # TODO
shift_vec = function(v) {
    len = length(v)
    return(c(v[2:len], v[1]))
}
spi_out = aperm(spi_out, shift_vec(1:length(dim(spi_out))))

#============= WRITE OUTPUT =============

# Get dimensions and variables to be cloned from the source file
# coordinates and grid_mapping are special variables that act as pointers to other variables/dimensions
nc_coords = nc_in %>% ncatt_get(var_in, 'coordinates')
nc_grid_mapping = nc_in %>% ncatt_get(var_in, 'grid_mapping')
vars2copy = NULL
if (nc_coords$hasatt) {
    flog.debug('Found "coordinates" attribute for variable %s', var_in)
    vars2copy = nc_coords$value %>% strsplit(' ') %>% unlist
}
if (nc_grid_mapping$hasatt){
    flog.debug('Found "grid_mapping" attribute for variable %s', var_in)
    vars2copy = c(vars2copy, nc_grid_mapping$value %>% strsplit(' ') %>% unlist)
}
dims2copy = c(nc_var_in_dims, lapply(vars2copy, function(v) nc_in %>% ncdim_get(v)) %>% unlist) %>% unique
flog.debug("Variables  to be cloned: {paste(vars2copy, collapse=', ')}" %>% glue)
flog.debug("Dimensions to be cloned: {paste(dims2copy, collapse=', ')}" %>% glue)

# Creating output file
flog.info('Initialising output file')
nc_var_out = ncvar_def(
    var_out, '1', nc_in$dim[nc_var_in_dims],
    longname='SPI index',
    shuffle = compress,
    compression = ifelse(compress, deflate_level, NA)
)
nc_out = fn_out %>% nc_create(c(list(nc_var_out), nc_in$var[vars2copy]), force_v4 = TRUE)

# Wrapper for ncatt_put to log in debug
ncatt_put = function(nc, varid, attname, attval, ...) {
    flog.debug('Setting %s attribute: %s = %s', ifelse(varid == 0, 'global', varid), attname, attval)
    ncdf4::ncatt_put(nc, varid, attname, attval, ...)
}

# Clone all the attributes from the relevant variables and dimensions in the source file
for (v in c('global attributes', dims2copy, vars2copy)) {
    if (v == 'global attributes') v = 0
    if (v == 0) {
        flog.debug('Getting global attributes')
    } else {
        flog.debug('Getting attributes for variable %s', v)
    }

    atts = nc_in %>% ncatt_get(v)
    for (a in names(atts)) {
        nc_out %>% ncatt_put(v, a, atts[[a]])
    }
}

# Update netCDF history
nc_h = ncatt_get(nc_in, 0, 'history')
nc_new_h = '{date()}: SPI index calculated by {program_name} version {version} ({gh_url})' %>% glue
if (nc_h$hasatt) {
    nc_out %>% ncatt_put(0, 'history', paste(nc_new_h, nc_h$value, sep='\n'))
} else {
    nc_out %>% ncatt_put(0, 'history', nc_new_h)
}

# Fill in cloned variables
for (v in vars2copy) {
    flog.debug('Copying values for cloned variable %s', v)
    nc_out %>% ncvar_put(v, nc_in %>% ncvar_get(v))
}
flog.debug('Closing input file')
nc_close(nc_in)

# Clone coordinates and grid_mapping attributes, if present
if (nc_coords$hasatt) nc_out %>% ncatt_put(var_out, 'coordinates', nc_coords$value)
if (nc_grid_mapping$hasatt) nc_out %>% ncatt_put(var_out, 'grid_mapping', nc_grid_mapping$value)

# Fill SPI variable with data and attributes
flog.info('Filling output file')
nc_out %>% ncvar_put(var_out, spi_out)
nc_out %>% ncatt_put(var_out, 'timescale', paste(ts, 'months'))
if (!is.null(ref_start)) nc_out %>% ncatt_put(var_out, 'ref_start', paste(ref_start, collapse='-'))
if (!is.null(ref_end  )) nc_out %>% ncatt_put(var_out, 'ref_end'  , paste(ref_end, collapse='-'))
nc_out %>% ncatt_put(var_out, 'max_spi', spi_max_thr)
nc_out %>% ncatt_put(var_out, 'na_thr', na_thr)
nc_out %>% ncatt_put(var_out, 'program', '{program_name} version {version} ({gh_url})' %>% glue)
nc_out %>% ncatt_put(var_out, 'R_version', p_ver(R) %>% as.character)
nc_out %>% ncatt_put(var_out, 'SPEI_pkg_version', p_ver(SPEI) %>% as.character)
nc_out %>% ncatt_put(var_out, 'ncdf4_pkg_version', p_ver(ncdf4) %>% as.character)

# Close
flog.debug('Closing output file')
nc_close(nc_out)
flog.info(' ### Goodbye! ### ')

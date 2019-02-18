#!/usr/bin/env Rscript
# Feb 2019
# Script to calculate SPI, new and better version

#============= INITIALIZATION =============

suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(optparse))

option_list = list(make_option(c("-t", "--timescale"),
                                type="integer",
                                default=12,
                                help="Timescale in months [default: %default]"),
                    make_option(c("-n", "--nthreads"),
                                type="integer",
                                default=NA,
                                help="Number of threads to be used, NA for automatic detection [default: %default]"),
                    make_option(c("-i", "--varin"),
                                type="character",
                                default=NA,
                                help="Variable to use from input file, NA for automatic detection [default: %default]"),
                    make_option(c("-o", "--varout"),
                                type="character",
                                default='SPI',
                                help="Variable to use in output file [default: %default]"),
                    make_option("--nathreshold",
                                type="double",
                                default=5/100,
                                help="Max amount of accepted NAs in a timeseries. If more, that cell point becomes NA. [default: %default]"),
                    make_option("--spinathreshold",
                                type="double",
                                default=1e10,
                                help="SPI output values greater than this in absolute value will become NA. [default: %default]"),
                    make_option("--refstart",
                                type="character",
                                default=NULL,
                                help="ref.start parameter to pass to SPEI:spi, a character in the YEAR-MON format. [default: %default]"),
                    make_option("--refend",
                                type="character",
                                default=NULL,
                                help="ref.end parameter to pass to SPEI:spi, a character in the YEAR-MON format. [default: %default]"),
                    make_option("--useraster",
                                action="store_true",
                                help="Use package 'raster' to read the file, then converts to 'stars', instead of using 'stars' directly. [default: %default]")
                    )
parser = OptionParser(
    usage = "%prog [options] INPUT OUTPUT",
    option_list=option_list,
    epilogue=paste0("This script calculates the SPI index for a given variable in a NetCDF file, and outputs a similar file. Input file MUST be monthly. Does everything in memory, so make sure your dataset fits in memory! No checks are performed, make sure you know what you are doing. For feature requests and bug reports, please write to:
    afantini@ictp.it

    ###################################################################################
    WARNING! CURRENT VERSION INVERTS LATITUDE DATA! Can be fixed with cdo invertlatdata
    ###################################################################################
    ")
)

#Gather input arguments
arguments = parse_args(parser, positional_arguments = 2) # , args = c('filein.nc', 'fileout.nc')) , args = c('--refend', '2005-12', 'prOK_mon_1981-2050.nc', 'test.nc'))
opt = arguments$options

suppressPackageStartupMessages(p_load(futile.logger))

flog.debug("Parsed input arguments")
fn_in = arguments$args[1]
flog.info("Input file: %s", fn_in)
fn_out = arguments$args[2]
flog.info("Output file: %s", fn_out)
timescale = opt$timescale
flog.info("Timescale: %d", timescale)
nthreads = opt$nthreads
var_in = opt$varin
var_out = opt$varout
na_thr = opt$nathreshold
spi_na_thr = opt$spinathreshold
nthreads = opt$nthreads
ref_start = opt$refstart
ref_end = opt$refend
useraster = isTRUE(opt$useraster)
flog.info("NA threshold: %f", na_thr)

notmonth = function(x) {
    x>12 | x<1
}

if (!is.null(ref_start)) {
    ref_start = as.numeric(strsplit(ref_start, '-')[[1]])
    if (length(ref_start) != 2L || notmonth(ref_start[2])) stop(flog.fatal('Did not understand the --refstart value: %d-%d', ref_start[1], ref_start[2], return=T))
    flog.info("Reference start: %d-%d", ref_start[1], ref_start[2])
}
if (!is.null(ref_end  )) {
    ref_end   = as.numeric(strsplit(ref_end  , '-')[[1]])
    if (length(ref_end) != 2L || notmonth(ref_end[2])) stop(flog.fatal('Did not understand the --refend value: %d-%d', ref_end[1], ref_end[2], return=T))
    flog.info("Reference end: %d-%d", ref_end[1], ref_end[2])
}
if ((!is.null(ref_start) || !is.null(ref_end)) && useraster) {
    stop(flog.fatal("Package 'raster' does not use dates properly, thus the current implementation of '--useraster' is incompatible with the 'refend' and 'refstart' options."))
}

#============= PACKAGE LOAD =============

flog.info("Loading packages")
suppressPackageStartupMessages(p_load(stars))
suppressPackageStartupMessages(p_load(ncdf4))
suppressPackageStartupMessages(p_load(furrr))
suppressPackageStartupMessages(p_load(SPEI))
suppressPackageStartupMessages(p_load(lubridate))
suppressPackageStartupMessages(p_load(magrittr))


#============= COMPUTE =============

flog.info("Reading data")
if (useraster) {
    suppressPackageStartupMessages(p_load(raster))
    if (is.na(var_in)) {
        b = brick(fn_in)
        nc_in = st_as_stars(b)
        var_in = b@data@zvar
        rm(b)
        names(nc_in) = var_in
    } else {
        nc_in = st_as_stars(brick(fn_in, varname = var_in))
    }
} else {
    if (is.na(var_in)) {
        nc_in = read_stars(fn_in)[1,,,]
        source('/home/netapp-clima-users1/afantini/R/functions/nc_guess_var.R')
        var_in = nc_guess_var(fn_in) # terrible way to do it...
    } else {
        nc_in = read_stars(fn_in, sub = var_in)
    }
}
flog.info("Input variable: %s", names(nc_in))

times = st_dimensions(nc_in)$time$values
if (isTRUE(is.null(times))) {
    times = st_dimensions(nc_in)$time$offset %>% as.POSIXlt + with(st_dimensions(nc_in)$time, (from:to)*delta)
    flog.warn('Applying workaround for time dimension. Assuming monthly data, no checks performed!')
    times = times[1] + months(1:length(times)-1)
}
day(times) = 15
times = as.Date(times)
#check periodicity
if( !isTRUE(all(times == times[1] + months(1:length(times)-1))) ) stop(flog.fatal('Incorrect monthly periodicity in the input file'))
start_ym = c(year(times[1]), month(times[1]))
end_ym = c(year(times[length(times)]), month(times[length(times)]))

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

if (is.na(nthreads)) {
    nthreads = parallel::detectCores()
}
if (nthreads > 1 ) {
    p_load(parallel)
    cluster = makeCluster(nthreads)
    flog.info("Using %i threads", nthreads)
} else {
    cluster = NULL
}

flog.info("Starting computation")
spi_res <- st_apply(nc_in, 1:2, calc_SPI, CLUSTER = cluster, ts = timescale, thr = na_thr, ref.s = ref_start, ref.e = ref_end, first_ym=start_ym, last_ym=end_ym)
flog.info("Ended computation")

shift_vec = function(v) {
    len = length(v)
    return(c(v[2:len], v[1]))
}
spi_res = aperm(spi_res, shift_vec(1:length(dim(spi_res))))
spi_res[abs(spi_res) > spi_na_thr] = NA

#============= WRITE =============

flog.info("Creating output file by copying input")
invisible(file.copy(from = fn_in, to = fn_out, overwrite = TRUE))
nc_out = nc_open(fn_out, write = TRUE)
nc_out = ncvar_rename(nc_out, var_in, var_out)

flog.info("Writing output file")
sf = ncatt_get(nc_out, var_out, 'scale_factor')
ao = ncatt_get(nc_out, var_out, 'add_offset')
if (isTRUE(ao$hasatt)) {
    spi_res = spi_res - ao$value
}
if (isTRUE(sf$hasatt)) {
    spi_res = spi_res / sf$value
}
ncvar_put(nc_out, var_out, vals = spi_res[[names(spi_res)]])

h = ncatt_get(nc_out, 0, 'history')$value
ncatt_put(nc_out, 0, 'history', paste0(format(Sys.time(), "%D %H:%M %Z: "), "SPI index calculated with SPI_new.R script from Adriano Fantini (afantini@ictp.it)\n ", h))

nc_close(nc_out)
flog.info("Done")

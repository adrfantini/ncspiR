#!/usr/bin/env Rscript

program_name = 'ncspiR.R'
description = 'Script to calculate the SPI index from a monthly netCDF file containing precipitation data'
further_description = 'Input file MUST be monthly. Does everything in memory, so make sure your dataset fits in memory!
    Very few checks are performed, so make sure you know what you are doing.'
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
    'futile.logger'
)

#============= INITIALIZATION =============

suppressPackageStartupMessages(library(pacman))
options("pac_update" = FALSE) # Do not try to update relevant packages on package install with p_load
suppressPackageStartupMessages(p_load(optparse))
suppressPackageStartupMessages(p_load(glue))

option_list = list(make_option(c("-t", "--timescale"),
                                type="integer",
                                default=12,
                                help="Timescale in months [default: %default]"),
                    make_option(c("-n", "--nthreads"),
                                type="integer",
                                default=NA,
                                help="Number of threads to be used, NA for automatic detection [default: %default]"),
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
                    make_option("--spinathreshold",
                                type="double",
                                default=1e10,
                                help="SPI output values greater than this (in absolute value) will become NA. [default: %default]"),
                    make_option("--refstart",
                                type="character",
                                default=NULL,
                                help="ref.start parameter to pass to SPEI:spi, a character in the YEAR-MON format. [default: %default]"),
                    make_option("--refend",
                                type="character",
                                default=NULL,
                                help="ref.end parameter to pass to SPEI:spi, a character in the YEAR-MON format. [default: %default]")
                    )
parser = OptionParser(
    usage = "%prog [options] INPUT OUTPUT",
    option_list=option_list,
    epilogue=glue("
    #================= DESCRIPTION =================#
    {program_name} version {version} from {contact}

    {description}
    {further_description}

    For feature requests and bug reports, please get in touch on github:
    {gh_url}/issues
    ")
)

#Gather input arguments
arguments = parse_args(parser, positional_arguments = 2) # , args = c('filein.nc', 'fileout.nc')) , args = c('--refend', '2005-12', 'prOK_mon_1981-2050.nc', 'test.nc'))
opt = arguments$options

suppressPackageStartupMessages(p_load(futile.logger))

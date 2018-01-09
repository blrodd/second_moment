#!/opt/antelope/5.6/bin/python

""" Python wrapper for McGuire 2017 MATLAB second moment program """

import os
import sys
import signal

signal.signal(signal.SIGINT, signal.SIG_DFL)
sys.path.append(os.environ['ANTELOPE'] + "/data/python")
sys.path.append(os.environ['ANTELOPE'] + "/contrib/data/python")

import os
import re
import sys
import glob
import stat
import json
import inspect
import logging
import csv
import subprocess

from time import sleep
from optparse import OptionParser

import antelope.stock as stock

# Set path of matlab script
code_sub_folder = 'matlab'
matlab_code_folder = '/Users/rrodd/bin/second_moment' + '/' + code_sub_folder
sys.path.append( matlab_code_folder )
matlab_code = matlab_code_folder + '/' + 'runSJFex.m'

#
# Get command-line arguments
#
usage = "\n\tUsage:\n"
usage += "\t\tsecond_moment -vdxw --nofig [-p parameter file] [-s select] [-r reject] [-f filter] [-t tw] database orid \n"

parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
        help="verbose output", default=False)

parser.add_option("-d", "--debug", action="store_true", dest="debug",
        help="debug output", default=False)

parser.add_option("-x", "--debug_plot", action="store_true", dest="debug_plot",
        help="debug plot output", default=False)

parser.add_option("-i", "--interactive", action="store_true", dest="interactive",
        help="run in interactive mode", default=False)

parser.add_option("--no_figure", action="store_true", dest="no_figure",
        help="save plots", default=False)

parser.add_option("-w", action="store_true", dest="window",
        help="run on active display", default=False)

parser.add_option("-e", "--egf", action="store", type="int", dest="egf",
        help="egf orid", default=-99)

parser.add_option("-p", action="store", type="string", dest="pf",
        help="parameter file", default='second_moment.pf')

parser.add_option("-s", "--select", action="store", type="string", dest="select",
        help="station select", default=".*")

parser.add_option("-r", "--reject", action="store", type="string", dest="reject",
        help="station reject", default="")

parser.add_option("-f", "--filter", action="store", type="string", dest="filter",
        help="filter", default=None)

parser.add_option("-t", "--time_window", action="store", type="string", dest="tw",
        help="time window", default=None)

(options, args) = parser.parse_args()


# parse parameter file
pf = stock.pfread( options.pf )
pf_file = stock.pffiles( options.pf )[0]

if not os.path.isfile(pf_file):
    sys.exit( 'Cannot find parameter file [%s]' % options.pf )

# matlab inversion parameters
loaddatafile = float(pf['loaddatafile'])
domeas = float(pf['domeasurement'])
pickt2 = float(pf['pickt2'])
doinversion = float(pf['doinversion'])
dojackknife = float(pf['dojackknife'])
azband = float(pf['azband'])
dobootstrap = float(pf['dobootstrap'])
nb = float(pf['nb'])
bconf = float(pf['bconf'])
niter = float(pf['niter'])

# set up folders
image_dir = pf['image_dir']
temp_dir = pf['temp_dir']

# remove files in directory if it exists or create directory (add later)
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

# create image directory if it does exist
if not os.path.exists(image_dir):
    os.makedirs(image_dir)

# matlab info
matlab_path = pf['matlab_path']
matlab_flags = pf['matlab_flags']
xvfb_path = pf['xvfb_path']
matlab_nofig = pf['matlab_nofig']

# egf selection criteria
loc_margin = float(pf['location_margin'])
dep_margin = float(pf['depth_margin'])
time_margin = float(pf['time_margin'])

# filter and time window
if not options.filter:
    options.filter = pf['filter']

if not options.tw:
    options.tw = pf['time_window']


# L-curve time duration maximum
stf_duration_criteria = float(pf['misfit'])

#
# Start of main script
#

# Set log level
loglevel = 'WARNING'
if options.verbose:
    loglevel = 'INFO'
if options.debug:
    loglevel = 'DEBUG'

# All modules should use the same logging function. We have
# a nice method defined in the logging_helper lib that helps
# link the logging on all of the modules.

try:
    from logging_helper import getLogger
except Exception,e:
    sys.exit('Problems loading logging lib. %s' % e)

# New logger object and set loglevel
logging = getLogger(loglevel=loglevel)
logging.info('loglevel=%s' % loglevel)

logging.info( "Start: %s %s" % ( 'second_moment',  stock.strtime( stock.now() ) )  )
logging.info( "Start: configuration parameter file %s" % options.pf  )
logging.info( " - Xvfb path: %s" % xvfb_path  )
logging.info( " - Matlab path: %s" % matlab_path  )
logging.info( " - Matlab flags: %s" % matlab_flags  )

# Set virtual display if needed
if not options.window and xvfb_path:
    """
    Open virtual display. Running Xvfb from Python
    """

    pid = os.getpid()
    cmd = '%s :%s -fbdir /var/tmp -screen :%s 1600x1200x24' % (xvfb_path, pid, pid)

    logging.info( " - Start virtual display: %s" % cmd  )

    xvfb = subprocess.Popen( cmd, shell=True)

    if xvfb.returncode:
        stdout, stderr = xvfb.communicate()
        logging.info( " - xvfb: stdout: %s" % stdout  )
        logging.info( " - xvfb: stderr: %s" % stderr  )
        sys.exit('Problems on %s ' % cmd )

    os.environ["DISPLAY"] = ":%s" % pid

    logging.info( " - xvfb.pid: %s" % xvfb.pid  )
    logging.info( " - $DISPLAY => %s" % os.environ["DISPLAY"]  )

#
# Run Matlab code
#
cmd = "%s -r \"verbose='%s'; debug='%s'; debug_plot='%s'; interactive='%s'; no_figure='%s' \
                ; image_dir='%s'; temp_dir='%s'; db='%s'; orid=%d; egf=%d; reject='%s'; select='%s'; filter='%s' \
                ; tw='%s'; misfit_criteria=%.2f; loc_margin=%.4f; dep_margin=%.2f \
                ; time_margin=%.1f; LOADDATAFILE=%d; DOMEAS=%d; PICKt2=%d; DOINVERSION=%d; DOJACKKNIFE=%d \
                ; AZBAND=%d; DOBOOTSTRAP=%d; NB=%d; BCONF=%.2f; NITER=%d\" < '%s'"  \
                % (matlab_path, options.verbose, options.debug, options.debug_plot, options.interactive \
                , options.no_figure, image_dir, temp_dir, args[0], int(args[1]), options.egf, options.reject \
                , options.select, options.filter, options.tw, stf_duration_criteria \
                , loc_margin, dep_margin, time_margin, loaddatafile, domeas, pickt2, doinversion \
                , dojackknife, azband, dobootstrap, nb, bconf, niter, matlab_code)

logging.info( " - Run Matlab script:"  )
logging.info( "   %s " % cmd  )

def execute(command):
    """Executes the command to run matlab script."""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()
 
    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        raise ProcessException(command, exitCode, output)

try:
    mcmd = execute(cmd)
except Exception,e:
    print " - Problem on: [%s] " % cmd
    print " - Exception %s => %s" % (Exception,e)

if options.verbose:
    logging.info( "Done: %s %s" % ( 'second_moment',  stock.strtime( stock.now() ) )  )


#
# Kill virtual display if needed
#
if not options.window:
    logging.info( "xvfb.kill: %s" % xvfb.terminate() )

def num(s, r=None):
    if r:
        return round(float(s), r)
    else: 
        try:
            return int(s)
        except ValueError:
            return float(s)


    

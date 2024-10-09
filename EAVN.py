#!/usr/bin/env ParselTongue
# A pipeline script for EAVN data reduction
# Yuzhu Cui, Baoqiang Lao (lbq19881213@gmail.com): Sep 2024

from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData, AIPSImage
#from AIPSTV import AIPSTV
#import Wizardry.AIPSData
from EAVN_aips_tasks import *
from EAVN_functions import *
import EAVN_functions
import EAVN_aips_tasks
#import crfuncs
import traceback

import copy, optparse, os, sys
import re, string, pprint, math
import time
from datetime import datetime
import logging

c = 2.99792458e8

# define our local functions. More general purpose functions can be found in
# EAVN_functions.py. The actual AIPS task calls can be found in EAVN_aips_tasks.py.

def def_name(name):
    # shortcut to getting the default name for output plot files
    plotname = get_plotname(output_prefix, name)
    return plotname

def getfile(extension):
    # get calibration files from indir assuming standard names.
    file = indir + '/' + experiment + '.' + extension
    if not os.path.isfile(file):
        file = indir + '/' + experiment + '.' + extension
        #assert(os.path.isfile(file)), file + ' does not exist!'
    return file

def checkin(control):
    # convert the control hash to a list of variables of the appropriate type.
    # Should probably add some more validity checking here. Make these global
    # for ease of reference - they should never change.
    global indir, refantnames, plotrefnames, fits, outdir
    global nfits, heads, solint, experiment
    global bpass_calibrators, phaseref_sources, target_sources, avg, plotavg
    global tmask, fitsdir, fits_file, nselfcal, doplot, msgkill
    global source_select, experiment_dir, output_prefix, nbp_table, fring_snr
    global disk, freqid, dopng, rmcalsour, fit_ant, interferometer
    global eop_path, TECU_model, tecdir

    EAVN_aips_tasks.AIPSTask.version = control['version']['default']
    EAVN_aips_tasks.task_versions = control['version']

    experiment = control['experiment'][0].lower()
    # remove the pass number from the experiment name if required
    experiment_dir = get_experiment_dir(experiment)

    interferometer = control['interferometer'][0] 

    refantnames = control['refant']
    for i in range(len(refantnames)):
        refantnames[i] = refantnames[i].upper()

    plotrefnames = control['plotref']
    for i in range(len(plotrefnames)):
        plotrefnames[i] = plotrefnames[i].upper()

    fitsdir = control.get('fitsdir', [False])[0]
    if (not fitsdir):
        fitsdir = archive_fitsdir(experiment_dir)

    outdir = control.get('outdir', [False])[0]
    if (not outdir):
        outdir = os.environ.get('OUT') + '/' + experiment_dir

    indir = control.get('indir', [False])[0]
    if (not indir):
        indir = os.environ.get('IN') + '/' + experiment_dir

    nfits = int(control.get('nfits', [False])[0])

    # heads defaults to 1
    heads = control.get('glue_pass', [1])
    for i in range(len(heads)):
        heads[i] = int(heads[i])
    if (len(heads) >= 4):
        raise """Too many passes specified (max = 4) - check 'glue_pass' in input file"""

    # solint defaults to the typical scan length on the phaseref source
    solint = float(control.get('solint', [0])[0])

    AIPS.userno = int(control['userno'][0])

    bpass_calibrators = control.get('bpass', [])
    for i in range(len(bpass_calibrators)):
        bpass_calibrators[i] = bpass_calibrators[i].upper()
    nbp_table = 0
    if bpass_calibrators:
        nbp_table = 1

    phaseref_sources = control.get('phaseref', [])
    for i in range(len(phaseref_sources)):
        phaseref_sources[i] = phaseref_sources[i].upper()

    target_sources = control.get('target', [])
    for i in range(len(target_sources)):
        target_sources[i] = target_sources[i].upper()

    rmcalsour = control.get('rmcalsour', [])
    for i in range(len(rmcalsour)):
        rmcalsour[i] = rmcalsour[i].upper()

    fit_ant = control['fit_ant'][0].lower()

    avg = float(control.get('avg', [0])[0])

    plotavg = float(control.get('plotavg', [0])[0])

    tmask = control['tmask']
    for i in range (len(tmask)):
        tmask[i] = int(tmask[i])
    if (len(tmask) == 1):
        tmask.append(999)

    fring_snr = float(control.get('fring_snr', [7])[0])

    freqid = int(control.get('freqid', [False])[0])

    # nselfcal defaults to 2 but 0 is a valid input
    nselfcal = int(control.get('sciter', [2])[0])

    doplot = int(control.get('doplot', [True])[0])
    if doplot < 0:
        doplot = 0

    dopng = int(control.get('dopng', [False])[0])

    disk = int(control.get('disk', [1])[0])

    msgkill = int(control.get('msgkill', [0])[0])

    source_select = control.get('sources', [])
    for i in range(len(source_select)):
        source_select[i] = source_select[i].upper()

    assert (len(phaseref_sources) == len(target_sources)), """unmatched number 
                    of phase reference and target sources! """

    assert (os.path.exists(fitsdir)), 'FITS directory does not exist: ' + \
            fitsdir

    fits_file = control.get('fits_file', [])

    if (not solint and not phaseref_sources):
        raise 'You must set either solint or phaseref_sources'

    output_prefix = outdir + '/' + experiment

    eop_path = control['eop_path'][0]
    TECU_model = control['TECU_model'][0]
    tecdir = control['tecdir'][0]

def EAVN_checkdata():
    # extract useful information from the dataset

    logger.info('*' *20)
    logger.info('Checking the data:')
    # Find the reference antennae numbers (this should be a prioritised list)
    refantlist = []
    for antenna in refantnames:
        refantlist.append(get_ant_num(uvdata, antenna))
    logger.info('refant number= %s' % refantlist)

    # and the reference antennae for plotting
    plotref = []
    for antenna in plotrefnames:
        plotref.append(get_ant_num(uvdata, antenna))
    logger.info('plotref number= %s' % plotref)
    #raw_input('hit return')

    uvdata.zap_table('AIPS PL', -1)
    assert(uvdata.table_highver('AIPS PL') == 0)

    #sources = uvdata.sources
    sutable = uvdata.table('SU', 1)
    sources = dict()
    for row in sutable:
        source = row.source.strip()
        if (not source_select) or (source in source_select):
            sources[source] = dict()
            sources[source]['id'] = row.id__no
            sources[source]['scanlen'] = list()

    assert(len(sources) > 0)

    # find out which sources are not targets. Note that selfcal_sources is a
    # superset of phaseref_sources
    for source in sources:
        if not source in (target_sources):
            selfcal_sources.append(source)

    assert(len(selfcal_sources) > 0), 'no fringe calibrators'

    for calibrator in (selfcal_sources):
        assert(calibrator in sources), ('fring calibrator ' + calibrator +  
            ' not in source list')
    for calibrator in bpass_calibrators:
        assert(calibrator in sources), 'bpass calibrator ' + calibrator + ' not in source list'
    for target in target_sources:
        assert(target in sources), 'target source ' + target + ' not in source list'

    logger.info('target sources= %s' % target_sources)
    logger.info('phaseref sources= %s' % phaseref_sources)
    logger.info('selfcal sources= %s' % selfcal_sources)

    stokes = uvdata.stokes
    global plotstokes
    assert(len(stokes) > 0)
    if len(stokes) == 1:
        plotstokes = stokes[0]
    else:
        plotstokes = 'HALF'

    polarizations = uvdata.polarizations
    polarizations.sort()
    assert(len(polarizations) > 0)

    chanwidth = uvdata.header.cdelt[2]  # in Hz
    nchan = uvdata.header.naxis[2]
    nif = uvdata.header.naxis[3]
    ifwidth = chanwidth*nchan

    logger.info('IF bandwidth= %s %s' % (ifwidth*1.e-6, ' MHz'))
    logger.info('channel width= %s %s'% (chanwidth*1.e-3, ' kHz'))
    logger.info('number channels per IF= %s' % nchan)
    logger.info('number IFs= %s' % nif)
    logger.info('antennas=i %s' % uvdata.antennas)

    logger.info('fring SNR cutoff= %s' % fring_snr)

    return (refantlist, sources, selfcal_sources, plotref, ifwidth,
            chanwidth, nchan, nif, stokes, plotstokes, polarizations)    


def save_fits_image(data, name):
    fitsoutfile = outdir + '/' + experiment + '_' + name + '.FITS'

    # Rename any old plot files before writing the new one
    save_old_file(fitsoutfile)
    runfittp(indata=data, outfile=fitsoutfile)

def EAVN_load(uvdata, msortdata, cltablemin, wtthreshhold, interferometer):
    # load the data then sort, vbglu and uvavg if required.

    # first zap any of the 'load' files from a previous run (in case the names
    # that are used have changed from the previous run due to e.g. different
    # number of heads)
    zap_old_data(uvdata)
    zap_old_data(msortdata)
    zap_old_data(uvcopdata)
    zap_old_data(uvavgdata)
    
    global nfits
    uvdata, headdata = load_data(uvdata, msortdata, experiment, fitsdir,
            fits_file, heads, nfits, cltablemin, wtthreshhold, interferometer)

    # and run VBGLU if necessary
    if len(heads) > 1:
        print >>PYPELOG, '\nNeed to VBGLU the data'

        uvdata = glue_data(headdata, vbgludata, fxpoldata)

    # select the freqid if necessary
    if freqid:
        runuvcop(uvdata, uvcopdata, freqid=freqid)
        runindxr(uvcopdata)
        uvdata = uvcopdata


    # now do averaging if requested
    if avg > 0.1e-4:
        runuvavg(indata=uvdata, outdata=uvavgdata, avg=avg)
        runindxr(uvavgdata)
        uvdata = uvavgdata

def EAVN_flag(uvdata):
    # Delete any old flag tables.
    uvdata.zap_table('AIPS FG', -1)

    try:
        #uvflg.infile = uvflg_file
        logger.info("flagging using file: %s" % uvflg_file)
        runuvflg(indata=uvdata, infile=uvflg_file)
    except:
        raise "uvflg failed, check " + uvflg_file + " is o.k.!"

    # run vlog if necessary
    doglobal, vlba_missing = vlba_logs(uvdata, vlog_out, vlbacal_file,
            ifwidth)

    # run uvflg for any VLBA antennas.
    if doglobal:
        vlba_uvflg_file = vlog_out + '.FLAG'
        try:
            logger.info("flagging VLBAs using file: %s" % vlba_uvflg_file)
            runuvflg(indata=uvdata, infile=vlba_uvflg_file)
        except:
            raise "uvflg failed, check " + vlba_uvflg_file + " is o.k.!"


    if os.path.exists(chflg_file):
        runuvflg(indata=uvdata, infile=chflg_file)
    else:
        # flag the edge channels if no chflag file given
        logger.info(chflg_file + " is missing")

        nflag = min(nchan//16, 8)
        logger.info('Flagging outer %d' % nflag + 'channels instead')
        reason = 'subband edge'
        runuvflg(indata=uvdata, bchan=0, echan=nflag, reason=reason)
        runuvflg(indata=uvdata, bchan=nchan+1-nflag, echan=nchan, reason=reason)

def EAVN_plot1():
    # scan list
    outfile = output_prefix + '.SCAN'
    save_old_file(outfile)        
    runlistr(indata=uvdata, optype='SCAN', outprint=outfile)

    # list of stations and number of vis.
    outfile = output_prefix + '.DTSUM'
    save_old_file(outfile)        
    rundtsum(indata=uvdata, outprint=outfile)

    if doplot:
        # possm plot of the autocorrelations 
        aparm = set_default_aparms('possm')
        # aparm(7) gives freq. labelling, aparm(8) gives autocorrelation
        aparm[7] = 1
        aparm[8] = 1
        #print 'aparm=',aparm

        uvdata.zap_table('AIPS PL', -1)
        runpossm(indata=uvdata, aparm=aparm, stokes=plotstokes)
        plot(uvdata, def_name('POSSM_AUTOCORR'), dopng=dopng)


        # vplot of the raw data
        bparm = set_default_bparms('vplot')
        uvdata.zap_table('AIPS PL', -1)
        for poln in polarizations:
            for incif in range(1, nif+1):
                error = runvplot(indata=uvdata, bparm=bparm, antenna=plotref,
                        bif=incif, eif=incif, docalib=0, stokes=poln+poln,
                        solint=plotavg)
                if (error):
                    print_vplot_err(incif, incif, poln+poln, uvdata)
        plot(uvdata, def_name('VPLOT_UNCAL'), dopng=dopng)

        # possm plot of the cross-correlations 
        aparm = set_default_aparms('possm')
        uvdata.zap_table('AIPS PL', -1)
        runpossm(indata=uvdata, aparm=aparm, stokes=plotstokes,
            antennas=plotref)
        plot(uvdata, def_name('POSSM_UNCAL'), dopng=dopng)

        # possm plot of the cross-hand cross-correlations 
        if len(stokes) > 2:
            aparm = set_default_aparms('possm')
            uvdata.zap_table('AIPS PL', -1)
            runpossm(indata=uvdata, aparm=aparm, stokes='LR', antennas=plotref)
            plot(uvdata, def_name('POSSM_CPOL'), dopng=dopng)

def EAVN_plot2a(uvdata):
    uvdata.zap_table('AIPS PL', -1)

    # plot the bandpass table
    if nbp_table:
        aparm = set_default_aparms('possm')
        #aparm[4] = 1.3
        #aparm[5:7] = [0, 0]
        aparm[8] = 2
        #print >>PYPELOG, 'aparm=', aparm
        runpossm(indata=uvdata, aparm=aparm, stokes=plotstokes, solint=0, bpver=1, gainuse=2)
        plot(uvdata, def_name('BANDPASS'), dopng=dopng)

    # plot the fringe solutions
    #uvdata.zap_table('AIPS PL', -1)
    #runsnplt(uvdata, 'CL', 3, 'PHAS')
    #plot(uvdata, def_name('FRING_PHAS'), dopng=dopng)

    #uvdata.zap_table('AIPS PL', -1)
    #runsnplt(uvdata, 'SN', 2, 'DELA')
    #plot(uvdata, def_name('FRING_DELAY'), dopng=dopng)

    #uvdata.zap_table('AIPS PL', -1)
    #runsnplt(uvdata, 'SN', 2, 'RATE')
    #plot(uvdata, def_name('FRING_RATE'), dopng=dopng)


def EAVN_ampcal(antab_file, uvdata):
    # run ACCOR, CLCAL, ANTAB and APCAL. Also VLOG for any EAVN antennas
    uvdata.zap_table('AIPS SN', -1)
    while uvdata.table_highver('AIPS CL') > 1:
        uvdata.zap_table('AIPS CL', 0)
    uvdata.zap_table('AIPS TY', -1)
    uvdata.zap_table('AIPS GC', -1)

    #run ACCOR --> SN 1
    if interferometer == 'EAVN':
       solint = -1
    elif interferometer == 'VLBA':
       solint = -2
    else:
       solint = 0 
    runaccor(uvdata, solint)    
    uvdata.zap_table('AIPS PL', -1)
    runsnplt(uvdata, 'SN', 1, 'AMP') 
    plot(uvdata, def_name('SN1'), dopng=dopng)

    #run CLCAL --> CL 2
    runclcal(gainver=1, gainuse=2, indata=uvdata, snver=1,
            refant=refantlist[0], interpol='self', doblank=1, dobtween=1,
            samptype='', inver=0) #refant=refantlist[0]
    uvdata.zap_table('AIPS PL', -1)
    aparm = [1,0,0,0,-180,180,0,0,0,0]
    runpossm(indata=uvdata, aparm=aparm, solint=-1, freqid=1, sources=[],
        timerang=[0,0,0,1], codetype='A&P',
        nplots=8, stokes='HALF', antennas=[], docalib=1, gainuse=2,
        doband=0, bpver=0, flagver=0, dotv=-1)
    plot(uvdata, def_name('CL2'), dopng=dopng)
   
    #run ANTAB
    logger.info('Antab file= %s' % antab_file)
    runantab(antab_file, uvdata)


    doglobal, sparm = global_exper(uvdata)
    if doglobal:
        # if a global then we must run ANTAB again, but add missing VLBA
        # antennas to SPARM in ANTAB (uggh!).

        antab_file = vlog_out + '.TSYS'
        logger.info('Antab file= %s' % antab_file)
        runantab(antab_file, uvdata, sparm=sparm)


    uvdata.zap_table('AIPS PL', -1)
    runsnplt(uvdata, 'TY', 1, 'TSYS')
    plot(uvdata, def_name('TSYS'), dopng=dopng)

    #Generate SN 2-4 table
    ant_id = int(get_ant_num(uvdata, fit_ant))
    len_ants = len(uvdata.antennas)
    runapcal(indata=uvdata, tyver=1, gcver=1, snver=2, freqid=1, opcode='grid', ant_id=ant_id, len_ants=len_ants)
    uvdata.zap_table('AIPS PL', -1)
    runsnplt(uvdata, 'SN', 2, 'AMP')
    plot(uvdata, def_name('SN2'), dopng=dopng)
    runapcal(indata=uvdata, tyver=1, gcver=1, snver=3, freqid=1, opcode='opac', ant_id=ant_id, len_ants=len_ants)
    uvdata.zap_table('AIPS PL', -1)
    runsnplt(uvdata, 'SN', 3, 'AMP')
    plot(uvdata, def_name('SN3'), dopng=dopng)
    runapcal(indata=uvdata, tyver=1, gcver=1, snver=4, freqid=1, opcode='lesq', ant_id=ant_id, len_ants=len_ants)
    uvdata.zap_table('AIPS PL', -1)
    runsnplt(uvdata, 'SN', 4, 'AMP')
    plot(uvdata, def_name('SN4'), dopng=dopng)
    # use 'box', dobtween and doblank so no sources get left out
    runclcal(gainver=2, gainuse=3, indata=uvdata, snver=3,
            refant=refantlist[0], interpol='self', doblank=1, dobtween=1,
            samptype='', inver=3)

    uvdata.zap_table('AIPS PL', -1)
    runsnplt(uvdata, 'CL', 3, 'AMP')
    plot(uvdata, def_name('GAIN'), dopng=dopng)

    # Do the parallactic angle correction.
    #runclcor(indata=uvdata, clcorprm=[0,1], opcode='PANG', gainver=2, gainuse=2)            

def EAVN_fring(fringdata, rmcalsour):
    # run fring (on averaged data set if necessary)

    aparm = []
    aparm[1:] = [0 for i in range(10)]
    aparm[6] = 2
    aparm[7] = 5
    dparm = []
    dparm[1:] = [0 for i in range(10)]
    dparm[2] = 100
    dparm[3] = 100
    #calsours = selfcal_sources
    #calsours.remove('RT-VIR') 
    rmcalsour_str = '-%s' % rmcalsour[0]
    calsours = [rmcalsour_str]
    runfring(indata=fringdata, snver=5, gainuse=3, refant=refantlist[0], 
            solint=0.5, calsour=calsours,
            aparm=aparm, dparm=dparm, snr=fring_snr)

    uvdata.zap_table('AIPS PL', -1)
    runsnplt(uvdata, 'SN', 5, 'AMP')
    plot(uvdata, def_name('SN5'), dopng=dopng)
    # apply the fring solutions CL3=>CL4
    calibrator = phaseref_sources
    target = target_sources
    #for (calibrator, target) in ( zip (phaseref_sources, target_sources) +
    #                zip (selfcal_sources, selfcal_sources)):
    logger.info('refant is %d' % refantlist[0])
    runclcal(indata=uvdata, opcode='CALI', interpol='2pt', snver=5,
                gainver=3, gainuse=4, calsour=[], sources=[],
                refant=refantlist[0], inver=5) #refantlist[0]
    uvdata.zap_table('AIPS PL', -1)
    aparm = [0,0,0,0,-180,180,0,0,0,0]
    runpossm(indata=uvdata, aparm=aparm, solint=-1, freqid=1, sources=calibrator[0],
        timerang=[0,0,0,1], codetype='A&P',
        nplots=7, stokes='HALF', antennas=[], docalib=1, gainuse=4,
        doband=0, bpver=1, flagver=0, dotv=-1)
    plot(uvdata, def_name('CL4'), dopng=dopng)


#Start main function
log_file = 'aips_{}.log'.format(datetime.now().strftime('%Y%m%d%H%M%S'))
logging.basicConfig(level=logging.INFO,#logging.DEBUG, #logging.INFO
                    filename=log_file,
                    datefmt='%Y/%m/%d %H:%M:%S',
                    format='%(asctime)s - %(name)s - %(levelname)s - %(lineno)d - %(module)s - %(message)s'
                    )
                   #handlers=[
                   #     logging.FileHandler("logfile.log"),
                   #     logging.StreamHandler(sys.stdout)
                   # ]
logger = logging.getLogger(__name__)
file_handler = logging.FileHandler(log_file)
logger.addHandler(file_handler)
logger.info('Start EAVN pipeline processing.')

usage = 'usage: ParselTongue %prog [options] template.inp'
parser = optparse.OptionParser(usage=usage, version='%prog 31DEC24')
(options, args) = parser.parse_args()
if len(args) != 1:
    parser.error("incorrect number of arguments")
today = time.asctime()   

control = parse_inp(args[0])

# check the inputs and re-type where necessary
checkin(control)

#if interferometer == 'VLBA' or interferometer == 'EVN':
#   uvflg_file = getfile('uvflg')
#   assert(os.path.isfile(uvflg_file)), uvflg_file + ' does not exist!'

if interferometer == 'EAVN' or interferometer == 'EVN':
   # get some necessary files assuming standard names
   antab_file = getfile('antab')
   assert(os.path.isfile(antab_file)), antab_file + ' does not exist!'

uvcopklass = 'FQ' + str(freqid)
uvname = experiment.upper()
uvdata = AIPSUVData(uvname, 'UVDATA', disk, 1)
msortdata = AIPSUVData(uvdata.name, 'MSORT', uvdata.disk, 1)
uvcopdata = AIPSUVData(uvdata.name, uvcopklass, uvdata.disk, 1)
uvavgdata = AIPSUVData(uvdata.name, 'UVAVG', uvdata.disk, 1)

# 1, Load and sort the data - average too if requested
if tmask[0] <= 1 <= tmask[1]:
    logger.info('Starting tmask 1: load and sort the data')
    if interferometer == 'EAVN':
       EAVN_load(uvdata, msortdata, 1./60., 0, interferometer)
    elif interferometer == 'VLBA':
       EAVN_load(uvdata, msortdata, 0.25, 0.7, interferometer)
    elif interferometer == 'EVN':
       EAVN_load(uvdata, msortdata, 1./60., 0, interferometer)
    else:
       raise ValueError('Please input the correct interferometer name: EAVN, VLAB or EVN!')   
    logger.info('Ending tmask 1')

# Always do this bit to make sure names are right at end of block. The order of
# the if statements here is significant.
if msortdata.exists():
    uvdata = msortdata

logger.info('uvdata= %s' % uvdata)

# check the data and extract useful information (freq, shape, etc.). Initialise
# some things as well
sources = []
selfcal_sources = []
plotref = []
nchan = None
nif = None
chanwidth = None
ifwidth = None
plotstokes = None


(refantlist, sources, selfcal_sources, plotref, ifwidth, chanwidth,
        nchan, nif, stokes, plotstokes, polarizations) = EAVN_checkdata()

# get the source list and scan details from the NX table
sources = read_nx(uvdata, sources)

if (not solint):
    # solint defaults to the typical scan length on the phase reference
    solint = getsolint(sources, phaseref_sources)

# Download TEC maps and EOPs
if interferometer == 'VLBA':
   (year, month, day)=get_observation_year_month_day(uvdata)
   num_days = get_number_days_observations(uvdata)

   doy=get_day_of_year(year, month, day)
    
   get_TEC(year,doy,TECU_model,tecdir)
   if not os.path.exists(eop_path):
      os.mkdir(eop_path)
   get_eop(eop_path)

   if num_days==2: get_TEC(year,doy+1,TECU_model, tecdir)

# 2, plot the raw data, scan list, summaries of entire data, integration time, visibilities on each baseline
if tmask[0] <= 2 <= tmask[1]:
    logger.info('Starting tmask 2: plot the data - vs time and frequency')
    table_vers(uvdata=uvdata, cl=1, sn=0, fg=0, bp=0)

    EAVN_plot1()
    logger.info('Ending tmask 2')

if interferometer == 'EAVN':
   # 3, Autocorrelation and amplitude calibration
   if tmask[0] <= 3 <= tmask[1]:
       logger.info('Starting tmask 3: autocorrelation and amplitude calibration')
       table_vers(uvdata=uvdata, cl=1, sn=0, fg=0, bp=0)
   
       EAVN_ampcal(antab_file, uvdata)
       logger.info('Ending tmask 3')
   
   # 4, Bandpass calibration
   if tmask[0] <= 4 <= tmask[1] and nbp_table:
       logger.info('Starting tmask 4: bandpass calibration')
       table_vers(uvdata=uvdata, cl=3, sn=3, fg=0, bp=0)
   
       runbpass(refant=refantlist[0],indata=uvdata, calsour=bpass_calibrators) #refant=refantlist[0]
       # Do the less time consuming plots now
       EAVN_plot2a(uvdata)
       logger.info('Ending tmask 4')
   
   # 5, Fringe fitting.
   if tmask[0] <= 5 <= tmask[1]:
       logger.info('Starting tmask 5: fringe fitting')
       table_vers(uvdata=uvdata, cl=3, sn=3, fg=0, bp=nbp_table)
       
       EAVN_fring(uvdata, rmcalsour)
       logger.info('Ending tmask 5')
   
   # 6, Split
   if tmask[0] <= 6 <= tmask[1]:
       logger.info('Starting tmask 6: split the calibrated data')
       table_vers(uvdata=uvdata, cl=4, sn=5, fg=0, bp=nbp_table)
       sources = list(sources)
       for source in sources:
           splitdata = AIPSUVData(source, 'SPLIT', uvdata.disk, 1)
           zap_old_data(splitdata)
   
       runsplit(sources=sources, indata=uvdata, gainuse=4, docalib=1,
               doband=nbp_table, bpver=nbp_table) #, outseq=1)
   
       logger.info('Ending tmask 6')
   
   # 7, Save the split data as fits
   if tmask[0] <= 7 <= tmask[1]:
       logger.info('Starting tmask 7: split the calibrated data')
       for source in sources:
           splitdata = AIPSUVData(source, 'SPLIT', uvdata.disk, 1)
           fitsoutfile = output_prefix + '_' + source + \
                               '.UVDATA.FITS'
           save_old_file(fitsoutfile)
           runfittp(indata=splitdata, outfile=fitsoutfile)
   
       logger.info('Ending tmask 7')

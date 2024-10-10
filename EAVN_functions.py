"""
Created on Tue Sep  3 21:19:08 2024

@author: Yuzhu Cui, Baoqiang Lao
"""

from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData, AIPSImage
import copy, optparse, os, sys
import re, string, pprint, math
import numpy as np
from EAVN_aips_tasks import *
import EAVN_aips_tasks
import logging

logger = logging.getLogger(__name__)

#default_ver = '31DEC24'
#alt_ver = '31DEC04'
#AIPSTask.version = default_ver

# set the output to be stdout by default. Scripts that call these routines can
# then override this.
FuncLog = sys.stdout
pid_list = []


def save_old_file(filename, save_old=False):
    # check if an external file exists and delete/rename if it does
    if os.path.isfile(filename):
        if save_old:
            oldfilename = filename
            uniq = 0
            while os.path.isfile(oldfilename):
                oldfilename = filename + '.' + str(uniq)
                uniq += 1
            os.rename(filename, oldfilename)
            output = "renaming " + str(filename) + " to " + str(oldfilename)
        else:
            os.remove(filename)
            output = "deleting " + str(filename)
        logger.info(output)


def create_image(uvdata, imgdata, niter, uvwtfn, dotv=-1, stokes='I',
        imsize=512, dopol=-1, robust=0, ccmin=-1e-5):
    # run imagr and kntr

    beamklass = stokes + 'BM001'
    beamdata = AIPSImage(imgdata.name, beamklass, imgdata.disk, imgdata.seq)
    zap_old_data(imgdata)
    zap_old_data(beamdata)

    cellsize = getcell(uvdata)
    print >>FuncLog, 'number pixels=', imsize
    print >>FuncLog, 'pixel size=', cellsize*1.e3, 'mas'
    nchan = uvdata.header.naxis[2]

    runimagr(indata=uvdata, source=uvdata.name, cellsize=cellsize, niter=niter,
            outdata=imgdata, uvwtfn=uvwtfn, nchan=nchan, dotv=dotv,
            stokes=stokes, imsize=imsize, dopol=dopol, ccmin=ccmin,
            robust=robust)


    #zap_old_data(beamdata)

    # make plot file of kntr map, also putting the kntr map on the TV if
    # working interactively:
    peakflux, rmsflux = plotmap(imgdata)
    if dotv > 0:
        peakflux, rmsflux = plotmap(imgdata, dotv=dotv)

    print >>FuncLog, 'peak flux=', peakflux
    print >>FuncLog, 'rms flux=', rmsflux
    print >>FuncLog, 'dynamic range=', peakflux/rmsflux

    return peakflux, rmsflux



def zap_old_data(uvdata):
    # delete AIPS data from catalogue
    if uvdata.exists():
        logger.info('zapping old data: aipsuvname(uvdata')
        uvdata.zap(force=True)
        #s = raw_input('return')
    else:
        logger.info('no old data to zap: aipsuvname(uvdata)')



def get_ant_num(uvdata, refant_name):
    # convert antenna name to number
    refant_num = 0
    for antenna in uvdata.antennas:
        refant_num +=  1
        if antenna.upper() == refant_name.upper():
            #print antenna + '=' + str(refant_num)
            break
    else:
        raise(RuntimeError, "Reference antenna " + refant_name + " not found!")

    return refant_num


def getcell(uvdata=None):
    # get cell size for image pixels
    uvmax = uvdata.keywords['UVPRAMAX']
    cellsize = (1./uvmax) * 206265.
    cellsize = cellsize/4.
    return cellsize


def table_vers(uvdata=None, cl=None, sn=None, fg=None, bp=None):
    '''make sure we have the correct version of various tables.'''

    #assert (cl != None, sn != None, fg != None, bp != None), '''set all tables
    #            in table_vers'''

    while cl != None and uvdata.table_highver('AIPS CL') > cl:
        uvdata.zap_table('AIPS CL', 0)
    while sn != None and uvdata.table_highver('AIPS SN') > sn:
        uvdata.zap_table('AIPS SN', 0)
    while fg != None and uvdata.table_highver('AIPS FG') > fg:
        uvdata.zap_table('AIPS FG', 0)
    while bp != None and uvdata.table_highver('AIPS BP') > bp:
        uvdata.zap_table('AIPS BP', 0)

    assert(fg == None or uvdata.table_highver('AIPS FG') == fg)
    assert(sn == None or uvdata.table_highver('AIPS SN') == sn)
    assert(cl == None or uvdata.table_highver('AIPS CL') == cl)
    assert(bp == None or uvdata.table_highver('AIPS BP') == bp)

def read_nx(uvdata, sources):
    # extract useful information from the nx table

    #global sources
    nxtable = uvdata.table('NX', 1)
    for row in nxtable:
        for source in sources:
            if sources[source]['id'] == row.source_id:
                # convert from days to minutes
                scanlen = np.multiply(row.time_interval, 24.*60.)
                sources[source]['scanlen'].append(scanlen)

    return sources

def getsolint(sources, phaseref_sources):
    # determine a reasonable solution interval
    #global solint
    calscans = list()
    for source in phaseref_sources:
        calscans = calscans + sources[source]['scanlen']

    solint = np.median(calscans)
    logger.info('median calibrator scan length = %s' % solint)
    if len(calscans) > 5:
        # add a couple of standard deviations to solint
        stdev = np.std(calscans)
        stderr = stdev/math.sqrt(len(calscans))
        solint = solint + stderr
        logger.info('stderr=%s' % stderr)
    #print 'solint=', solint
    return solint

def parse_inp(filename):
    # form a hash of parameter names and values parsed from an input file.
    # Don't worry about parameter types (lists, strings, etc.) as we sort that
    # out after
    INPUTFILE = open(filename, "r")
    control = dict()

    # a few useful regular expressions
    newline = re.compile(r'\n')
    space = re.compile(r'\s')
    char = re.compile(r'\w')
    comment = re.compile(r'#.*')

    # parse the input file assuming '=' is used to separate names from values
    for line in INPUTFILE:
        if char.match(line):
            line = comment.sub(r'', line)
            line = line.replace("'", '')
            (param, value) = line.split('=')

            param = newline.sub(r'', param)
            param = param.strip()
            param = space.sub(r'', param)

            value = newline.sub(r'', value)
            value = value.strip()
            value = space.sub(r'', value)
            valuelist = value.split(',')
            #print param,'=',valuelist
            control[param] = valuelist
    # Move task-version parameters into a single dictionary
    version_dict = dict()
    if 'version' in control:
        version_dict['default'] = control['version'][0]
        del control['version']
    else:
        version_dict['default'] = EAVN_aips_tasks.AIPSTask.version
    for param in control.keys():
        if param.startswith('version'):
            version_dict[param.split('_')[1]] = control[param][0]
            del control[param]
    control['version'] = version_dict
    return control

def plot_closure(uvdata, imgdata, source, nselfcal):
    # subtract model from data and plot closure phases of result

    splitdata = AIPSUVData(uvdata.name, 'SPLIT', uvdata.disk, 2)
    zap_old_data(splitdata)

    runsplit(sources=[source], indata=uvdata, gainuse=nselfcal+1,
            doband=0, outseq=2)

    uvsubdata = AIPSUVData(uvdata.name, 'UVSUB', uvdata.disk, 1)
    zap_old_data(uvsubdata)
    runuvsub(indata=splitdata, in2data=imgdata)

    uvdata.zap_table('AIPS PL', -1)
    runclplt(indata=uvsubdata)
    return uvsubdata

def plotmap(imgdata, dotv=-1):
    # produce a contour plot of a map

    peakflux, rmsflux = get_dyn_range(imgdata)

    #imgdata.zap_table('AIPS PL', -1)

    runkntr(indata=imgdata, rmsflux=rmsflux, peakflux=peakflux, dotv=dotv)

    return peakflux, rmsflux


def get_dyn_range(imgdata):
    x = imgdata.header.naxis[0]
    y = imgdata.header.naxis[1]
    (peakflux, pixstd) = runimean(imgdata, [0, 0], [x // 4, y // 4])
    rmsflux = 3. * pixstd
    if (rmsflux < 1.e-5):
        rmsflux = 1.e-5
    if (rmsflux > peakflux):
        rmsflux = peakflux/4.

    return peakflux, rmsflux

def archive_fitsdir(experiment):
    # get the location of fits files from a given experiment on the JIVE
    # archive.

    #archive_base = '/jop15_1/archive/exp/'
    archive_base = os.environ.get('PIPEFITS')
    if (not archive_base):
        raise """environment variable $PIPEFITS is not set - please set or
                specify fitsdir in the input file"""
    exper_date_list = os.listdir(archive_base)
    for exper_date in exper_date_list:
        if experiment.upper() in exper_date:
            fitsdir = archive_base + '/' + exper_date + '/fits/'
            break
    else:
        raise 'Experiment directory not found on archive'

    return fitsdir

def get_experiment_dir(experiment):
    # remove the pass number from the end of the experiment name. This should
    # give the name of the experiment directories
    pass_number = re.compile('_\w+$')
    experiment = pass_number.sub(r'', experiment)
    return experiment

def get_nfits(fits_base):
    # check the fits files in the given fits directory and with the correct
    # name to determine number of fits files that must be loaded

    all_fitsfiles = os.listdir(os.path.dirname(fits_base))

    nfits=0
    for fits in all_fitsfiles:
        if os.path.split(fits_base)[1] in fits:
            nfits = nfits+1

    return nfits


def selfcal_map(uvdata, source, nsc, solint, refant, ncc=300, prefix='',
        imgseq=1, dopng=False):

    # iterate on selfcal and mapping steps

    for i in range(nsc):
        selfcal_solint = solint/((i//2)+1)
        imgdata = AIPSImage(uvdata.name, 'ICL001', uvdata.disk, imgseq)
        if (i%2 == 0):
            solmode = 'P'
            plottype = 'PHAS'
        else:
            solmode = 'A&P'
            plottype = 'AMP'

        print >>FuncLog, "selfcal iteration ", i+1, \
                ' using solmode= ', solmode, ' and source model ', \
                aipsuvname(imgdata)
        print >>FuncLog, "selfcal solint=", selfcal_solint


        runcalib(indata=uvdata, in2data=imgdata, docalib=2, ncomp=ncc,
                snr=5, solint=selfcal_solint, calsour=uvdata.name,
                solmode=solmode, refant=refant, snver=0)
        sntable = uvdata.table_highver('AIPS SN')
        if (solmode == 'A&P'):
            runsnsmo(indata=uvdata, refant=refant, invers=sntable)

        runclcal(indata=uvdata, opcode='CALI', interpol='2PT',
                snver=sntable, gainver=0,
                gainuse=0, calsour=[uvdata.name], sources=[], refant=refant)

        uvdata.zap_table('AIPS PL', -1)
        error = runsnplt(indata=uvdata, inext='SN', invers=0, optype=plottype)
        if error:
            print >>FuncLog, '\nPype warning, SNPLT failed for SN= ', i+1, \
                    ', source= ', source
            print >>FuncLog, '\nPype warning, aborting the mapping for'  \
                    ' source: ', source
            break

        plotname = get_plotname(prefix, source + '_CALIB_' + plottype +
                str(imgseq))

        plot(uvdata, plotname, dopng=dopng)
        # get a text version of the SN table if it's amp cal
        if solmode == 'A&P':
            ampfile = (prefix + '_' + source +
                    '_CALIB_AMP' + str(i+1) + '.TXT')
            save_old_file(ampfile)
            runlistr(indata=uvdata, inext='SN', inver=0,
                    optype='GAIN', outprint=ampfile)

        # next round of imaging
        imgseq += 1
        imgdata = AIPSImage(uvdata.name, 'ICL001', uvdata.disk, imgseq)
        print >>FuncLog, 'Creating image:', aipsuvname(imgdata)
        create_image(uvdata, imgdata, ncc, 'N', dotv=-1)
        plotname = get_plotname(prefix, source + '_ICLN_' + str(imgseq))
        plot(imgdata, plotname, dopng=dopng)


def plot(data, plotname, dopng=False):
    # standard plotting routine

    # Rename any old plot files before writing the new one
    save_old_file(plotname)
    logger.info('Creating plot file: %s' % plotname)
    runlwpla(indata=data, outfile=plotname)

    #plot_base = os.path.splitext(plotname)[0]
    logger.info('Converting %s to pdf.' % plotname)
    if dopng:
        logger.info('Converting %s to png.' % plotname)

    # convert from PS. Conversion can be quite slow so do in parallel while the
    # pipeline continues.
    global pid_list
    pid_list = pid_cleanup(pid_list)
    pid = os.fork()
    if pid == 0:
        plot_convert(plotname, dopng)
        os._exit(0)

    pid_list.append(pid)
    #print >>FuncLog, 'pid_list = ', pid_list

    # purge the plot files - this should be done again before starting the next
    # plot as plot files could get left around if program is aborted while,
    # say, POSSM is running
    #data.zap_table('AIPS PL', -1)


def plot_convert(plotname, dopng):
    # convert the postscript files to pdf (and png in some cases). Could use
    # 'convert', but ps2pdf seems to produce smaller files with the same
    # quality.

    pdfname = os.path.splitext(plotname)[0] + '.pdf'
    save_old_file(pdfname)
    os.system('ps2pdf ' +  plotname + ' ' + pdfname)
    #print 'Finished converting to pdf: ', pdfname
    #pid = os.spawnv(os.P_NOWAIT, 'ps2pdf', ['plotname', 'pdfname'])

    # and convert to png if requested.
    if dopng:
        pngname = os.path.splitext(plotname)[0] + '.png'
        save_old_file(pngname)
        #print >>FuncLog, 'Converting plot file to png: ', pngname
        os.system('convert ' + plotname + ' ' + pngname)
        #print 'Finished converting to png: ', pngname
        #pid = os.spawnv(os.P_NOWAIT, 'convert', ['plotname', 'pngname'])

    # delete the postscript
    logger.info('Deleting the postscript file: %s' % plotname)
    os.remove(plotname)

def pid_cleanup(pid_list):
    # clean up child processes, and return list of child processes still
    # running
    pid_notdone = []
    for pid in pid_list:
        #print >>FuncLog, 'pid=', pid
        wait_pid, status = os.waitpid(pid, os.WNOHANG)
        #print >>FuncLog, 'wait_pid, status=', wait_pid, status
        if wait_pid != pid:
            pid_notdone.append(pid)

    return pid_notdone


def get_plotname(prefix, name):
    # form standard name for plot files.
    plotname =  prefix + '_' + name + '.PS'
    return plotname

def vlba_logs(uvdata, vlog_out, vlbacal_file, ifwidth):
    # if a global then we must run VLOG
    doglobal, vlba_missing = global_exper(uvdata)
    if doglobal:
        assert(os.path.isfile(vlbacal_file)), vlbacal_file + ' does not exist!'

        print >>FuncLog, 'vlog file= ', vlbacal_file
        # delete any old versions of VLOG output
        vlog_cleanup(vlog_out)

        runvlog(indata=uvdata, infile=vlbacal_file, outfile=vlog_out,
                ifwidth=ifwidth)

        # vlog always leaves the dataset in READ status.
        uvdata.clrstat()

    return doglobal, vlba_missing

def global_exper(uvdata):
    # Check if its a global and record missing VLBA antennas to pass to ANTAB
    # (uggh!). GBT counts as a VLBA antenna.

    vlba_ants = ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC','GB']
    vlba_missing = []
    doglobal = None
    for vlba in vlba_ants:
        if vlba in uvdata.antennas:
           doglobal = 1
        else:
           vlba_missing.append(vlba)

    return doglobal, vlba_missing

def vlog_cleanup(vlog_out):
    vlog_files = [ 'FLAG', 'MKIII', 'PCAL', 'SCAN', 'TSYS', 'WX' ]
    for file in vlog_files:
        save_old_file(vlog_out + '.' + file)



def glue_data(headdata, vbgludata, fxpoldata=None):
    '''Take a list of files and vbglu them. Also run FXPOL if necessary.'''

    # sort the headdata files so that the largest is the first input to VBGLU
    headdata.sort(lambda x, y: len(y) - len(x))
    headpols = []
    for headnum in range(len(headdata)):
        if headdata[headnum]:
            headpols.append( headdata[headnum].polarizations )

    print >>FuncLog, 'VBGLU input files= ', headdata
    print >>FuncLog, 'pols for each head= ', pprint.pformat(headpols)

    runvbglu(indata=headdata[0], in2data=headdata[1], in3data=headdata[2],
            in4data=headdata[3], outdata=vbgludata)
    runindxr(vbgludata)
    uvdata = vbgludata

    # if the 'heads' were two pols of the same IFs, then need to run FXPOL
    if headpols[0] != headpols[1]:
        bandpol = "*(" + headpols[0][0] + headpols[1][0] + ")"
        print >>FuncLog, '\nVBGLU inputs appear to have different pols'
        print >>FuncLog, 'Running FXPOL with bandpol= ', bandpol
        runfxpol(indata=vbgludata, outdata=fxpoldata, bandpol=bandpol)
        runindxr(fxpoldata)
        uvdata = fxpoldata

    return uvdata


def load_data(uvdata, msortdata, experiment, fitsdir, fits_file=None, heads=[],
        nfits=None, cltablemin = 1./60., wtthreshhold=0, interferometer='EAVN'):
    uvname = uvdata.name
    #head_klass = []
    headdata = ['' for i in range(4)]
    i = -1
    for headnum in heads:
        i = i+1
        headname = uvname
        if (len(heads) > 1):
            headname += '_' + str(headnum)
        uvdata = AIPSUVData(headname, 'UVDATA', uvdata.disk, 1)
        logger.info('uvdata.exists = %s' % uvdata.exists())
        zap_old_data(uvdata)
        head_klass = 'UVDATA'

        if not fits_file:
            fits_file_base = experiment
            if (fits_file_base[-2:-1] != '_'):
                fits_file_base = fits_file_base + '_1'
            infile = fitsdir + '/' + fits_file_base + '_' + str(headnum) + '.IDI'
        else:
            infile = fitsdir + '/' + fits_file[i]

        if (not nfits):
            nfits = get_nfits(infile)
        if (not nfits):
            raise 'No fits files were found matching: ' + infile

        logger.info('Loading %d pieces of %s' % (nfits, infile))
        runfitld(infile=infile, outdata=uvdata, nfits=nfits, cltablemin=cltablemin, wtthreshhold=wtthreshhold, interferometer=interferometer)

        # sort the data if necessary
        if (not uvdata.header.sortord == 'TB'):
            logger.info('Data not in TB order. Need to MSORT the data')
            msortdata = AIPSUVData(headname, 'MSORT', uvdata.disk, 1)
            zap_old_data(msortdata)
            runmsort(indata=uvdata, outdata=msortdata)
            head_klass = 'MSORT'
            uvdata = msortdata

        # run indxr now we are sure the data are 'TB' sorted
        runindxr(uvdata, interferometer)

        headdata[i] = AIPSUVData(uvname + '_' + str(headnum),
                    head_klass, uvdata.disk, uvdata.seq)

    return uvdata, headdata

# Get the day of year from the Year, month, day
def get_day_of_year(year, month, day):
    """Get the day-of-year integer from the year/month/day

    year   I  YYYY
    month  I  MM starting from January == 1
    day    I  DD starting from 1 == 1
    """
    day_of_year_list = [0,31,59,90,120,151,181,212,243,273,304,334]
    doy = day_of_year_list[month-1] + day
    if(month>2):
        if((year&0x3)==0):
            if((year % 100 != 0) or (year % 400 == 0)):
                doy = doy+1
    return doy

# Get the day of year from the Year, month, day for the start of observations
def get_observation_year_month_day(uvdata):
    date_string = uvdata.header.date_obs
    date_list = date_string.split('-')
    year = int(date_list[0])
    month = int(date_list[1])
    day = int(date_list[2])
    return (year, month, day)

# Get the day of year from the Year, month, day for the start of observations
def get_observation_day_of_year(uvdata):
    (year, month, day) = get_observation_year_month_day(uvdata)
    return get_day_of_year(year, month, day)

# Get the number of days of observations
def get_number_days_observations(uvdata):
    """Find out how long the observations are.  Assumes the data has been sorted

    """
    nx_table = uvdata.table('AIPS NX', 0)
    num_rows = len(nx_table)
    # The NX table is FORTRANish
    # Check the first observation
    first_row = nx_table[0]
    num_days = int(nx_table[num_rows - 1]['time'] + 1)
    #start_time = first_row['TIME'][0] - 0.5 * first_row['TIME INTERVAL'][0]
    #if(start_time < 0.0):
    #    raise (RuntimeError, "Starting observation before starting date: Run FXTIM!")
    ## Now check the last scan
    #last_row = nx_table[num_rows-1]
    #end_time = last_row['TIME'][0] + 0.5 * last_row['TIME INTERVAL'][0]
    return num_days #int(math.ceil(end_time))

# Download TEC maps
def get_TEC(year,doy,TECU_model,tecdir):
    year=str(year)[2:4]
    if doy<10:
        doy='00'+str(doy)
    elif doy<100:
        doy='0'+str(doy)
    else:
        doy=str(doy)
    name=TECU_model+doy+'0.'+year+'i'
    if os.path.exists(name):
        logger.info('File already there.')
    else:
        #path='ftp://cddis.gsfc.nasa.gov/gps/products/ionex/20'+year+'/'+doy+'/'
        path='ftp://gdc.cddis.eosdis.nasa.gov/gnss/products/ionex/20'+year+'/'+doy+'/'
        #os.popen(r'wget -t 30 -O '+name+'.Z '+path+name+'.Z')
        os.system('curl --insecure -O --ftp-ssl '+path+name+'.Z')
        os.system('uncompress -f '+name+'.Z')
        os.system('mv %s %s' % (name, tecdir))
        

# Download EOP file
def get_eop(eop_path):
    if os.path.exists(eop_path+'usno_finals.erp'):
        #age = (time.time() - os.stat(eop_path+'usno_finals.erp')[8])/3600
        logger.info('usno_finals.erp exists, not downloaded.')
    else:
        os.system('curl --insecure -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp')   #wget ftp://cddis.gsfc.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno500_finals.erp http://gemini.gsfc.nasa.gov/solve_save ftp://ftp.lbo.us/pub/staff/wbrisken/EOP
        os.system('curl --insecure -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno500_finals.erp')
        os.system('mv usno500_finals.erp '+eop_path+'usno_finals2.erp')
        os.system('mv usno_finals.erp '+eop_path+'usno_finals.erp')

def check_table_version_exists(uvdata, table, version):
    """Check whether a table version exists, as 'AIPS CL', 'AIPS BP', etc.

    if version == 0, then check for ANY of that table, 
    """
    # If version is zero, then accept any of that table
    if(version == 0):
        for t in uvdata.tables:
            if(t[1] == table): break
        else: return 0
        return 1
    # Otherwise, check for the version number too
    for t in uvdata.tables:
        if((t[0] == version) and (t[1] == table)): break
    else: return 0
    return 1

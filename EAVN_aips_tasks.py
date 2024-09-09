# -*- coding: utf-8 -*-
"""
Created on Tue Sep 3 08:06:36 2024

@author: Yuzhu Cui, Baoqiang Lao
"""
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from AIPSTV import AIPSTV

import copy, optparse, os, sys
import re, string, pprint, math
import numpy as np
import logging

logger = logging.getLogger(__name__)

AIPSTask.version = '31DEC24' 

#task_versions = {}

#def versioned_task(taskname):
#    if taskname in task_versions:
#        return AIPSTask(taskname, version=task_versions[taskname])
#    else:
#        return AIPSTask(taskname)

FuncLog = sys.stdout
my_tv = AIPSTV()


def runfitld(outdata, infile, nfits=1):

    """Must set outdata, infile"""
    #assert (outdata != None, infile != None), "set indata, infile in runfitld"
    logger.info('Task FITLD  (release of 31DEC24) begins') 
    fitld = AIPSTask('fitld')
    fitld.outdata = outdata
    fitld.ncount = nfits
    fitld.douvcomp = -1
    fitld.doconcat = -1
    fitld.clint = 1./60.
    fitld.digicor = -1
    try:
        fitld.datain = infile
    except AttributeError:
        fitld.infile = infile
    logger.info('fits file = ' + fitld.datain)
    fitld()
    logger.info('Task FITLD appears to have ended successfully')


def runmsort(indata, outdata):
    '''Must set indata, outdata'''
    #assert(indata != None, outdata != None)
    logger.info('Task MSORT  (release of 31DEC24) begins')
    msort = AIPSTask('msort')
    msort.indata = indata
    msort.sort = 'TB'
    msort.outdata = outdata
    msort()   
    logger.info('Task MSORT appears to have ended successfully')

def runindxr(indata):
    logger.info('Task INDXR  (release of 31DEC24) begins')
    indxr = AIPSTask('indxr')
    indxr.indata = indata
    indxr.cparm[1:] = [0, 0, 1./60.,0]
    indxr()   
    logger.info('Task INDXR appears to have ended successfully')


def runuvcop(indata, outdata, freqid=0, bif=0, eif=0, bchan=0, echan=0,
        uvrange=[], timerang=[], antennas=[], baseline=[]):
    """Must set indata, outdata"""
    logger.info('Task UVCOP  (release of 31DEC24) begins')
    uvcop = AIPSTask('uvcop')
    uvcop.indata = indata
    uvcop.outdata = outdata
    uvcop.bif = bif
    uvcop.eif = eif
    uvcop.freqid = freqid
    uvcop.bchan = bchan
    uvcop.echan = echan
    uvcop.uvrange[1:] = uvrange
    uvcop.timerang[1:] = timerang
    uvcop.antennas[1:] = antennas
    uvcop.baseline[1:] = baseline
    uvcop()
    logger.info('Task UVCOP appears to have ended successfully')

def runlistr(indata, outprint, optype=''):#, inver=0, freqid=1,
        #docalib=0, inext='SN', stokes='HALF', sources=[], dparm=[]):
    """Must set indata, outprint"""
    #assert (indata != None, outprint != None)
    logger.info('Task LISTR  (release of 31DEC24) begins')
    listr = AIPSTask('listr')
    listr.indata = indata
    listr.outprint = outprint
    listr.optype = optype
    #listr.inext = inext
    #listr.inver = inver
    #listr.stokes = stokes
    #listr.freqid = freqid
    #listr.sources[1:] = sources
    #listr.bif = 0
    #listr.eif = 0
    #listr.bchan = 1
    #listr.echan = 0
    #listr.docalib = docalib
    #listr.gainuse = 0
    #listr.flagver = 0
    #listr.dparm[1:] = dparm[1:]
    #listr.dparm[3] = 6
    #listr.docrt = -1
    listr()
    logger.info('Task LISTR appears to have ended successfully')

def rundtsum(indata, outprint=''):
    '''Must set indata'''
    logger.info('Task DTSUM  (release of 31DEC24) begins')
    #assert (indata != None)
    dtsum = AIPSTask('dtsum')
    dtsum.indata = indata
    dtsum.aparm[1] = 1.0
    dtsum.outprint = outprint
    dtsum.docrt = -1
    dtsum()    
    logger.info('Task DTSUM appears to have ended successfully')

def runpossm(indata, aparm, solint=-1, freqid=1, sources=[],
        timerang=[0,0,0,1], codetype='A&P',
        nplots=9, stokes='HALF', antennas=[], docalib=0, gainuse=0,
        doband=0, bpver=0, flagver=1, dotv=-1):

    """Must set indata, aparm"""

    #assert (aparm != None, indata != None), "set indata, aparm in runpossm!"

    #print 'possm antennas=', antennas
    #raw_input('press return')

    possm = AIPSTask('possm')

    possm.stokes = stokes
    possm.indata = indata
    possm.solint = solint
    possm.freqid = freqid
    if (type(sources) == type('string')):
       possm.source[1] = sources
    else:
       possm.source[1:] = sources
    #possm.sources[1:] = sources
    possm.timerang[1:] = timerang
    possm.antennas[1:] = antennas
    possm.docalib = docalib
    possm.gainuse = gainuse
    possm.flagver = flagver
    possm.doband = doband
    possm.bpver = bpver
    possm.aparm[1:] = aparm[1:]
    possm.codetype = codetype
    possm.nplots = nplots
    possm.dotv = dotv
    possm()

def set_default_aparms(task):
    # simple function to set default aparms for various tasks. Returns default
    # values as a list, which can then be changed before calling the function
    # that runs the task. Note that the bparm elements match those in AIPS - do
    # not pass element zero to the tasks (e.g. use possm.aparm[1:] = aparm[1:])
    aparm = []
    aparm[1:] = [0 for i in range(10)]
    if task == 'possm':
        aparm[1] = 0
        aparm[2] = 1
        aparm[5] = -180
        aparm[6] = 180
        aparm[9] = 1

    if task == 'fring':
        aparm[1] = 3
        aparm[6] = 2
        # aparm[7] (the snr cutoff) is set in the call to fring
        #aparm[7] = 11
        aparm[9] = 1

    return aparm

def set_default_bparms(task):
    # cf aparms above
    bparm = []
    bparm[1:] = [0 for i in range(10)]
    if task == 'vplot':
        bparm[2] = -1
        bparm[3] = 1
        bparm[8] = -180
        bparm[9] = 180

    return bparm   

def runlwpla(indata, outfile, inver=None, plver=1):

    '''Must set indata, outfile'''
    #assert (indata != None, outfile != None)

    lwpla = AIPSTask('lwpla')
    lwpla.indata = indata
    lwpla.outfile = outfile
    lwpla.lpen = 1.
    lwpla.dparm[6] = 4
    lwpla.dparm[8] = 9
    lwpla.plver = plver
    if (inver == None):
        inver = indata.table_highver('AIPS PL')
    lwpla.inver = inver
    lwpla()


def runvplot(indata, bparm, antenna=[], in2data=None, freqid=1, bchan=1, 
        echan=4096, bif=1, eif=0, sources=[], docalib=1, gainuse=0, flagver=0,
        doband=-1, bpver=0, solint=0, nplots=4, ncomp=0, nmaps=0, stokes='',
        dotv=-1):

    """Must set indata, bparm"""
    #assert (indata != None, bparm != None  and antenna !=
    #        None), """bparms or antenna not set in runvplot!"""
    logger.info('Task VPLOT  (release of 31DEC24) begins')
    vplot = AIPSTask('vplot')
    vplot.indata = indata
    if is_aipsdata(in2data):
        vplot.in2data = in2data
    vplot.ncomp[1] = ncomp
    vplot.nmaps = nmaps
    vplot.sources[1:] = sources
    vplot.freqid = freqid 
    vplot.bchan = bchan
    vplot.echan = echan 
    vplot.bif = bif
    vplot.eif = eif 
    vplot.stokes = stokes
    vplot.antenna[1:] = antenna
    vplot.docalib = docalib
    vplot.gainuse = gainuse
    vplot.flagver = flagver
    vplot.doband = doband
    vplot.doebar = -1
    vplot.solint = solint 
    vplot.bparm[1:] = bparm[1:]
    vplot.nplots = nplots 
    vplot.ltype = 3 
    vplot.dotv = dotv
    vplot.grchan = 1
    error = 0
    try:
        vplot()
        error = 0
    except:
        error = 1
    logger.info('Task VPLOT appears to have ended successfully')
    return error   

def is_aipsdata(aipsdata):
    '''Check whether the passed object has the valid attributes for an AIPS
    data object'''

    got_attr = False
    if (hasattr(aipsdata, 'name') and hasattr(aipsdata, 'disk') and
            hasattr(aipsdata, 'seq') and hasattr(aipsdata, 'klass') ):
        got_attr = True

    return got_attr

def runaccor(indata):
    logger.info('Task ACCOR  (release of 31DEC24) begins')
    accor = AIPSTask('ACCOR')
    accor.indata = indata
    #accor.timer[1:] = [0]
    accor.solint = -1
    accor()
    logger.info('Task ACCOR appears to have ended successfully')

def runclcal(snver, indata, refant, gainver=0, gainuse=0, opcode='CALI',
        interpol='self', calsour=[], sources=[], samptype='',
        doblank=0, dobtween=0, inver=0):
    # function to run clcal. Don't allow default snver.
    """Must set snver, refant, indata"""

    #assert (snver != None, refant != None, indata != None), """missing snver,
    #        refant or indata in runclcal"""
    logger.info('Task CLCAL  (release of 31DEC24) begins')
    clcal = AIPSTask('clcal')
    clcal.indata = indata
    #clcal.subarray = 1
    clcal.calsour[1:] = calsour
    clcal.sources[1:] = sources
    #clcal.opcode = opcode
    clcal.interpol = interpol
    #clcal.intparm = 0.00001
    clcal.samptype = samptype
    #clcal.doblank = doblank
    #clcal.dobtween = dobtween
    clcal.refant = refant
    clcal.snver = snver
    clcal.inver = inver
    clcal.gainver = gainver
    clcal.gainuse = gainuse
    clcal()
    logger.info('Task CLCAL appears to have ended successfully')

def runsnplt(indata, inext='SN', invers=0, optype='PHAS', dotv=-1, antennas=[],
            bif=0, eif=0, sources=[]):

    '''Must set indata'''
    #assert(indata != None)
    logger.info('Task SNPLT  (release of 31DEC24) begins')
    if dotv > 0:
        my_tv.clear()

    snplt = AIPSTask('snplt')
    snplt.indata = indata
    snplt.inext = inext
    snplt.invers = invers
    snplt.optype = optype
    snplt.antennas[1:] = antennas
    snplt.sources[1:] = sources
    snplt.bif = bif
    snplt.eif = eif
    snplt.nplots = 10
    snplt.factor = 0.5
    #snplt.doebar = -1.
    snplt.cutoff = 1e-6
    snplt.dotv = dotv
    try:
        snplt()
        error = 0
    except:
        error = 1
    logger.info('Task SNPLT appears to have ended successfully')
    return error

def runapcal(indata, tyver=1, gcver=1, snver=1, freqid=1, opcode='grid', ant_id):
    """Must set indata"""
    #assert (indata != None)
    logger.info('Task APCAL  (release of 31DEC24) begins')
    apcal = AIPSTask('apcal')
    apcal.source = []
    apcal.antenna = []
    apcal.dofit[1:] = -1
    apcal.aparm[1] = 1.3
    apcal.dofit[ant_id] = 1
    apcal.indata = indata
    #apcal.freqid = freqid
    apcal.tyver = tyver
    apcal.gcver = gcver
    apcal.snver = snver
    apcal.opcode = opcode
    #apcal.prtlev = 1
    apcal()
    logger.info('Task APCAL appears to have ended successfully')

def runantab(antab_file, indata, tyver=1, gcver=1, offset=0., sparm=[]):

    """Must set indata, antab_file"""
    #assert (antab_file != None, indata != None), '''set antab_file, indata in
    #        runantab'''
    logger.info('Task ANTAB  (release of 31DEC24) begins')
    antab = AIPSTask('antab')
    try:
        antab.calin = antab_file
    except AttributeError:
        antab.infile = antab_file
    antab.indata = indata
    #antab.tyver = tyver
    #antab.gcver = gcver
    #antab.offset = 0.
    #antab.sparm[1:] = sparm
    antab()
    logger.info('Task ANTAB appears to have ended successfully')
   
def runbpass(refant, indata, calsour):

    """Must set refant, calsour, indata"""
    #assert (refant != None, calsour != None, indata != None), '''set
    #        refantlist in runbpass'''
    logger.info('Task BPASS  (release of 31DEC24) begins')
    bpass = AIPSTask('bpass')
    bpass.indata = indata
    bpass.calsour[1:] = calsour
    #bpass.freqid = 1
    #bpass.bif = 1
    #bpass.eif = 0
    #bpass.flagver= 1
    #bpass.subarray = 1
    bpass.docalib = 1
    bpass.gainuse = 3
    bpass.solint = 0
    bpass.refant = refant
    bpass.soltype = 'l1r'
    #bpass.bpver = 1
    #bpass.smooth[1] = 1
    bpass.bpassprm[1] = 1
    bpass.bpassprm[5] = 1
    bpass.bpassprm[9] = 1
    bpass.bpassprm[10] = 6
    bpass()   
    logger.info('Task BPASS appears to have ended successfully')

def runfring(indata, snver, solint, aparm, dparm, refant, freqid=1, 
        gainuse=0, calsour=[], snr):

    """Must set indata, snver, solint, refantlist, aparm, dparm"""
    #assert (indata != None, snver != None, solint != None, refantlist != None,
    #        aparm != None, dparm != None), """must set indata, snver, solint
    #        and refantlist for runfring"""
    logger.info('Task FRING  (release of 31DEC24) begins')
    fring = AIPSTask('fring')

    # calsour has limit of 30 sources. Split calsour into groups of
    # max. 30 and fring each group separately
    #bsource = -1
    #esource = -1
    #while bsource == -1 or esource < len(calsour)-1:
    #    bsource = esource+1
    #    esource = split_list(calsour, bsource, len(fring.calsour)-1)
    #    print >>FuncLog, 'Fringing sources ', str(bsource+1), ' to ', str(esource+1)

    fring.indata = indata
    fring.outdata = indata
    #fring.freqid = freqid
    fring.gainuse = gainuse
    fring.bpver = 1
    fring.doband = 1
    #fring.flagver = flagver
    fring.docalib = 1
    #fring.subarray = 1
    fring.refant = refant
    #fring.search[1:] = refantlist[1:]
    fring.snver = snver
    fring.solint = solint
    fring.calsour[1:] = calsour
    fring.aparm[1:] = aparm[1:]
    fring.aparm[7] = snr
    fring.dparm[1:] = dparm[1:]
    fring()
    logger.info('Task FRING appears to have ended successfully')

def runsplit(indata, sources=[], gainuse=0, doband=0, bpver=0, outseq=0,
        docalib=2):

    """Must set indata"""
    #assert (indata != None)
    split = AIPSTask('split')

    # 'sources' has limit of 30 sources. Split into groups of
    # max. 30 and split each group separately
    #bsource = -1
    #esource = -1
    #while bsource == -1 or esource < len(sources)-1:
    #    bsource = esource+1
    #    esource = split_list(sources, bsource, len(split.sources)-1)
    #    print >>FuncLog, 'Splitting sources ', str(bsource+1), ' to ', str(esource+1)

    split.indata = indata
    split.outclass = 'SPLIT'
    split.outdisk = indata.disk
    split.outseq = outseq
    #split.subarray = 1
    split.bchan = 17
    split.echan = 240
    split.sources[1:] = sources#[bsource:esource+1]
    split.docalib = docalib
    split.gainuse = gainuse
    split.doband = doband
    split.aparm[1] = 1
    #split.aparm[4] = 1
    split()


def runfittp(indata, outfile):

    """Must set indata, outfile"""
    #assert (outfile != None), 'set outfile in runfittp'

    fittp = AIPSTask ('fittp')
    fittp.indata = indata
    #fittp.doall = -1
    #fittp.intype = ''
    #fittp.outtape = 1
    try:
        fittp.dataout = outfile
    except AttributeError:
        fittp.outfile = outfile
    fittp()    

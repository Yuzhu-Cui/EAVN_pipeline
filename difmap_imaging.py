import subprocess
import time
import logging
from utils.uvfits import *
import os
import matplotlib.pyplot as plt
import optparse,re
from astropy.io import fits
from astropy import wcs
import aplpy 
import numpy as np
import matplotlib
from matplotlib.colors import SymLogNorm
from matplotlib.patches import Ellipse
from matplotlib.ticker import FormatStrFormatter
import pandas as pd

logfile = open("difmap.log", "w")

def pix2word(header,xy,core_xy):
    x,y = xy
    x = header['CDELT1']*3.6E6*(x-core_xy[0])
    y = header['CDELT2']*3.6E6*(y-core_xy[1])
    return x[0],y[0]

#find the lower and upper limits of X-axis and Y-axis in WCS 'mas' units for plotting (left, right, bottom, top):
def pixelwin2maswin(header,core_xy,win=()):
    if len(win)==4:
        x0, x1, y0, y1 = win
    else:
        x0, y0 = 0, 0
        x1, y1 = header['naxis1'], header['naxis2']
    x0, y0 = pix2word(header,(x0,y0),core_xy)
    x1, y1 = pix2word(header,(x1,y1),core_xy)
    win_mas = x0, x1, y0, y1
    return win_mas



def get_latest_log(log_dir):
    latest_log = None
    latest_mtime = 0
    for filename in os.listdir(log_dir):
        if filename.endswith('.log'):
            file_path = os.path.join(log_dir, filename)
            if os.path.getmtime(file_path) > latest_mtime:
                latest_mtime = os.path.getmtime(file_path)
                latest_log = filename
    return latest_log

def extract_content_between_markers(log_file, start_marker, end_marker):
    pattern = f'{re.escape(start_marker)}(.*?){re.escape(end_marker)}'
    with open(log_file, 'r') as file:
        content = file.read()
        matches = re.search(pattern, content, re.DOTALL)
        if matches:
            return matches.group(1)
        else:
            return None

def plot_cleanimage(fitsfile, xyrange):
    hdu= fits.open(fitsfile)[0]
    if len(hdu.data.shape) == 4:
       data = hdu.data[0,0,:,:]
    else:
       data = hdu.data[:,:]
    img = data 
    w1=wcs.WCS(hdu.header, naxis=2)
    header = hdu.header
    cra = header['CRVAL1'] #RA in deg
    cdec = header['CRVAL2']
    #print(img.shape)
    core_x, core_y = w1.wcs_world2pix([[cra,cdec]],0).transpose()
    core_xy = [core_x, core_y]
    win_mas = pixelwin2maswin(header,core_xy)
    
    rms =1.3* np.median(abs(img-np.median(img)))
    
    fig, ax = plt.subplots(figsize=(7,5))

    #norm=SymLogNorm(linthresh=10*rms)
    #print(win_mas)
    im = ax.imshow(img, vmin=-0.1,vmax=np.max(img), 
          extent=win_mas, origin='lower', cmap='jet', interpolation='none')

    levs_positive = 3*rms*np.array([1,np.sqrt(2),2,np.sqrt(2)*2,4,4*np.sqrt(2),8,8*np.sqrt(2),16,16*np.sqrt(2),32,32*np.sqrt(2),64,64*np.sqrt(2),128,128*np.sqrt(2),256,256*np.sqrt(2)])
    levs_negative = 3*rms*np.array([-1])   
    
    ax.contour(img,levs_positive,extent=win_mas,linestyles='solid',linewidths=1,colors='w')
    ax.contour(img,levs_negative,extent=win_mas,linestyles='dashed',linewidths=0.5,colors='r')
    ax.set_xlabel('Relative RA (mas)',fontsize=20)
    ax.set_ylabel('Relative Dec (mas)',fontsize=20,labelpad=0.01)
    ax.set_aspect('equal')

    #Set the colorbar
    cbar = fig.colorbar(im,aspect=25)
    cbar.ax.minorticks_off() ## IMPORTANT
    cbar.ax.tick_params('both',direction='in',right=True,top=True,which='both',labelsize = 20)
    #cbar.set_ticks([1,2,4,8,16,25])
    #cbar.set_ticklabels([1,2,4,8,16,25])
    #cbar.ax.set_position([0.88,0.08,0.2,0.9])
    cbar.set_label('Jy/beam',fontsize=20)

    #Specify tick label (value) font size
    ax.tick_params(axis = 'both', which = 'major', labelsize = 20)
    
    #set the format of the major tick label
    #It should be similar to C-programming format syntex.
    #Use 'd' for decimal, 'f' for float, and 's' for string valued ticks.
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))


    #Set the values ranges for x and y axes:
    k= 20
    xyrange = [int(xyrange[1]), int(xyrange[0]), int(xyrange[2]), int(xyrange[3])]
    ax.set_xlim(xyrange[0],xyrange[1])
    ax.set_ylim(xyrange[2],xyrange[3])
    
    #Add beam:
    b = (xyrange[0]-1.5*header['BMAJ']*3.6E6,xyrange[2]+1.5*header['BMAJ']*3.6E6,header['BMAJ']*3.6E6,header['BMIN']*3.6E6,header['BPA'])
    ellp = Ellipse(xy=(xyrange[1]+b[2]/2,xyrange[2]+b[3]*(1.2)),width=b[2],height=b[3],angle=90-b[4],ec='k',facecolor='grey')
    ax.add_artist(ellp)
    
    plt.tick_params(labelsize=20)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]

    fig.savefig(fitsfile.replace('.fits', '.png'), bbox_inches = 'tight')

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
    return control


def uvaver(inputfile):
    if os.path.exists(os.path.splitext(inputfile)[0] + "_uvaver.uvf"):
       print('%s already exists!' % (os.path.splitext(inputfile)[0] + "_uvaver.uvf"))
    else:
       difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
       difmap.stdin.write(("obs " + inputfile + "\n").encode())
       difmap.stdin.write(("select ll" + "\n").encode())
       difmap.stdin.write(("uvaver 30,true" + "\n").encode())
       difmap.stdin.write(("save " + os.path.splitext(inputfile)[0] + "_uvaver" + "\n").encode())
       difmap.stdin.write("exit\n\n".encode())
       #Wait for difmap exit
       difmap.communicate()
    uvaver_file = os.path.splitext(inputfile)[0] + "_uvaver.uvf"
    return uvaver_file


def check_data(uvfile):
    print('Start  checking data ... ...')
    uv = UVFits(uvfile)

    uv.u = uv.u_raw * uv.freq / 1e6 # multiply by frequency to get uv distance in wavelengths. 1e6 to convert to Megalambdas
    uv.v = uv.v_raw * uv.freq / 1e6
    uv.r = np.sqrt(uv.u**2 + uv.v**2)

    #print(len(uv.u))
    bl= uv.baselines
    suba = uv.subarrays
    ant1 = uv.ant1_inds
    ant2 = uv.ant2_inds
    ant1_n = []
    ant2_n = []
    for ind in range(uv.gcount):
        ant1_name, ant2_name = uv.get_ant_by_ind(ind)
        ant1_n.append(ant1_name)
        ant2_n.append(ant2_name)
    ant1_n = np.array(ant1_n)
    ant2_n = np.array(ant2_n)
    ant_n_all = list(set(list(set(ant1_n)) + list(set(ant2_n))))
    print('Antennas: ', ant_n_all)
    fig, ax = plt.subplots(2,1, figsize=(8,12))
    colors = ['b', 'k', 'r', 'g', 'm', 'y', 'c', 'pink', 'olive', 'tan']
    markers = ['o', '.', '^', '*', 'v', 's', 'D', '1', '2', '3']
    amp_means = []
    amp_medians = []
    phase_means = []
    phase_medians = []
    for m,ant in enumerate(ant_n_all):
        ant1_bl, = np.where(ant1_n == ant)
        ant2_bl, = np.where(ant2_n == ant)
        ant_bl = np.hstack((ant1_bl,ant2_bl))
        ant_bl_final = np.array(list(set(ant_bl)))
        #print(uv.r[ant_bl_final].shape, uv.amplitudes[ant_bl_final].shape)
        bl = uv.r[ant_bl_final]
        amp =  uv.amplitudes[ant_bl_final]
        amp_means.append(np.mean(amp))
        amp_medians.append(np.median(amp))
        phase = uv.phases[ant_bl_final]
        phase_means.append(np.mean(phase))
        phase_medians.append(np.median(phase))
        bl = bl.reshape((1, bl.size)).repeat(8, axis=0).reshape((bl.size * 8,))
        amp = amp.flatten()
        phase = phase.flatten() /np.pi * 180.0
        ax[0].scatter(bl, amp, marker=markers[m], label=ant, edgecolors=colors[m], color='') #, alpha=0.5)
        ax[1].scatter(bl, phase, marker=markers[m], label=ant, edgecolors=colors[m], color='') #, alpha=0.5) 
       #ax.scatter(uv.r[ant_bl_final], uv.amplitudes[ant_bl_final][:,0], '.', label=ant, color=colors[m])
    ax[0].set_ylabel('Correlated flux density [Jy]')
    ax[0].set_xlabel(r'Baseline length projection [M$\lambda$]')
    ax[1].set_ylabel('Phases [deg]')
    ax[1].set_xlabel(r'Baseline length projection [M$\lambda$]')    
    ax[0].legend()
    ax[1].legend()
    #print(os.path.basename(uvfile))
    plt.savefig('./data/radplot_check_%s.png' % os.path.basename(uvfile))
    plt.close(fig)
    amp_level = np.mean(np.array(amp_means))
    amp_diff = abs(np.array(amp_means) - amp_level)
    ant_id, = np.where(amp_diff==np.max(amp_diff))
    #print(ant_id)
    #print(amp_means, amp_medians)
    print('Abnormal antenna is: ', ant_n_all[ant_id[0]])
    print('Check done!')
    return ant_n_all[ant_id[0]]

def flag_data(inputfile, antn):
    if os.path.exists(os.path.splitext(inputfile)[0] + "_flag.uvf"):
       print('%s already exists!' % (os.path.splitext(inputfile)[0] + "_flag.uvf"))
    else:
       difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
       difmap.stdin.write(("obs " + inputfile + "\n").encode())
       difmap.stdin.write(("select ll" + "\n").encode())
       difmap.stdin.write(("flag %s" % antn + "\n").encode())
       difmap.stdin.write(("save " + os.path.splitext(inputfile)[0] + "_flag" + "\n").encode())
       difmap.stdin.write("exit\n\n".encode())
       #Wait for difmap exit
       difmap.communicate()
    flag_file = os.path.splitext(inputfile)[0] + "_flag.uvf"
    return flag_file


def difmap_imaging(vis_file, out_dir, flag_ant,
                   clean_gain,
                   clean_sigma,
                   map_size,
                   pixel_size,
                   clean_win_file, target_name, xyrange):
    print('Start imaging processing ... ...')    	
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
    print('read data')
    difmap.stdin.write(("obs " + vis_file + "\n").encode())
    print('select data')
    difmap.stdin.write(("select ll" + "\n").encode()) 
    # time average
    print('uv average')
    difmap.stdin.write(("uvave 30,true" + "\n").encode())
    # start a point source model
    difmap.stdin.write(("startmod" + "\n").encode())
    difmap.stdin.write(("save %s/kava2_startmod" % out_dir + "\n").encode())
    difmap.stdin.write(("mapcolor color" + "\n").encode())
    difmap.stdin.write(("mapsize %d,%f" % (map_size,pixel_size) + "\n").encode())
    print('flag bad antenna')
    difmap.stdin.write(("flag %s" % flag_ant + "\n").encode())
    difmap.stdin.write(("save %s/kava2_flagtianma" % out_dir + "\n").encode())
    # uniform weight
    difmap.stdin.write(("uvw 0,-1" + "\n").encode())
    difmap.stdin.write(("save %s/kava2_phase" % out_dir + "\n").encode())
    difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    print('clean ...')
    difmap.stdin.write(("clean 1000,%f" % clean_gain+ "\n").encode())
    #peakwin 1.5;
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
    difmap.stdin.write(("save  %s/be_gscaleclean" % out_dir + "\n").encode())
    print('gscale ...')
    difmap.stdin.write(("gscale" + "\n").encode())
    difmap.stdin.write(("clrmod true" + "\n").encode())
    #peakwin 1.5;
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
    difmap.stdin.write(("save %s/kava2_gscaleclean" % out_dir + "\n").encode())
    for i in [120, 60, 20, 10, 6, 3, 1, 0.5]:
        difmap.stdin.write(("selfcal true,true,%s" % i + "\n").encode())
        difmap.stdin.write(("clrmod true" + "\n").encode())
        difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
        difmap.stdin.write(("save %s/kava2_selfcal%sclean" % (out_dir, i) + "\n").encode())
    difmap.stdin.write(("print imstat(rms)" + "\n").encode())
    difmap.stdin.write(("unflag %s" % flag_ant + "\n").encode())
    difmap.stdin.write(("selfant ,true,1" + "\n").encode())
    difmap.stdin.write(("selfant tia,fal,1" + "\n").encode())
    difmap.stdin.write(("save %s/tian2_fixkava" % out_dir + "\n").encode())
    difmap.stdin.write(("uncal false,false,true" + "\n").encode())
    difmap.stdin.write(("save %s/tian2_recovertianma" % out_dir + "\n").encode())
    #difmap.stdin.write(("rmod tian2_recovertianma.mod" + "\n").encode())
    #clean-selfcal calibration:
    
    #difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    difmap.stdin.write(("clean 1000,%f" % clean_gain+ "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())

    difmap.stdin.write(("gscale" + "\n").encode())
    difmap.stdin.write(("clrmod true" + "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
    difmap.stdin.write(("save %s/tian2_flagisggscale" % out_dir + "\n").encode())


    for i in [120, 60, 30, 15, 6, 3, 1, 0.5]:
        difmap.stdin.write(("selfcal true,true,%s" % i + "\n").encode())
        difmap.stdin.write(("clrmod true" + "\n").encode())
        difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
        difmap.stdin.write(("save %s/tian2_flagisgselfcal%sclean" % (out_dir, i) + "\n").encode())
 

    #test open all antenna
    difmap.stdin.write(("selfant ,false,1" + "\n").encode())
    difmap.stdin.write(("save %s/tian2_openkava" % out_dir + "\n").encode())
    #difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    difmap.stdin.write(("clean 1000,%f" % clean_gain+ "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())

    difmap.stdin.write(("gscale" + "\n").encode())
    difmap.stdin.write(("clrmod true" + "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
    difmap.stdin.write(("save %s/tian2_opengscale" % out_dir + "\n").encode())


    for i in [120, 60, 30, 15, 6, 3, 1, 0.5]:
        difmap.stdin.write(("selfcal true,true,%s" % i + "\n").encode())
        difmap.stdin.write(("clrmod true" + "\n").encode())
        difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
        difmap.stdin.write(("save %s/tian2_openselfcal%sclean" % (out_dir, i) + "\n").encode())

    difmap.stdin.write(("print imstat(rms)" + "\n").encode())
    difmap.stdin.write(("delwin" + "\n").encode())
    difmap.stdin.write(("device /NULL" + "\n").encode())
    difmap.stdin.write(("mapl" + "\n").encode())
    if xyrange[0] != '':
       xyrange = '%s,%s,%s,%s' % (xyrange[0], xyrange[1], xyrange[2], xyrange[3])
       difmap.stdin.write(("xyrange %s" % xyrange + "\n").encode())

    difmap.stdin.write(("mapcolor color" + "\n").encode())
    difmap.stdin.write(("cmul=3*imstat(rms)" + "\n").encode())
    difmap.stdin.write(("levs=-1,1.00,1.41,2.00,2.83,4.00,5.65,7.99,11.30,15.98,22.60,31.95,45.18,63.88,90.33,127.73,180.61,255.38,361.11,510.61,722.00,1020.91,1443.57" + "\n").encode())
    #Save image
    difmap.stdin.write(("device %s/%s.ps/vcps" % (out_dir,target_name) + "\n").encode())
    difmap.stdin.write(("mapl cln" + "\n").encode())   
 
    difmap.stdin.write(("print \"MARKING_STRING\";print peak(flux);print imstat(rms);print cmul;print imstat(bmin);print imstat(bmaj);print imstat(bpa);print \"END_MARKING\"" + "\n").encode())
    difmap.stdin.write(("save %s/%s" % (out_dir, target_name) + "\n").encode())
    difmap.stdin.write(("exit\n\n").encode())
    difmap.communicate()
    print('imaging done!')

def difmap_imaging_vlba(vis_file, out_dir,
                   clean_gain,
                   clean_sigma,
                   map_size,
                   pixel_size,
                   clean_win_file, target_name):
    print('Start imaging processing ... ...')
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
    print('read data')
    difmap.stdin.write(("obs " + vis_file + "\n").encode())
    print('select data')
    difmap.stdin.write(("select ll" + "\n").encode())
    # time average
    print('uv average')
    difmap.stdin.write(("uvave 30,true" + "\n").encode())
    # start a point source model
    difmap.stdin.write(("startmod" + "\n").encode())
    difmap.stdin.write(("save %s/kava2_startmod" % out_dir + "\n").encode())
    difmap.stdin.write(("mapcolor color" + "\n").encode())
    difmap.stdin.write(("mapsize %d,%f" % (map_size,pixel_size) + "\n").encode())
    # uniform weight
    difmap.stdin.write(("uvw 0,-1" + "\n").encode())
    difmap.stdin.write(("save %s/kava2_phase" % out_dir + "\n").encode())
    difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    print('clean ...')
    difmap.stdin.write(("clean 1000,%f" % clean_gain+ "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())

    difmap.stdin.write(("gscale" + "\n").encode())
    difmap.stdin.write(("clrmod true" + "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
    difmap.stdin.write(("save %s/tian2_opengscale" % out_dir + "\n").encode())


    for i in [120, 60, 30, 15, 6, 3, 1, 0.5]:
        difmap.stdin.write(("selfcal true,true,%s" % i + "\n").encode())
        difmap.stdin.write(("clrmod true" + "\n").encode())
        difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal;  until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma, clean_sigma) + "\n").encode())
        difmap.stdin.write(("save %s/tian2_openselfcal%sclean" % (out_dir, i) + "\n").encode())

    difmap.stdin.write(("print imstat(rms)" + "\n").encode())
    difmap.stdin.write(("delwin" + "\n").encode())
    difmap.stdin.write(("device /NULL" + "\n").encode())
    difmap.stdin.write(("mapl" + "\n").encode())
    difmap.stdin.write(("mapcolor color" + "\n").encode())
    difmap.stdin.write(("cmul=3*imstat(rms)" + "\n").encode())
    difmap.stdin.write(("levs=-1,1.00,1.41,2.00,2.83,4.00,5.65,7.99,11.30,15.98,22.60,31.95,45.18,63.88,90.33,127.73,180.61,255.38,361.11,510.61,722.00,1020.91,1443.57" + "\n").encode())
    #Save image
    difmap.stdin.write(("device %s/%s.ps/vcps" % (out_dir,target_name) + "\n").encode())
    difmap.stdin.write(("mapl cln" + "\n").encode())
    difmap.stdin.write(("print \"MARKING_STRING\";print peak(flux);print imstat(rms);print cmul;print imstat(bmin);print imstat(bmaj);print imstat(bpa);print \"END_MARKING\"" + "\n").encode())
    difmap.stdin.write(("save %s/%s" % (out_dir, target_name) + "\n").encode())
    difmap.stdin.write(("exit\n\n").encode())
    difmap.communicate()
    print('imaging done!')


def difmap_modeling_vlba(visfile,
                          clean_sigma,
                          map_size,
                          pixel_size,
                          observation_length,
                          #model_sigma,
                          #phi,
                          #major_axis,
                          #minor_axis,
                          #model_type,
                          #model_iter,
                          #max_jet_component_number,
                          clean_win_file,
                          out_dir, target_name):
    print('Start imaging processing ... ...')
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
    difmap.stdin.write(("obs %s" % visfile + "\n").encode())
    difmap.stdin.write(("select ll" + "\n").encode())
    difmap.stdin.write(("mapsize %d,%f" % (map_size,pixel_size) + "\n").encode())
    difmap.stdin.write(("startmod \"\",1" + "\n").encode())     
    difmap.stdin.write(("uvw 2,0" + "\n").encode())
    #difmap.stdin.write(("peakwin 1.5" + "\n").encode())
    difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    difmap.stdin.write(("clean 100,0.05" + "\n").encode())
    difmap.stdin.write(("selfcal false,false,0" + "\n").encode())

    #if snr > clean_sigma more calibration steps will be done
    #clean-selfcal calibration:
    difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal; until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma,clean_sigma) + "\n").encode())
    difmap.stdin.write(("uvw 2,-1" + "\n").encode())
   
    #clean-selfcal calibration:
    difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal; until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma,clean_sigma) + "\n").encode())

    difmap.stdin.write(("uvw 0,-1" + "\n").encode())

    #clean-selfcal calibration:
    difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > %d) repeat; clean; selfcal; until(peak(flux,max)/imstat(rms) < %d) end if" % (clean_sigma,clean_sigma) + "\n").encode())
    difmap.stdin.write(("gscale" + "\n").encode())

    #here one should use tlpot to determine the true value of actual_time

    #observation time domain calibration:
    difmap.stdin.write(("float observation_length; observation_length = %f; rwin %s; float signal_to_noise_p; float signal_to_noise_a; signal_to_noise_p = peak(flux,max)/imstat(rms); \n \
    repeat; \n \
        selfcal true,true,observation_length; \n \
        signal_to_noise_a = peak(flux,max)/imstat(rms); \n \
        if(peak(flux,max)/imstat(rms) > %d) \n \
           if(signal_to_noise_a <= signal_to_noise_p) \n \
              clean; selfcal; \n \
              observation_length=observation_length/2; \n \
              signal_to_noise_a = signal_to_noise_p; \n \
           end if \n \
           if(signal_to_noise_a > signal_to_noise_p) \n \
              clrmod false,true; \n \
              observation_length=observation_length/2; \n \
              signal_to_noise_a = signal_to_noise_p; \n \
           end if \n \
        else \n \
           observation_length=observation_length/2; \n \
        end if \n \
    until(observation_length < 2)" % (observation_length, clean_win_file, clean_sigma) + "\n").encode())
    
    difmap.stdin.write(("selfcal true,true,0" + "\n").encode())
    difmap.stdin.write(("delwin" + "\n").encode())
    difmap.stdin.write(("clean 1000,0.01" + "\n").encode())	
    difmap.stdin.write(("device /NULL" + "\n").encode())
    difmap.stdin.write(("mapl" + "\n").encode())	
    difmap.stdin.write(("cmul=3*imstat(rms)" + "\n").encode())	

    #Save image
    difmap.stdin.write(("device %s/%s.ps/vcps" % (out_dir, target_name)  + "\n").encode())

    difmap.stdin.write(("mapl cln" + "\n").encode())
    difmap.stdin.write(("print \"MARKING_STRING\";print peak(flux);print imstat(rms);print cmul;print imstat(bmin);print imstat(bmaj);print imstat(bpa);print \"END_MARKING\"" + "\n").encode())
    difmap.stdin.write(("save %s/%s" % (out_dir, target_name) + "\n").encode())
    difmap.stdin.write(("exit\n\n").encode())
    difmap.communicate()
    print('imaging done!')


# main function
usage = 'usage: python3 %prog [options] imaging.inp'
parser = optparse.OptionParser(usage=usage, version='%prog difmap2.5q')
(options, args) = parser.parse_args()
if len(args) != 1:
    parser.error("incorrect number of arguments")



control = parse_inp(args[0])
outdir = control.get('outdir', [False])[0]
poitsourceuvfits  = control.get('poitsourceuvfits', [])[0]
visfile  = control.get('visfile', [])[0]
cleanwinfile  = control.get('cleanwinfile', [])[0]
targetsource = control.get('targetsource', [])[0]
cleangain = float(control.get('clean_gain', [])[0])
cleansigma = int(control['cleansigma'][0])
mapsize = int(control['mapsize'][0])
pixelsize = float(control.get('pixelsize', [0])[0])
interferometer = control.get('interferometer', [])[0]
xyrange  = control.get('xyrange', [])


if interferometer == 'EAVN':
   #1. strong point source processing
   #uv average 
   poitsource_uvfits = poitsourceuvfits
   uvaver_file = uvaver(poitsource_uvfits)
   
   #check antenna baseline data
   checked_flag_ant = check_data(uvaver_file)#poitsource_uvfits)
   
   #flag data test
   #antn = "TIA"
   #flag_file = flag_data(uvaver_file, antn)
   
   #print(flag_file)
   #
   #check_data(flag_file)
   
   # target source imaging
   vis_file = visfile
   flag_ant = checked_flag_ant 
   clean_win_file = cleanwinfile 
   target_source = targetsource
   difmap_imaging(vis_file=vis_file, out_dir=outdir, flag_ant=flag_ant,
                      clean_gain = cleangain,
                      clean_sigma=cleansigma,
                      map_size=mapsize,
                      pixel_size=pixelsize,
                      clean_win_file=cleanwinfile, target_name=targetsource, xyrange=xyrange)
   
   #os.system('ps2pdf ./data/%.ps')
   plot_cleanimage('%s/%s.fits' % (outdir, targetsource), xyrange)
   logfile.close()
   log_file = "difmap.log"
   START_MARKER = "MARKING_STRING"
   END_MARKER = "END_MARKING"
   extracted_content = extract_content_between_markers(log_file, START_MARKER, END_MARKER)
   if extracted_content:
       #imstat(rms);print cmul;print imstat(bmin);print imstat(bmaj);print imstat(bpa)
       outputs = extracted_content.split('\n')
       print('clean map info:')
       print('peak: %s Jy/beam' % outputs[1])
       print('rms: %s Jy/beam' % outputs[2])
       print('cmul: %s Jy/beam' % outputs[3])
       print('bmaj: %s mas' % outputs[5])
       print('bmin: %s mas' % outputs[4])
       print('bpa: %s deg' % outputs[6])
       data = {'file_name': [os.path.basename(vis_file)],
        'Source': [target_source],
        'Peak(Jy/beam)': [float(outputs[1])],
        'rms(Jy/beam)': [float(outputs[2])],
        'cmul(Jy/beam)': [float(outputs[3])],
        'bmaj(mas)': [float(outputs[5])],
        'bmin(mas)': [float(outputs[4])],
        'bpa(deg)': [float(outputs[6])],
        'bad_ant': [flag_ant]}
 
       df = pd.DataFrame(data)
       df.to_csv('%s.csv' % vis_file, index=False)
   else:
       print('No content found between markers.')

if interferometer == 'VLBA' or interferometer == 'EVN':
   # target source imaging
   vis_file = visfile
   clean_win_file = cleanwinfile
   target_source = targetsource

   difmap_imaging_vlba(vis_file=vis_file, out_dir=outdir,
                      clean_gain = cleangain,
                      clean_sigma=cleansigma,
                      map_size=mapsize,
                      pixel_size=pixelsize,
                      clean_win_file=cleanwinfile, target_name=targetsource)

   #difmap_modeling_vlba(visfile=vis_file,
   #                       clean_sigma=cleansigma,
   #                       map_size=mapsize,
   #                       pixel_size=pixelsize,
   #                       observation_length=120,
   #                       #model_sigma,
   #                       #phi,
   #                       #major_axis,
   #                       #minor_axis,
   #                       #model_type,
   #                       #model_iter,
   #                       #max_jet_component_number,
   #                       clean_win_file=cleanwinfile,
   #                       out_dir=outdir, target_name=target_source)

   #os.system('ps2pdf ./data/%.ps')
   plot_cleanimage('%s/%s.fits' % (outdir, targetsource), xyrange)
   logfile.close()
   log_file = "difmap.log"
   START_MARKER = "MARKING_STRING"
   END_MARKER = "END_MARKING"
   extracted_content = extract_content_between_markers(log_file, START_MARKER, END_MARKER)
   if extracted_content:
       #imstat(rms);print cmul;print imstat(bmin);print imstat(bmaj);print imstat(bpa)
       outputs = extracted_content.split('\n')
       print('clean map info:')
       print('peak: %s Jy/beam' % outputs[1])
       print('rms: %s Jy/beam' % outputs[2])
       print('cmul: %s Jy/beam' % outputs[3])
       print('bmaj: %s mas' % outputs[5])
       print('bmin: %s mas' % outputs[4])
       print('bpa: %s deg' % outputs[6])
       data = {'file_name': [os.path.basename(vis_file)],
        'Source': [target_source],
        'Peak(Jy/beam)': [float(outputs[1])],
        'rms(Jy/beam)': [float(outputs[2])],
        'cmul(Jy/beam)': [float(outputs[3])],
        'bmaj(mas)': [float(outputs[5])],
        'bmin(mas)': [float(outputs[4])],
        'bpa(deg)': [float(outputs[6])],
        }

       df = pd.DataFrame(data)
       df.to_csv('%s.csv' % vis_file, index=False)
   else:
       print('No content found between markers.')


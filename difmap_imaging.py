import subprocess
import time
import logging
from utils.uvfits import *
import os
import matplotlib.pyplot as plt

logfile = open("difmap.log", "w")

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
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    colors = ['b', 'k', 'r', 'g', 'm', 'y', 'c']
    markers = ['o', '.', '^', '*', 'v', 's', 'D']
    amp_means = []
    amp_medians = []
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
        ax.scatter(bl, amp[:,0], marker=markers[m], label=ant, edgecolors=colors[m], color='') #, alpha=0.5)
        #ax.scatter(uv.r[ant_bl_final], uv.amplitudes[ant_bl_final][:,0], '.', label=ant, color=colors[m])
    ax.set_ylabel('Correlated flux density [Jy]')
    ax.set_xlabel(r'Baseline length projection [M$\lambda$]')
    plt.legend()
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
                   clean_sigma,
                   map_size,
                   pixel_size,
                   clean_win_file, target_name):
    print('Start imaging processing ... ...')    	
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
    difmap.stdin.write(("obs " + vis_file + "\n").encode())
    difmap.stdin.write(("select ll" + "\n").encode()) 
    difmap.stdin.write(("uvave 30,true" + "\n").encode())
    difmap.stdin.write(("startmod" + "\n").encode())
    difmap.stdin.write(("save %s/kava2_startmod" % out_dir + "\n").encode())
    difmap.stdin.write(("mapcolor color" + "\n").encode())
    difmap.stdin.write(("mapsize %d,%f" % (map_size,pixel_size) + "\n").encode())
    difmap.stdin.write(("flag %s" % flag_ant + "\n").encode())
    difmap.stdin.write(("save %s/kava2_flagtianma" % out_dir + "\n").encode())
    #difmap.stdin.write(("uvw 0,-1" + "\n").encode())
    difmap.stdin.write(("save %s/kava2_phase" % out_dir + "\n").encode())
    difmap.stdin.write(("gscale" + "\n").encode())
    difmap.stdin.write(("save %s/kava2_gscaleclean" % out_dir + "\n").encode())
    for i in [120, 60, 20, 10, 6, 3, 1, 0.5]:
        difmap.stdin.write(("selfcal true,true,%s" % i + "\n").encode())
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
    difmap.stdin.write(("rwin %s" % clean_win_file + "\n").encode())
    difmap.stdin.write(("clean 1000,0.01" + "\n").encode())
    difmap.stdin.write(("if(peak(flux,max)/imstat(rms) > 3) repeat; peakwin 1.5; clean; selfcal;  until(peak(flux,max)/imstat(rms) < 3) end if" + "\n").encode())

    difmap.stdin.write(("gscale" + "\n").encode())
    difmap.stdin.write(("save %s/tian2_flagisggscale" % out_dir + "\n").encode())


    for i in [120, 60, 30, 15, 6, 3, 1, 0.5]:
        difmap.stdin.write(("selfcal true,true,%s" % i + "\n").encode())
        difmap.stdin.write(("save %s/tian2_flagisgselfcal%sclean" % (out_dir, i) + "\n").encode())
 

    difmap.stdin.write(("print imstat(rms)" + "\n").encode())
    difmap.stdin.write(("delwin" + "\n").encode())
    difmap.stdin.write(("device /NULL" + "\n").encode())
    difmap.stdin.write(("mapl" + "\n").encode())
    difmap.stdin.write(("mapcolor color" + "\n").encode())
    difmap.stdin.write(("cmul=3*imstat(rms)" + "\n").encode())
    #Save image
    difmap.stdin.write(("device %s/%s.ps/vps" % (out_dir,target_name) + "\n").encode())
    difmap.stdin.write(("mapl cln" + "\n").encode())    
    difmap.stdin.write(("print \"MARKING_STRING\";print peak(flux);print imstat(rms);print cmul;print imstat(bmin);print imstat(bmaj);print imstat(bpa);print \"END_MARKING\"" + "\n").encode())
    difmap.stdin.write(("save %s/%s" % (out_dir, target_name) + "\n").encode())
    difmap.stdin.write(("exit\n\n").encode())
    difmap.communicate()
    print('imaging done!')

#1. strong point source processing
#uv average 
poitsource_uvfits = './data/a17078a_1219+044.UVDATA.FITS'
uvaver_file = uvaver(poitsource_uvfits)

#check antenna baseline data
flag_ant = check_data(uvaver_file)#poitsource_uvfits)

#flag data test
#antn = "TIA"
#flag_file = flag_data(uvaver_file, antn)

#print(flag_file)
#
#check_data(flag_file)

# target source imaging
vis_file = './data/a17078a_M87.UVDATA.FITS'
flag_ant = 'TIA'
clean_win_file = './data/kava2_selfcal0.5clean.win'
target_source = 'M87' 
difmap_imaging(vis_file=vis_file, out_dir='./data', flag_ant=flag_ant,
                   clean_sigma=1,
                   map_size=2048,
                   pixel_size=0.02,
                   clean_win_file=clean_win_file, target_name=target_source)


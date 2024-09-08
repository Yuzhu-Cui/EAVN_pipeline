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
    uv = UVFits(uvfile)

    uv.u = uv.u_raw * uv.freq / 1e6 # multiply by frequency to get uv distance in wavelengths. 1e6 to convert to Megalambdas
    uv.v = uv.v_raw * uv.freq / 1e6
    uv.r = np.sqrt(uv.u**2 + uv.v**2)

    print(len(uv.u))
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
    print(ant_n_all)
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    colors = ['k', 'b', 'r', 'g', 'm', 'y', 'c']
    for m,ant in enumerate(ant_n_all):
        ant1_bl, = np.where(ant1_n == ant)
        ant2_bl, = np.where(ant2_n == ant)
        ant_bl = np.hstack((ant1_bl,ant2_bl))
        ant_bl_final = np.array(list(set(ant_bl)))
        print(uv.r[ant_bl_final].shape, uv.amplitudes[ant_bl_final].shape)
        bl = uv.r[ant_bl_final]
        amp =  uv.amplitudes[ant_bl_final]
        ax.scatter(bl, amp[:,0], marker='.', label=ant, color=colors[m]) #, alpha=0.5)
        #ax.scatter(uv.r[ant_bl_final], uv.amplitudes[ant_bl_final][:,0], '.', label=ant, color=colors[m])
    ax.set_ylabel('Correlated flux density [Jy]')
    ax.set_xlabel(r'Baseline length projection [M$\lambda$]')
    plt.legend()
    #print(os.path.basename(uvfile))
    plt.savefig('./data/radplot_check_%s.png' % os.path.basename(uvfile))
    plt.close(fig)

def flag_data(inputfile):
    if os.path.exists(os.path.splitext(inputfile)[0] + "_flag.uvf"):
       print('%s already exists!' % (os.path.splitext(inputfile)[0] + "_flag.uvf"))
    else:
       difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
       difmap.stdin.write(("obs " + inputfile + "\n").encode())
       difmap.stdin.write(("select ll" + "\n").encode())
       difmap.stdin.write(("flag TIA" + "\n").encode())
       difmap.stdin.write(("save " + os.path.splitext(inputfile)[0] + "_flag" + "\n").encode())
       difmap.stdin.write("exit\n\n".encode())
       #Wait for difmap exit
       difmap.communicate()
    flag_file = os.path.splitext(inputfile)[0] + "_flag.uvf"
    return flag_file


def difmap_imaging(vis_file,
                   output_name,
                   clean_sigma,
                   map_size,
                   pixel_size,
                   observation_length):

    float clean_sigma; clean_sigma  = {clean_sigma}
    float map_size; map_size = {map_size}
    float pixel_size; pixel_size = {pixel_size}
    float observation_length; observation_length = {observation_length}
    
    float signal_to_noise_p;
    float signal_to_noise_a;
    	
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
    difmap.stdin.write(("obs " + vis_file + "\n").encode())
    difmap.stdin.write(("select ll" + "\n").encode()) 
    difmap.stdin.write(("mapsize %d %d" % (map_size,pixel_size) + "\n").encode())
    difmap.stdin.write(("startmod %s,%d" % ("",1) + "\n").encode()) 
    difmap.stdin.write(("uvw 2,0" + "\n").encode())
    difmap.stdin.write(("peakwin 1.5" + "\n").encode())
    difmap.stdin.write(("clean 100,0.05" + "\n").encode())
    difmap.stdin.write(("selfcal false,false,0" + "\n").encode())
    
    #if snr > clean_sigma more calibration steps will be done
    
    #clean-selfcal calibration:
    #if(peak(flux,max)/imstat(rms) > clean_sigma)
    #	repeat;\
    #		peakwin 1.5; clean; selfcal
    #	until(peak(flux,max)/imstat(rms) < clean_sigma)
    #end if
    #
    #uvw 2,-1
    #
    #!clean-selfcal calibration:
    #if(peak(flux,max)/imstat(rms) > clean_sigma)
    #	repeat;\
    #		peakwin 1.5; clean; selfcal
    #	until(peak(flux,max)/imstat(rms) < clean_sigma)
    #end if
    #
    #uvw 0,-1
    #
    #!clean-selfcal calibration:
    #if(peak(flux,max)/imstat(rms) > clean_sigma)
    #	repeat;\
    #		peakwin 1.5; clean; selfcal
    #	until(peak(flux,max)/imstat(rms) < clean_sigma)
    #end if
    #
    #gscale
    #
    #!here one should use tlpot to determine the true value of actual_time
    #
    #!observation time domain calibration:
    #signal_to_noise_p = peak(flux,max)/imstat(rms)
    #repeat;\
    #	selfcal true,true,observation_length
    #	signal_to_noise_a = peak(flux,max)/imstat(rms)
    #	if(peak(flux,max)/imstat(rms) > clean_sigma)
    #		if(signal_to_noise_a <= signal_to_noise_p)
    #			peakwin 1.5; clean; selfcal
    #		observation_length=observation_length/2
    #			signal_to_noise_a = signal_to_noise_p
    #		end if
    #		if(signal_to_noise_a > signal_to_noise_p)
    #			clrmod false,true
    #		observation_length=observation_length/2
    #			signal_to_noise_a = signal_to_noise_p
    #		end if
    #	else
    #		observation_length=observation_length/2
    #	end if
    #until(observation_length < 2)
    #
    #selfcal true,true,0
    #
    #delwin
    #clean 1000,0.01
    #device /NULL
    #mapl
    #cmul=3*imstat(rms)
    #
    #!Save image
    #!device {obj}.ps/vcps
    #device {obj}.ps/vps
    #mapl clean, false 
    #
    #!get clean image and beam statistics
    #print "MARKING_STRING"
    #print peak(flux)
    #print imstat(rms)
    #print cmul
    #print imstat(bmin)
    #print imstat(bmaj)
    #print imstat(bpa)
    #print "END_MARKING"
    #
    #!save clean map into a fits file
    #save {obj}
    #
    #!quit from difmap
    #quit'''.format(obj_file=visibility_file,
    #			obj=output_name,
    #			clean_sigma=clean_sigma,
    #			map_size=map_size,
    #			pixel_size=pixel_size,
    #			observation_length=observation_length));
    #
    #fn.close();


#uv average
input_uvfits = './data/a17078a.cal.fits'
uvaver_file = uvaver(input_uvfits)

#check antenna baseline data
check_data(uvaver_file)

#flag data
flag_file = flag_data(uvaver_file)

#print(flag_file)
#
check_data(flag_file)




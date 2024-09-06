import subprocess
import time
import logging
from utils.uvfits import *
import os

logfile = open("difmap.log", "w")

def uvaver(inputfile):
    difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
    difmap.stdin.write(("obs " + inputfile + "\n").encode())
    difmap.stdin.write(("select ll" + "\n").encode())
    difmap.stdin.write(("uvaver 30,true" + "\n").encode())
    difmap.stdin.write(("save " + os.path.splitext(inputfile)[0] + "_uvaver" + "\n").encode())
    difmap.stdin.write("exit\n".encode())


def check_data(uvfile):
    uv = UVFits(uvfile)

    uv.u = uv.u_raw * uv.freq / 1e6 # multiply by frequency to get uv distance in wavelengths. 1e6 to convert to Megalambdas
    uv.v = uv.v_raw * uv.freq / 1e6
    uv.r = np.sqrt(uv.u**2 + uv.v**2)

    print(len(uv.v))
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

    #fig, ax = plt.subplots(1,1, figsize=(8,6))
    #ax.plot(uv.r, uv.amplitudes, 'o', label='radplot')
    #ax.set_ylabel('Correlated flux density [Jy]')
    #ax.set_xlabel(r'Baseline length projection [M$\lambda$]')
    #plt.savefig('lisakov.png')
    #plt.close(fig) 

#uv average
inputfile = './data/a17078a.cal.fits'
uvaver(inputfile)

check_data(inputfile)

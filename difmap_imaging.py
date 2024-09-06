import subprocess
import time
import logging
from utils.uvfits import *

logfile = open("difmap.log", "w")

def check_data(uvfile):
    uv = UVFits(uvfile)
    
    uv.u = uv.u_raw * uv.freq / 1e6 # multiply by frequency to get uv distance in wavelengths. 1e6 to convert to Megalambdas 
    uv.v = uv.v_raw * uv.freq / 1e6
    uv.r = np.sqrt(uv.u**2 + uv.v**2)
    
    print(uv.baselines)
    print(uv.subarrays, uv.ant1_inds, uv.ant2_inds)

    #fig, ax = plt.subplots(1,1, figsize=(8,6))
    #ax.plot(uv.r, uv.amplitudes, 'o', label='radplot')
    #ax.set_ylabel('Correlated flux density [Jy]')
    #ax.set_xlabel(r'Baseline length projection [M$\lambda$]')
    #plt.savefig('lisakov.png')
    #plt.close(fig) 
#test    
#inputfile = "./data/a17078a.cal.fits"
#difmap = subprocess.Popen("difmap", stdin=subprocess.PIPE, stdout=logfile)
#difmap.stdin.write(("obs " + inputfile + "\n").encode())
#difmap.stdin.write("exit\n\n".encode())

check_data('./data/a17078a.cal.fits')

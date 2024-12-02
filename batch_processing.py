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
import Wizardry.AIPSData
import numpy as np
import pandas as pd

M87_name_list = ['M87', '3C274', '1228+126']


def delete_temp():
    if os.path.exists('tmp.txt'):
        os.remove('tmp.txt')
    if os.path.exists('tmp1.txt'):
        os.remove('tmp1.txt')
    if os.path.exists('tmp2.txt'):
        os.remove('tmp2.txt')

def get_baselineLen_sort(uvdata):
    an_tab = uvdata.table('AN', 1)
    baselineLen = {}
    ant_names = []
    for m, row in enumerate(an_tab):
        baselineLen["%d" % m] = []
        for row2 in an_tab:
            xsep = row.stabxyz[0] - row2.stabxyz[0]
            ysep = row.stabxyz[1] - row2.stabxyz[1]
            zsep = row.stabxyz[2] - row2.stabxyz[2]
            bl = np.sqrt((xsep * xsep) + (ysep * ysep) + (zsep * zsep))
            baselineLen["%d" % m].append(bl)
        baselineLen["%d" % m] = np.sum(baselineLen["%d" % m])
        #ant_names.append(row.anname)
    baselineLen_pd = pd.DataFrame.from_dict(
        baselineLen, orient='index')
    baselineLen_pd.reset_index(drop=False, inplace=True)
    baselineLen_pd.rename(
        columns={'index': 'ant_num', 0: 'baselineLen'}, inplace=True)
    baselineLen_pd['ant_num'] = baselineLen_pd['ant_num'].astype(
        int) + 1
    baselineLen_pd.sort_values(
        by=['baselineLen'], ascending=True, inplace=True)
    baselineLen_sort = list(baselineLen_pd['ant_num'])
    print(list(baselineLen_pd['baselineLen']))
    #print(ant_names)
    return baselineLen_sort

def get_refant_id(uvdata):
    centralAntennas = get_baselineLen_sort(uvdata)
    # Removing temp files:
    delete_temp()
    # Running DTSUM task to get visibility data summary for each antenna:
    dtsum = AIPSTask('DTSUM')
    dtsum.default()
    dtsum.indata = uvdata
    dtsum.docrt = -3  # To suppress header information while writing.
    dtsum.outprint = 'tmp.txt'
    dtsum.go()
    # Modifying the visibility data summary text file:
    with open('tmp.txt', 'r') as infile:
        with open('tmp1.txt', 'w') as outfile:
            prtline = False
            for line in infile:
                if '---' in line:
                    prtline = True
                if prtline:
                    outfile.write(line)
    lines = open('tmp1.txt', 'r').readlines()
    # open('tmp_test2.txt', 'w').writelines(lines[1:-1])
    with open('tmp2.txt', 'w') as outfile:
        for line in lines[1:-1]:
            # line=line.split('|')[1]
            line = line.replace(' | ', ' ', 1)
            outfile.write(line)
    # Read ASCII table file as a dataframe:
    df1 = pd.read_csv('tmp2.txt', sep='\s+', na_filter=False, header=None)
    # Making a column as the dataframe row indices:
    df2 = df1.set_index(keys=df1[0], drop=True, append=False, inplace=False)
    df2.drop(df2.columns[[0]], axis=1, inplace=True)
    # Modifying the dataframe element-wise:
    for i in range(len(df2)):
        for j in range(len(df2)):
            if df2.iloc[i, j] == 0:
                df2.iloc[i, j] = df2.iloc[j, i]
    # Averaging each antenna column's visibilities:
    antMeanVisibilities = df2.mean(axis=0).values
    # Taking care of any completely flagged antenna:
    df3 = pd.DataFrame()
    df3[0] = df1[0]
    df3[1] = antMeanVisibilities
    # Number of available antennas:
    N_ant = len(df3)
    # Removing temp files:
    delete_temp()
    # Find the 3 central antennas with maximum number of visibilities out of first 5 central antennas:
    x = []
    y = []
    for i in range(5):
        j = centralAntennas[i]
        k = 0.0
        for a in range(len(df3)):
            if df3.iloc[a, 0] == j:
                k = df3.iloc[a, 1]
        if k == 0.0:
            j = np.nan
        x.append(j)
        y.append(k)
    x = np.array(x)
    x1 = x[~np.isnan(x)]
    refant_id = x1[:3].tolist()
    return refant_id

data_dir = './data'

interferometers = ['EAVN', 'VLBA']

AIPS.userno = 1001
# https://astrogeo.org/sol/rfc/rfc_2024b/rfc_2024b_cat.txt
rfc_csv = pd.read_csv('./rfc_2024b_catalogue/rfc_2024b_cat.csv')
IVS_names = rfc_csv['IVS_name'].values
IAU_names = rfc_csv['IAU_name'].values

for interferometer in interferometers:
    fs = os.listdir('%s/%s/raw' % (data_dir, interferometer))
    for f in fs:
        if f.endswith('.fits') or f.endswith('.uvfits') or f.endswith('.idifits'):
           fits_file = f
           if interferometer == 'EAVN':
               antab_file = '%s/%s/raw/%s.antab' % (data_dir, interferometer, os.path.splitext(fits_file)[0]) 
               if not os.path.exists(antab_file):
                  download_data_dir = '%s/%s/%s' %  (data_dir, interferometer, os.path.splitext(fits_file)[0])
                  os.system("cat %s/*VERA* > %s" % (download_data_dir, antab_file))
                  os.system("cat %s/*.ANTAB.1 >> %s" % (download_data_dir, antab_file))
                  os.system("sed -i \"s/INDEX = 'L1:16'/INDEX = 'L1:8'/g\" %s" % antab_file)
                  os.system("sed -i \"s/INDEX = 'L1','L2:16'/INDEX = 'L1:8','L5:8'/g\" %s" % antab_file)
                  os.system("sed -i \"s/INDEX='L1','L2','L3','L4'/INDEX='L1:4','L5:8','L9:12','L13:16'/g\" %s" % antab_file)
                  os.system("sed -i \"s/INDEX='L1:7:2','R2:8:2'/INDEX ='L1:4','L5:8'/g\" %s" % antab_file)
                  os.system("sed -i \"s/'R1','R2','R3','R4'/'L5','L6','L7','L8'/g\" %s" % antab_file)
                  os.system("sed -i \"s/INDEX='L1','L2:16'/INDEX='L1:4','X'/g\" %s" % antab_file)
                  os.system("sed -i \"s/'R1','L2','L3','L4'/INDEX='X','L1:4','X','X'/g\" %s" % antab_file)
                  os.system("sed -i \"s/INDEX='L1:1'/INDEX='L1:16'/g\" %s" % antab_file)
                  os.system("sed -i \"s/T6/TIA/g\" %s" % antab_file)
                  
           #print(fits_file) 
           experiment = os.path.splitext(fits_file)[0]
           print(experiment)
           disk = 1
           uvname = experiment.upper()
           uvdata = AIPSUVData(uvname, 'UVDATA', disk, 1)
           msortdata = AIPSUVData(uvdata.name, 'MSORT', uvdata.disk, 1)
           nfits = 1 
           fitsdir = '%s/%s/raw' % (data_dir, interferometer)
           cltablemin = 1/60.
           heads = [1]
           wtthreshhold = 0
           uvdata, headdata = load_data(uvdata, msortdata, experiment, fitsdir,
                   [fits_file], heads, nfits, cltablemin, wtthreshhold, interferometer)


           # Always do this bit to make sure names are right at end of block. The order of
           # the if statements here is significant.
           if msortdata.exists():
              uvdata = msortdata
           
           #antenna_names = get_antenna_names(uv=uvdata, version = 1)
           #print(antenna_names)
           print(uvdata.antennas)
           for source in uvdata.sources:
               if source in M87_name_list:
                  target = source               
           print(uvdata.sources) 
           print(uvdata.header['telescop'])
           #ant_table = uvdata.table('AN', 1)
           bl_sort = get_baselineLen_sort(uvdata)
           print(bl_sort)
           refant_id = get_refant_id(uvdata)
           flux_s = []
           flux_c = []
           flux_x = []
           flux_u = []
           flux_k = []
           category = []
           for source in uvdata.sources:
               if source in IVS_names:
                  rfc_csv_new = rfc_csv[rfc_csv['IVS_name']==source]
                  flux_s.append(float(rfc_csv_new['S-band_total_flux'].values[0]))
                  flux_c.append(float(rfc_csv_new['C-band_total_flux'].values[0]))   
                  flux_x.append(float(rfc_csv_new['X-band_total_flux'].values[0]))
                  flux_u.append(float(rfc_csv_new['U-band_total_flux'].values[0]))
                  flux_k.append(float(rfc_csv_new['K-band_total_flux'].values[0]))
                  category.append(rfc_csv_new['Category'].values[0])
               elif source in IAU_names:
                  rfc_csv_new = rfc_csv[rfc_csv['IAU_name']==source]
                  flux_s.append(float(rfc_csv_new['S-band_total_flux'].values[0]))
                  flux_c.append(float(rfc_csv_new['C-band_total_flux'].values[0]))
                  flux_x.append(float(rfc_csv_new['X-band_total_flux'].values[0]))
                  flux_u.append(float(rfc_csv_new['U-band_total_flux'].values[0]))
                  flux_k.append(float(rfc_csv_new['K-band_total_flux'].values[0]))
                  category.append(rfc_csv_new['Category'].values[0])
               else:
                  flux_s.append(np.nan)
                  flux_c.append(np.nan)
                  flux_x.append(np.nan)
                  flux_u.append(np.nan)
                  flux_k.append(np.nan)
                  category.append('')
           print(flux_s)
           print(flux_c)
           print(flux_x)
           print(flux_u)
           print(flux_k)  
           print(category)
           flux_s = np.array(flux_s)
           flux_c = np.array(flux_c)
           flux_x = np.array(flux_x)
           flux_u = np.array(flux_u)
           flux_k = np.array(flux_k)
           flux_all = np.vstack((flux_s, flux_c, flux_x, flux_u, flux_k))
           print(flux_all)
           flux_all[flux_all<=0.0] = 0.0
           print(flux_all)
           flux_max = np.nanmax(flux_all, axis=0)
           print(flux_max)
           cali_id = np.where(flux_max == np.nanmax(flux_max))[0]
           print(cali_id)
           bpass_calibrator = uvdata.sources[cali_id[0]]
           print(bpass_calibrator)
           # https://astrogeo.org/sol/rfc/rfc_2024b/rfc_2024b_cat.txt
           #
           with open('template.inp', 'r') as file:
               lines = file.readlines()
 
           new_filename = '%s/%s/%s/%s.inp' % (data_dir, interferometer, experiment, experiment)
           param_names = ['experiment', 'userno', 'refant', 
                          'plotref', 'bpass', 'phaseref',
                           'target', 'fitsdir', 'indir',
                           'outdir', 'fits_file', 'interferometer',
                           'eop_path', 'tecdir', 'fit_ant']
           fitsdir = '%s/%s/raw' % (data_dir, interferometer)
           outdir = '%s/%s/%s' % (data_dir, interferometer, experiment)
           if not os.path.exists(outdir):
              os.makedirs(outdir)
           tecdir = '%s/%s/%s/' % (data_dir, interferometer, experiment)
           if 'TIA' in uvdata.antennas:
              fit_ant = 'TIA'
           else:
              fit_ant = ''
           new_values = [experiment, '1001', uvdata.antennas[refant_id[0]-1], 
                         uvdata.antennas[refant_id[0]-1], bpass_calibrator, 
                         bpass_calibrator, target, fitsdir, fitsdir,
                         outdir, fits_file, interferometer, tecdir,
                         tecdir,fit_ant 
                         ] 
           with open(new_filename, 'w') as file:
               for line in lines:
                   for index, param_name in enumerate(param_names): 
                       if line.startswith(param_name):
                           items = line.split('=')
                           items[-1] = str(new_values[index])
                           line = ' = '.join(items)
                           line = '%s\n' % line 
                   file.write(line)
           os.system('ParselTongue EAVN.py %s' % new_filename)
           if interferometer == 'EAVN':
              with open('imaging.inp', 'r') as file:
                  lines_img = file.readlines()
           else:
              with open('imaging_vlba.inp', 'r') as file:
                  lines_img = file.readlines()

           new_filename_img = '%s/%s/%s/imaging_%s.inp' % (data_dir, interferometer, experiment, experiment)
           param_names_img = ['poitsourceuvfits', 'visfile', 'targetsource', 'outdir']
           new_values_img = ['%s/%s/%s/%s_1219+044.UVDATA.FITS' % (data_dir, interferometer, experiment,experiment),
                         '%s/%s/%s/%s_%s.UVDATA.FITS' % (data_dir, interferometer, experiment, experiment, target),
                         '%s_%s' % (experiment, target), outdir]

           with open(new_filename_img, 'w') as file:
               for line in lines_img:
                   for index, param_name in enumerate(param_names_img):
                       if line.startswith(param_name):
                           items = line.split('=')
                           items[-1] = str(new_values_img[index])
                           line = ' = '.join(items)
                           line = '%s\n' % line
                   file.write(line)
 
           os.system('python3 difmap_imaging.py %s' % new_filename_img)




            

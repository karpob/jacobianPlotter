#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import argparse, os, sys, h5py,glob
import numpy as np
from lib.graphics.profile import plotContour, plotLines,plotContourLabelIdx
from lib.graphics.linePlots import basicLine
    
def readKmatrixFile(filename):
    
    channels = []
    levels = []
    pressures = []
    jacobians = []
    with open(filename) as f:
        comment = f.readline()
        lat,lon = f.readline().split()
        comment = f.readline()
        lines = f.readlines()
        for l in lines:
            c,lv,p,jac = l.split()
            if(int(c) not in channels): channels.append(int(c))
            if(int(lv) not in levels): levels.append(int(lv))
            if(float(p) not in pressures): pressures.append(float(p))
            jacobians.append(float(jac))
    jacobians = np.asarray(jacobians).reshape((len(channels),len(pressures)))
    return np.asarray(channels), float(lat), float(lon), np.asarray(levels), np.asarray(pressures), jacobians
       
def readGsiCrtmKmatrix(path):
    files = glob.glob( os.path.join(path,'kmatrix*.txt') )
    d = {}
    for f in files:
        print("Reading: {}".format(f))
        if('ATM' in os.path.basename(f)):
            c, lat, lon, lev, pres, tJac = readKmatrixFile(f)
        elif('WV' in os.path.basename(f)):
            c, lat, lon, lev, pres, qJac = readKmatrixFile(f)
    d['T'] = tJac
    d['Q'] = qJac
    d['lat'] = lat
    d['lon'] = lon
    d['pres'] = pres
    d['channels'] = c       
    return d


if __name__ == "__main__":

    pathToThisScript = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser( description = 'plot jacobians from crtm kmatrix call in gsi')
    parser.add_argument('--path', help = 'path to ascii', required = True, dest = 'path')
    parser.add_argument('--outpath', help = 'path to outut', required = False, default=pathToThisScript,  dest = 'outpath')
    a = parser.parse_args()

    d = readGsiCrtmKmatrix(a.path)
    chans = d['channels'] 
    h5 = h5py.File(os.path.join('etc','cris-fsr_wavenumbers.h5'),'r')
    subsetChannelNumbers = np.asarray(h5['idxBufrSubset']).astype('int')
    wavenumbers2211 = np.asarray(h5['wavenumbers'])
    wavenumbers431 = wavenumbers2211[subsetChannelNumbers-1]
    wavenumbersToShow = wavenumbers431[np.asarray(chans)-1]
    idxLw, = np.where( wavenumbersToShow < 1095.1)
    idxSw,  = np.where( wavenumbersToShow > 2154.9)
    sensitivities = ['T','Q']
    wv = wavenumbersToShow
    profileNames = ['fun']
    for i,n in enumerate(profileNames): 
        for s in sensitivities:
            sValGsi = d[s]
            pressure = d['pres']
            maxS = sValGsi.max().max()
            minS = sValGsi.min().min()
            symMaxS = max(abs(minS),abs(maxS))
            symMinS = -1.0*symMaxS


            maxS2 = sValGsi[idxSw].max().max()
            minS2 = sValGsi[idxSw].min().min()
            symMaxS2 = max(abs(minS2),abs(maxS2))
            symMinS2 = -1.0*symMaxS2

            maxS3 = sValGsi[idxLw].max().max()
            minS3 = sValGsi[idxLw].min().min()
            symMaxS3 = max(abs(minS3),abs(maxS3))
            symMinS3 = -1.0*symMaxS3
 

           
            plotContour(wv, pressure, sValGsi,\
                        'Wavenumber [cm$^{-1}$]','Pressure [hPa]','Jacobian [dTb/dx]',\
                        'CRTM {} Jacobian'.format(s),\
                        os.path.join(a.outpath,'{}_crtm_lwsw.png'.format( s.lower() )),\
                        zlim = [symMinS, symMaxS], figureResolution=300 ) 
            wvTrunc = []
            for w in wv[idxSw]:
                wvTrunc.append('{0:.3f}'.format(w))
            plotContourLabelIdx(wvTrunc, pressure, sValGsi[idxSw],\
                        'Wavenumber [cm$^{-1}$]','Pressure [hPa]','Jacobian [dTb/dx]',\
                        'CRTM {} Jacobian SWIR'.format(s),\
                        os.path.join(a.outpath,'{}_crtm_sw.png'.format( s.lower() )),\
                        zlim = [symMinS2,symMaxS2],\
                        figureResolution=300 ) 
            wvTrunc = []
            for w in wv[idxLw]:
                wvTrunc.append('{0:.3f}'.format(w))
            plotContourLabelIdx(wvTrunc, pressure, sValGsi[idxLw],\
                        'Wavenumber [cm$^{-1}$]','Pressure [hPa]','Jacobian [dTb/dx]',\
                        'CRTM {} Jacobian LWIR'.format(s),\
                        os.path.join(a.outpath,'{}_crtm_lw.png'.format( s.lower() )),\
                        zlim = [symMinS3,symMaxS3],\
                        figureResolution=300 ) 
            for ii in range(idxLw.shape[0]):
                idxPair = np.array([idxLw[ii],idxSw[ii]]) 
                plotLines ( sValGsi[idxPair], pressure, "Jacobian [dTb/dx]", \
                       "Pressure [hPa]",(wv[idxPair]).tolist(),\
                       'Jacobians', os.path.join(a.outpath,'{}_crtm_sub_{}_{}_lines.png'.format( s.lower(),wv[idxLw[ii]],wv[idxSw[ii]] )),\
                        ylim =[100,1000],figureResolution = 300 )

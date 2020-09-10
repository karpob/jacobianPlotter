#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import argparse, os, sys, h5py,glob
import numpy as np
from lib.graphics.profile import plotContour, plotLines,plotContourLabelIdx
from lib.graphics.linePlots import basicLine
    
def readGsiTlmFile(filename):
    tLev = []
    tPres = []
    tJac = []
    qLev = []
    qPres = []
    qJac = []

    with open(filename) as f:
        comment = f.readline()
        channel,lat,lon = f.readline().split()
        comment = f.readline()
        buf = f.readline()
        while('model' not in buf):
            l,p,jac = buf.split()
            tLev.append(int(l))
            tPres.append(float(p))
            tJac.append(float(jac))
            buf = f.readline()
        buf = f.readline()
        while(len(buf) != 0):
            l,p,jac = buf.split()
            qLev.append(int(l))
            qPres.append(float(p))
            qJac.append(float(jac))
            buf = f.readline()
           
    return int(channel), float(lat), float(lon), np.asarray(tLev), np.asarray(tPres), np.asarray(tJac),np.asarray(qLev), np.asarray(qPres), np.asarray(qJac)
          
def readGsiTLMs(path):
    files = glob.glob( os.path.join(path,'jacobian*ch*.txt') )
    d = {}
    d['Tv'] = {}
    d['Moisture'] = {}
    for f in files:
        print("Reading: {}".format(f))
        c, lat, lon, lev, pres, tJac ,_ ,_ ,qJac = readGsiTlmFile(f)
        d['Tv'][c] = tJac
        d['Moisture'][c] = qJac
    d['lat'] = lat
    d['lon'] = lon
    d['pres'] = pres       
    return d
       


if __name__ == "__main__":
    pathToThisScript = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser( description = 'plot jacobians from setuprad in gsi')
    parser.add_argument('--path', help = 'path to ascii', required = True, dest = 'path')
    parser.add_argument('--outpath', help = 'path to outut', required = False, default=pathToThisScript,  dest = 'outpath')
    a = parser.parse_args()



    d = readGsiTLMs(a.path)
    chans = list(d['Tv'].keys())
    chans.sort()
    Tv = np.zeros([len(chans),len(d['pres'])])
    Moisture = np.zeros([len(chans),len(d['pres'])])
    i = 0
    for i,c in enumerate(chans):
        Tv[i,:] = d['Tv'][c][:]
        Moisture[i,:] = d['Moisture'][c][:]
    h5 = h5py.File(os.path.join('etc','cris-fsr_wavenumbers.h5'),'r')
    subsetChannelNumbers = np.asarray(h5['idxBufrSubset']).astype('int')
    wavenumbers2211 = np.asarray(h5['wavenumbers'])
    wavenumbers431 = wavenumbers2211[subsetChannelNumbers-1]
    wavenumbersToShow = wavenumbers431[np.asarray(chans)-1]
    idxLw, = np.where( wavenumbersToShow < 1000.0)
    idxSw,  = np.where( wavenumbersToShow > 1000.0)
    d['Tv'] = Tv
    d['Moisture'] = Moisture
    sensitivities = ['Tv','Moisture']
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

            maxS3 = sValGsi[idxSw].max().max()
            minS3 = sValGsi[idxSw].min().min()
            symMaxS3 = max(abs(minS3),abs(maxS3))
            symMinS3 = -1.0*symMaxS3

            plotContour(wv, pressure, sValGsi,\
                        'Wavenumber [cm$^{-1}$]','Pressure [hPa]','Jacobian [dTb/dx]',\
                        'GSI {} Jacobian'.format(s),\
                        os.path.join(a.outpath,'{}_gsi_lwsw.png'.format( s.lower() )),\
                        zlim = [symMinS, symMaxS], figureResolution=300 ) 
            wvTrunc = []
            for w in wv[idxSw]:
                wvTrunc.append('{0:.3f}'.format(w))
            plotContourLabelIdx(wvTrunc, pressure, sValGsi[idxSw],\
                        'Wavenumber [cm$^{-1}$]','Pressure [hPa]','Jacobian [dTb/dx]',\
                        'GSI {} Jacobian SWIR'.format(s),\
                        os.path.join(a.outpath,'{}_gsi_sw.png'.format( s.lower() )),\
                        zlim = [symMinS2,symMaxS2],\
                        figureResolution=300 ) 
            wvTrunc = []
            for w in wv[idxLw]:
                wvTrunc.append('{0:.3f}'.format(w))
            plotContourLabelIdx(wvTrunc, pressure, sValGsi[idxLw],\
                        'Wavenumber [cm$^{-1}$]','Pressure [hPa]','Jacobian [dTb/dx]',\
                        'GSI {} Jacobian LWIR'.format(s),\
                        os.path.join(a.outpath,'{}_gsi_lw.png'.format( s.lower() )),\
                        zlim = [symMinS3,symMaxS3],\
                        figureResolution=300 ) 
            for ii in range(idxLw.shape[0]):
                idxPair = np.array([idxLw[ii],idxSw[ii]]) 
                plotLines ( sValGsi[idxPair], pressure, "Jacobian [dTb/dx]", \
                       "Pressure [hPa]",(wv[idxPair]).tolist(),\
                       'Jacobians', os.path.join(a.outpath,'{}_gsi_sub_{}_{}_lines.png'.format( s.lower(),wv[idxLw[ii]],wv[idxSw[ii]] )),\
                        ylim =[100,1000],figureResolution = 300 )

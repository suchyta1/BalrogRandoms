#!/usr/bin/env python

import numpy as np
import numpy.lib.recfunctions as recfunctions
import os
import sys
import esutil
import healpy as hp

import matplotlib.pyplot as plt


def GetData(version='v3-combined', band='i', method='FITS', catalogdir=os.path.join(os.environ['GLOBALDIR'],'DBFits'), truth='truth', sim='matched', des='des', nosim='nosim'):
    fitscats = [truth, sim, nosim, des]

    # Read from FITS file
    cats = []
    if method.upper()=='FITS':
        for cat in fitscats:
            cat = esutil.io.read(os.path.join(catalogdir, version, '{0}-{1}.fits'.format(cat,band)))
            cats.append(cat)



    #Get data from DB, possibly implement eventually if it's useful
    if method.upper()=='DB':
        pass


    return cats


def PointMap(data, band='i', x='alphawin_j2000', y='deltawin_j2000', ax=None, plotkwargs={}, downfactor=None, downsize=None, title=None):
    x = '{0}_{1}'.format(x,band)
    y = '{0}_{1}'.format(y,band)
  
    if downfactor is not None:
        size = len(data) / downfactor
        keep = np.random.choice(len(data), size=size, replace=False)
    elif downsize is not None:
        keep = np.random.choice(len(data), size=downsize, replace=False)
    else:
        keep = np.ones(len(data), dtype=np._bool)

    if ax is None:
        fig, ax = plt.subplots()

    ax.scatter(data[x][keep],data[y][keep], **plotkwargs)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if title is not None:
        ax.set_title(title)
    return len(data[keep])

def DistributionPlot(data, col, bins, plotkwargs={}):
    center = 1



if __name__=='__main__': 

    band = 'i' 
    truth, sim, nosim, des = GetData(version='v3-combined', band=band, method='FITS', catalogdir=os.environ['DBFITS'])
    #truth, sim, nosim, des = GetData(version='sva1v3_3', band=band, method='FITS')

    fig, axarr = plt.subplots(1,2, figsize=(12,6))
    npoints = PointMap(sim, band=band, downfactor=100, plotkwargs={'lw':0, 's':0.2}, ax=axarr[0], title='Balrog')
    npoints = PointMap(des, band=band, downsize=npoints, plotkwargs={'lw':0, 's':0.2}, ax=axarr[1], title='DES')
    #plt.show()
    fig.suptitle('Raw')
    plt.subplots_adjust(top=0.85)
    plt.savefig('plots/simple-map.png')

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
    hist, bins = np.histogram(data[col], bins)
    dbin = np.diff(bins)
    n = hist / (dbin * np.sum(hist))
    c = (bins[0:-1]+bins[1:])/2.0
    return n, c


def PlotBD(sim, des, col, bins, ax, band='i', plotkwargs={}, title=None):
    col = '{0}_{1}'.format(col,band)
    plotkwargs['color'] = 'red'
    s, c = DistributionPlot(sim, col, bins=bins)
    ax.plot(c, np.log10(s), label='Balrog', **plotkwargs)

    plotkwargs['color'] = 'blue'
    d, c = DistributionPlot(des, col, bins=bins)
    ax.plot(c, np.log10(d), label='DES', **plotkwargs)

    ax.legend(loc='best')
    ax.set_xlabel(col)
    ax.set_ylabel(r'$\log P$')
    if title is not None:
        ax.set_title(title)


def BorisPlot(file):
    pass



if __name__=='__main__': 

    band = 'i' 
    #truth, sim, nosim, des = GetData(version='v3-combined', band=band, method='FITS', catalogdir=os.environ['DBFITS'])
    truth, sim, nosim, des = GetData(version='v3-combined', band=band, method='FITS')
    #truth, sim, nosim, des = GetData(version='sva1v3_3', band=band, method='FITS')

    fig, axarr = plt.subplots(1,2, figsize=(12,6))
    npoints = PointMap(sim, band=band, downfactor=100, plotkwargs={'lw':0, 's':0.2}, ax=axarr[0], title='Balrog')
    npoints = PointMap(des, band=band, downsize=npoints, plotkwargs={'lw':0, 's':0.2}, ax=axarr[1], title='DES')
    fig.suptitle('Raw')
    plt.subplots_adjust(top=0.85)
    #plt.savefig('plots/simple-map.png')
    plt.show()

    '''
    #fig, axarr = plt.subplots(1,3, figsize=(12,6))
    fig, axarr = plt.subplots(1,1, figsize=(6,6))
    bins=np.arange(20,26,0.1)
    #PlotBD(sim, des, 'mag_auto', bins, axarr[0], band=band, plotkwargs={}, title='All')
    PlotBD(sim, des, 'mag_auto', bins, axarr, band=band, plotkwargs={}, title='All')
    #plt.savefig('plots/mag-i-all.png')
    plt.show()
    '''

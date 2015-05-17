#!/usr/bin/env python

import numpy as np
import numpy.lib.recfunctions as recfunctions
import os
import sys
import esutil
import healpy as hp
import pyfits

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



def BorisNative(file):
    data = esutil.io.read(file)
    return data

def BorisAsMap(file, bnside=4096, nest=False):
    data = esutil.io.read(file)
    pix = data['PIXEL']
    value = data['SIGNAL']
    
    map = np.zeros(hp.nside2npix(bnside))
    map[:] = hp.UNSEEN
    map[pix] = value
    return map


def RaDec2Healpix(ra, dec, nside, nest=False):
    phi = np.radians(ra)
    theta = np.radians(90.0 - dec)
    hpInd = hp.ang2pix(nside, theta, phi, nest=nest)
    return hpInd


def GetBorisPlot(map, boris, cat, ra='alphawin_j2000_i', dec='deltawin_j2000_i', borisnside=4096, bins=np.array([0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35])):
    hps = RaDec2Healpix(cat[ra],cat[dec], nside=borisnside)
    nfpix = len(np.unique(hps))
    navg = len(hps) / float(nfpix)
    #relerr_navg2 = 1.0/navg
    relerr_navg2 = 1.0/len(hps)

    bvalue = map[hps]
    avg = np.average(boris['SIGNAL'])
    newvalue = bvalue/avg

    n = np.zeros(len(bins)-1) 
    relerr_n = np.zeros(len(bins)-1)
    v = (bins[1:]+bins[:-1])/2.0
    for i in range(len(bins)-1):
        cut = (newvalue > bins[i]) & (newvalue < bins[i+1])
        num = np.sum(cut)
        nn = num / float( len(np.unique(hps[cut])) )
        n[i] = nn / navg
        relerr_n[i] = np.sqrt(relerr_navg2 + 1.0/num) * n[i]
        print np.sum(cut), relerr_navg2,  1.0/num, relerr_n[i]

    return v,n,relerr_n


def Compare2Boris(sim, des, map, boris, title=None, bins=np.array([0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35])):
    fig, ax = plt.subplots()
    vs, s, se = GetBorisPlot(map, boris, sim, ra='alphawin_j2000_i', dec='deltawin_j2000_i', bins=bins)
    BorisPlot(vs, s, se, ax=ax, plotkwargs={'c':'red', 'label':'Balrog'})
    vd, d, de = GetBorisPlot(map, boris, des, ra='alphawin_j2000_i', dec='deltawin_j2000_i', bins=bins)
    BorisPlot(vd, d, de, ax=ax, plotkwargs={'c':'blue', 'label':'DES'})

    ax.legend(loc='best')
    ax.set_ylabel(r'$n/\bar{n}$')
    ax.set_xlabel(r'$Q/\bar{Q}$')
    if title is not None:
        ax.set_title(title)


def BorisPlot(v,n,e, ax=None, title=None, plotkwargs={}):

    if ax is None:
        fig, ax = plt.subplots()
    
    ax.errorbar(v, n, yerr=e, **plotkwargs)
    print e
    
    if title is not None:
        ax.set_title()



if __name__=='__main__': 

    band = 'i' 
    truth, sim, nosim, des = GetData(version='v3-combined', band=band, method='FITS', catalogdir=os.environ['DBFITS'])
    #truth, sim, nosim, des = GetData(version='v3-combined', band=band, method='FITS')
    #truth, sim, nosim, des = GetData(version='sva1v3_3', band=band, method='FITS')

    sim = sim[ (sim['mag_auto_i'] > 20) & (sim['mag_auto_i'] < 25) ]
    des = des[ (des['mag_auto_i'] > 20) & (des['mag_auto_i'] < 25) ]


    #file = 'nside4096_oversamp4/SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_FWHM__mean.fits'
    #t = 'PSF FWHM'

    #file = 'nside4096_oversamp4/SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_maglimit__.fits.gz'
    #t = 'Mag Limit'

    #file = 'nside4096_oversamp4/SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_SKYBRITE_coaddweights_mean.fits.gz'
    #t = 'Sky Brightness'

    #file = 'nside4096_oversamp4/SVA1_IMAGE_SRC_band_g_nside4096_oversamp4_SKYSIGMA_coaddweights_mean.fits.gz'
    #t = 'Sky Sigma'

    file = 'nside4096_oversamp4/SVA1_IMAGE_SRC_band_i_nside4096_oversamp4_AIRMASS__mean.fits.gz'
    t = 'Airmass'

    map = BorisAsMap(file)
    boris = BorisNative(file)
    #Compare2Boris(sim, des, map, boris, title=r'%s %s-band'%(t,band))
    Compare2Boris(sim, des, map, boris, title=r'%s %s-band'%(t,band), bins=np.array([0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04]))
    plt.savefig('plots/boris-%s-%s.png'%(t,band))


    '''
    fig, axarr = plt.subplots(1,2, figsize=(12,6))
    npoints = PointMap(sim, band=band, downfactor=100, plotkwargs={'lw':0, 's':0.2}, ax=axarr[0], title='Balrog')
    npoints = PointMap(des, band=band, downsize=npoints, plotkwargs={'lw':0, 's':0.2}, ax=axarr[1], title='DES')
    fig.suptitle('Raw')
    plt.subplots_adjust(top=0.85)
    #plt.savefig('plots/simple-map.png')
    plt.show()
    '''

    '''
    #fig, axarr = plt.subplots(1,3, figsize=(12,6))
    fig, axarr = plt.subplots(1,1, figsize=(6,6))
    bins=np.arange(20,26,0.1)
    #PlotBD(sim, des, 'mag_auto', bins, axarr[0], band=band, plotkwargs={}, title='All')
    PlotBD(sim, des, 'mag_auto', bins, axarr, band=band, plotkwargs={}, title='All')
    #plt.savefig('plots/mag-i-all.png')
    plt.show()
    '''

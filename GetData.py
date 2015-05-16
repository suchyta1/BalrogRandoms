#!/usr/bin/env python

import numpy as np
import numpy.lib.recfunctions as recfunctions
import os
import sys
import esutil


if __name__=='__main__': 
    version = 'sva1v3'
    bands = ['g','r','i','z','Y']
    outdir = 'Data'
   

    DBselect = {'table': version,
                'des': 'sva1_coadd_objects',
                'bands': bands,
                'truth': ['balrog_index', 'mag', 'ra', 'dec', 'objtype' ],
                'sim': ['mag_auto', 'flux_auto', 'fluxerr_auto', 'flags', 'spread_model', 'spreaderr_model', 'class_star', 'mag_psf', 'alphawin_j2000', 'deltawin_j2000', 'flux_radius', 'kron_radius', 'ellipticity', 'mag_aper_4', 'fwhmpsf_image', 'fwhmpsf_world']
               }

    truth, sim, nosim, des = DBfunctions.GetAllViaTileQuery(select)
    if MPI.COMM_WORLD.Get_rank()==0:
        outdir = os.path.join(outdir, version)
        if not os.path.exists(dir):
            os.makedirs(dir)
        esutil.io.write(os.path.join(dir,'truth-%s.fits'%(band)), truth, clobber=True)
        esutil.io.write(os.path.join(dir,'matched-%s.fits'%(band)), matched, clobber=True)
        esutil.io.write(os.path.join(dir,'nosim-%s.fits'%(band)), nosim, clobber=True)
        esutil.io.write(os.path.join(dir,'des-%s.fits'%(band)), des, clobber=True)

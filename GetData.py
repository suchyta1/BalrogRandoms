#!/usr/bin/env python

import numpy as np
import numpy.lib.recfunctions as recfunctions
import os
import sys
import esutil
import DBfunctions
import query
from mpi4py import MPI


if __name__=='__main__': 

    band = sys.argv[1]
    DBselect = query.Select()
    DBselect['bands'] = [band]
   
    truth, sim, nosim, des = DBfunctions.GetAllViaTileQuery(DBselect)
    if MPI.COMM_WORLD.Get_rank()==0:
        dir = os.path.join(DBselect['outdir'], DBselect['table'])
        if not os.path.exists(dir):
            os.makedirs(dir)
        esutil.io.write(os.path.join(dir,'truth-{0}.fits'.format(band)), truth, clobber=True)
        esutil.io.write(os.path.join(dir,'matched-{0}.fits'.format(band)), sim, clobber=True)
        esutil.io.write(os.path.join(dir,'nosim-{0}.fits'.format(band)), nosim, clobber=True)

        print len(des)
        esutil.io.write(os.path.join(dir,'des-{0}.fits'.format(band)), des, clobber=True)

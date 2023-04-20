#!/usr/bin/env python
import numpy as np
import pyETS
from pfs.utils.coordinates.CoordTransp import CoordinateTransform as ctrans


def readTargets(fname):
    print (fname)
    with open(fname) as f:
        ll = f.readlines()
        ids = []
        ras = []
        decs = []
        times = []
        pris = []
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_, ra, dec, time, pri = (tt[0], float(tt[1]), float(tt[2]),
                                           float(tt[3]), int(tt[4]))
                ids.append(id_)
                ras.append(ra)
                decs.append(dec)
                times.append(time)
                pris.append(pri)
    return ids, ras, decs, times, pris


# Temporary, very crude method to convert Ra/Dec pairs to x/y coordinates
# on the focal plane. To be replaced by the official functionality once
# available.
# All input angles are expected in degrees.
def radec2pos(ras, decs, raTel=None, decTel=None, posang=0.,
              time="2016-04-03T08:00:00Z"):
    if raTel is None:
        raTel = np.average(ras)
    if decTel is None:
        decTel = np.average(decs)
    tmp = np.array([ras, decs])
    tmp = ctrans(xyin=tmp,
                 mode="sky_pfi", pa=posang,
                 cent=np.array([raTel, decTel]).reshape((2,1)),
                 time=time)
    return tmp[0, :] + 1j*tmp[1, :]


if __name__ == '__main__':
    # get a data structure containing the idealized cobras
    cobras = pyETS.getAllCobras()

    # Parse a target file and return the quantities of interest
    ids, ras, decs, times, pris = readTargets("data/ets_test_data.dat")
    pos = radec2pos(ras, decs)

    # get a list of targets, and a list of Cobras that can observe them
    visibility_map = pyETS.getVis(pos, cobras)

    # perform target assignment using the "draining" algorithm, and return the
    # list of assigned targets and which cobras were used to observe them.
    res = pyETS.getObs(pos, times, pris, cobras, "draining_closest")

    print ("TargetID   Cobra  X         Y          RA         Dec")
    for idx, val in res.items():
        print ("%s %6i %10.5f %10.5f %10.5f %10.5f" % (ids[idx],
              val+1, pos[idx].real, pos[idx].imag, ras[idx],
              decs[idx]))

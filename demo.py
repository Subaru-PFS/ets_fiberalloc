import numpy as np
import pyETS
import pycconv

def readTargets(fname):
    with open(fname) as f:
        ll=f.readlines()
        ids=[]
        ras=[]
        decs=[]
        times=[]
        pris=[]
        for l in ll[1:]:
            if not l.startswith("#"):
                tt = l.split()
                id_,ra,dec,time,pri = tt[0],float(tt[1]),float(tt[2]),float(tt[3]),int(tt[4])
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
    return pycconv.cconv(ras,decs,raTel,decTel,posang,time)

# get a data structure containing the idealized cobras
cobras = pyETS.getAllCobras()

# Parse a target file and return the quantities of interest
ids,ras,decs,times,pris = readTargets("data/ets_test_data.dat")
pos = radec2pos(ras,decs)

# get a list of targets, and a list of Cobras that can observe them
visibility_map=pyETS.getVis(pos,cobras)

# perform target assignment using the "draining" algorithm, and return the list
# of assigned targets and which cobras were used to observe them.
res=pyETS.getObs(pos,times,pris,cobras,"draining_closest")

print "TargetID   Cobra  X         Y          RA         Dec"
for i in range(len(res.keys())):
    idx = res.keys()[i]
    print "%s %6i %10.5f %10.5f %10.5f %10.5f" % (ids[idx], res.values()[i]+1, pos[idx].real, pos[idx].imag, ras[idx], decs[idx])

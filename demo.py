from pyETS import *

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
    return ids,ras,decs,times,pris

# Parse a target file and return the quantities of interest
ids,ras,decs,times,pris = readTargets("cosmology.dat")

# Initialize the ETS class with target data, telescope pointing and observation time
ets=ETShelper(ids,ras,decs,times,pris,getAllCobras(),33.5, -4.,0,"2016-04-03T08:00:00Z")

# get a list of targets, and a list of Cobras that can observe them
visibility_map=ets.getVis()

# perform target asignment using the "draining" algorithm, and return the list
# of assigned targets and which cobras were used to observe them.
res=ets.getObservation("draining")

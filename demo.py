import numpy as np
import pyETS

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

def radec2vec3(ra,dec):
    r2=ra*np.pi/180.
    d2=dec*np.pi/180.
    return np.array([np.cos(d2)*np.cos(r2),np.cos(d2)*np.sin(r2),np.sin(d2)])

# Temporary, very crude method to convert Ra/Dec pairs to x/y coordinates
# on the focal plane. To be replaced by the official functionality once
# available.
# All input angles are expected in degrees.
def radec2pos(ras, decs, raTel=None, decTel=None, raUp=None, decUp=None):
    # convert Ra/Dec to 3-vectors
    ras=np.asarray(ras)
    decs=np.asarray(decs)
    loc=np.empty((len(ras),3),dtype=np.float64)
    for i in range (len(ras)):
        loc[i,:]=radec2vec3(ras[i],decs[i])
    # convert telescope pointing to 3-vector
    # if not specified, use average direction of all targets
    if raTel is not None and decTel is not None:
        vtel=radec2vec3(raTel,decTel)
    else:
        vtel=np.average(loc,axis=0)
        vtel/=np.linalg.norm(vtel)
    # convert "up" direction to 3-vector
    # if not specified, use (0,0,1)
    if raUp is not None and decUp is not None:
        vup=radec2vec3(raUp,decUp)
    else:
        vup=np.array([0.,0.,1.])
    # compute normalized x and y axes
    vx=np.cross(vtel,vup)
    vx/=np.linalg.norm(vx)
    vy=np.cross(vx,vtel)
    # convert target positions to x/y
    rad2mm=180/np.pi*320 #very crude approximation
    pos=(-np.dot(loc,vx) + 1j*np.dot(loc,vy))*rad2mm
    return pos

# get a data structure containing the idealized cobras
cobras = pyETS.getAllCobras()

# Parse a target file and return the quantities of interest
ids,ras,decs,times,pris = readTargets("cosmology.dat")
pos = radec2pos(ras,decs)

# get a list of targets, and a list of Cobras that can observe them
visibility_map=pyETS.getVis(pos,cobras)

# perform target assignment using the "draining" algorithm, and return the list
# of assigned targets and which cobras were used to observe them.
res=pyETS.getObs(pos,times,pris,cobras,"draining_closest")


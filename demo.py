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

def radec2pos(ras, decs, raTel=None, decTel=None, raUp=None, decUp=None):
    ras=np.asarray(ras)*np.pi/180.
    decs=np.asarray(decs)*np.pi/180.
    loc=np.empty((len(ras),3),dtype=np.float64)
    loc[:,0]=np.cos(decs)*np.cos(ras)
    loc[:,1]=np.cos(decs)*np.sin(ras)
    loc[:,2]=np.sin(decs)
    vavg=np.average(loc,axis=0)
    vavg/=np.linalg.norm(vavg)
    xvec=np.cross(vavg,np.array([0.,0.,1.]))
    xvec/=np.linalg.norm(xvec)
    yvec=np.cross(vavg,xvec)
    rad2mm=180/np.pi*320 #very crude approximation
    pos=(-np.dot(loc,xvec) + 1j*np.dot(loc,yvec))*rad2mm
    return pos

# Parse a target file and return the quantities of interest
ids,ras,decs,times,pris = readTargets("cosmology.dat")
pos = radec2pos(ras,decs)

# get a list of targets, and a list of Cobras that can observe them
visibility_map=pyETS.getVis(pos,pyETS.getAllCobras())

# perform target asignment using the "draining" algorithm, and return the list
# of assigned targets and which cobras were used to observe them.
res=pyETS.getObs(pos,times,pris,pyETS.getAllCobras(),"draining_closest")
print res


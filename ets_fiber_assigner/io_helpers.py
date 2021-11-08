import numpy as np
from collections import defaultdict
from pfs.datamodel.utils import calculate_pfsDesignId


def readPfsDesign(pfsDesignId, pfsDesignDirectory="."):
    import pfs.datamodel
    return pfs.datamodelPfsDesign.read(pfsDesignId, pfsDesignDirectory)


def inputParamsFromPfsDesign(pfsDesignId, pfsDesignDirectory):
    import pfs.datamodel
    from .netflow import Telescope
    design = pfs.datamodel.PfsDesign.read(pfsDesignId, pfsDesignDirectory)
    return Telescope(design.raBoresight, design.decBoresight, design.posAng,
                     300.)  # placeholder until we know where to get the observation time


def generatePfsDesign(vis, tp, tel, tgt, classdict):
    import pfs.datamodel

    tdict = defaultdict(int)
    N = len(vis.items())
    ra = []
    dec = []
    pfiNominal = []
    fiberId = []
    objId = []
    targetType = []

    for tidx, cidx in vis.items():
        tdict[tgt[tidx].targetclass] += 1
        ra.append(tgt[tidx].ra)
        dec.append(tgt[tidx].dec)
        fiberId.append(cidx)
        objId.append(tgt[tidx].ID)
        pfiNominal.append([ tp[tidx].real, tp[tidx].imag ])
        targetType.append(classdict[tgt[tidx].targetclass ] )

    # Compute designId from the fiberId, ra and dec for each
    # target. Using datamodel utility function
    pfsDesignId = pfs.datamodel.utils.calculate_pfsDesignId(
        np.asarray(fiberId), np.asarray(ra), np.asarray(dec))

    d = dict(pfsDesignId=pfsDesignId,
        raBoresight=tel._ra,
        decBoresight=tel._dec,
        posAng=tel._posang,
        arms="",
        fiberId=fiberId,
        tract=[np.nan] * N,
        patch=["nan,np.nan"] * N,
        ra=ra,
        dec=dec,
        catId=[np.nan] * N,
        objId=objId,
        targetType=targetType,
        fiberStatus=[1] * N,
        fiberFlux=[[np.nan]] * N,
        psfFlux=[[np.nan]] * N,
        totalFlux=[[np.nan]] * N,
        fiberFluxErr=[[np.nan]] * N,
        psfFluxErr=[[np.nan]] * N,
        totalFluxErr=[[np.nan]] * N,
        filterNames=[['g']] * N,
        pfiNominal=pfiNominal,
        guideStars=None)

    return pfs.datamodel.PfsDesign(**d)


def writePfsDesign(pfsDesignDirectory, vis, tp, tel, tgt, classdict):
    import pfs.datamodel

    pfsDesign = generatePfsDesign(vis, tp, tel, tgt, classdict)
    filename = pfs.datamodel.PfsDesign.fileNameFormat % (pfsDesign.pfsDesignId)
    pfsDesign.write(dirName=pfsDesignDirectory, fileName=filename)

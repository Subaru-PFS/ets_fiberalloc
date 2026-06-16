import numpy as np

def getBench():
    import os
    from ics.cobraOps.Bench import Bench
    from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach
    os.environ["PFS_INSTDATA_DIR"] = "/home/martin/codes/pfs_instdata"
    cobraCoach = CobraCoach(
        loadModel=True, trajectoryMode=True, rootDir="/home/martin/codes/ets_fiberalloc/")
    black_dot_radius_margin = 1.65 # default is 1.0
    bench = Bench(cobraCoach, blackDotsMargin=black_dot_radius_margin)
    print("Number of cobras:", bench.cobras.nCobras)
    return bench

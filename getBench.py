import numpy as np

def getBench():
    import os
    from ics.cobraOps.Bench import Bench
    from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach
    from ics.cobraOps.BlackDotsCalibrationProduct import BlackDotsCalibrationProduct
    os.environ["PFS_INSTDATA_DIR"] = "/home/martin/codes/pfs_instdata"
    cobraCoach = CobraCoach(
        loadModel=True, trajectoryMode=True, rootDir="/home/martin/codes/efa/")
    
    # Get the calibration product
    calibrationProduct = cobraCoach.calibModel
    
    # Fix the phi and tht angles for some of the cobras
    wrongAngles = calibrationProduct.phiIn == 0
    calibrationProduct.phiIn[wrongAngles] = -np.pi
    calibrationProduct.phiOut[wrongAngles] = 0
    calibrationProduct.tht0[wrongAngles] = 0
    calibrationProduct.tht1[wrongAngles] = (2.1 * np.pi) % (2 * np.pi)
    print(f"Number of cobras with wrong phi and tht angles: {np.sum(wrongAngles)}")
    
    # Check if there is any cobra with too short or too long link lengths
    tooShortLinks = np.logical_or(
        calibrationProduct.L1 < 1, calibrationProduct.L2 < 1)
    tooLongLinks = np.logical_or(
        calibrationProduct.L1 > 5, calibrationProduct.L2 > 5)
    print(f"Number of cobras with too short link lenghts: {np.sum(tooShortLinks)}")
    print(f"Number of cobras with too long link lenghts: {np.sum(tooLongLinks)}")
        
    # Load the black dots calibration file
    calibrationFileName = os.path.join(
        os.environ["PFS_INSTDATA_DIR"],"data/pfi/dot", "black_dots_mm.csv")
    blackDotsCalibrationProduct = BlackDotsCalibrationProduct(calibrationFileName)
    
    # Create the bench instance
    black_dot_radius_margin = 1.65 # default is 1.0
    bench = Bench(cobraCoach, blackDotsCalibrationProduct,
                  blackDotsMargin=black_dot_radius_margin)
    print("Number of cobras:", bench.cobras.nCobras)
    return bench

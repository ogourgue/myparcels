import numpy as np

#############################################
# Compute probability density horizontally. #
#############################################

def compute_hprob(LON, LAT, lon, lat):

    # Attention: This function assumes regular spacing for grid coordinates.

    # Horizontal grid size in degrees (assumes regular spacing).
    dLON = (LON[-1] - LON[0]) / (len(LON) - 1)
    dLAT = (LAT[-1] - LAT[0]) / (len(LAT) - 1)

    # Initialize probability density.
    P = np.zeros((len(LAT), len(LON)))

    # Compute grid indices for each particle positions.
    idlon = np.round((lon - LON[0]) / dLON).astype(int)
    idlat = np.round((lat - LAT[0]) / dLAT).astype(int)
    idlon = np.clip(idlon, 0, len(LON) -1)
    idlat = np.clip(idlat, 0, len(LAT) -1)

    # Compute number of particles per grid cell.
    for i in range(len(lon)):
        P[idlat[i], idlon[i]] += 1

    # Divide by total number of particles to obtain probability.
    P /= len(lon)

    # Divide by cell surface area to obtain probability density.
    P /= (1852 * 60) * (1852 * 60) * np.cos(LAT.reshape(-1, 1) * np.pi / 180)

    return P
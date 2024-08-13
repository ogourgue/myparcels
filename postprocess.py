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

    # Loop over horizontal grid cells.
    for i in range(len(LON)):
        for j in range(len(LAT)):
            ind0 = lon > LON[i] - .5 * dLON
            ind1 = lon < LON[i] + .5 * dLON
            ind2 = lat > LAT[j] - .5 * dLAT
            ind3 = lat < LAT[j] + .5 * dLAT
            ind01 = np.logical_and(ind0, ind1)
            ind23 = np.logical_and(ind2, ind3)
            ind = np.logical_and(ind01, ind23)
            P[j, i] = np.mean(ind)

    # Divide by cell surface area to obtain probability density.
    for i in range(len(LAT)):
        P[:, i] /= (1852 * 60) * (1852 * 60) * np.cos(LAT[i] * np.pi / 180)

    return P
import numpy as np

#########################################################
# Compute vertical velocity from horizontal components. #
#########################################################


def compute_w(depth, lat, lon, bathy, u, v, zeta):

    # East neighbor of zonal velocity.
    ue = np.zeros(u.shape)
    ue[:, :, :, :-1] = u[:, :, :, 1:]   # Domain.
    ue[:, :, :, -1] = u[:, :, :, -1]    # East boundary.
    ue[np.isnan(ue)] = -u[np.isnan(ue)] # Land.

    # West neighbor of zonal velocity.
    uw = np.zeros(u.shape)
    uw[:, :, :, 1:] = u[:, :, :, :-1]   # Domain.
    uw[:, :, :, 0] = u[:, :, :, 0]      # West boundary.
    uw[np.isnan(uw)] = -u[np.isnan(uw)] # Land.

    # North neighbor of meridional velocity.
    vn = np.zeros(u.shape)
    vn[:, :, :-1, :] = v[:, :, 1:, :]   # Domain.
    vn[:, :, -1, :] = v[:, :, -1, :]    # North boundary.
    vn[np.isnan(vn)] = -v[np.isnan(vn)] # Land.

    # South neighbor of meridional velocity.
    vs = np.zeros(u.shape)
    vs[:, :, 1:, :] = v[:, :, :-1, :]   # Domain.
    vs[:, :, 0, :] = v[:, :, 0, :]      # South boundary.
    vs[np.isnan(vs)] = -v[np.isnan(vs)] # Land.

    # Horizontal grid size in degrees.
    dx = np.zeros(u.shape[-2:]) + (lon[-1] - lon[0]) / (len(lon) - 1)
    dy = np.zeros(u.shape[-2:]) + (lat[-1] - lat[0]) / (len(lat) - 1)

    # Horizontal grid size in meters.
    dx *= 1852 * 60 * np.cos(lat * np.pi / 180).reshape(-1, 1)
    dy *= 1852 * 60

    # Vertical coordinates at cell edges.
    zs = np.zeros((u.shape[0], len(depth) + 1, len(lat), len(lon)))
    zs[:, 0, :, :] = zeta                                       # Water surface.
    for i in range(len(depth) - 1):
        zs[:, i + 1, :, :] = .5 * (depth[i] + depth[i + 1])      # Water column.
    zs[:, -1, :, :] = bathy   # Bottom.

    # Vertical grid size between cell edges.
    dzs = np.zeros(u.shape)
    dzs = zs[:, 1:, :, :] - zs[:, :-1, :, :]

    # Lateral contributions.
    dudx = (ue - uw) / (2 * dx)
    dvdy = (vn - vs) / (2 * dy)
    dws = dzs * (dudx + dvdy)

    # Initialize vertical velocity.
    w = np.zeros(u.shape) + np.nan

    # Loop over horizontal cells.
    for i in range(len(lon)):
        for j in range(len(lat)):

            # Initialize vertical velocity at cell edges.
            ws = np.zeros((u.shape[0], len(depth) + 1))

            # If not land.
            if np.mean(np.isnan(u[:, :, j, i])) < 1:

                # Vertical velocity is zero at the surface.
                # Loop over other vertical cell edges.
                for k in range(1, len(depth) + 1):

                    # If above bottom.
                    if np.mean(np.isnan(u[:, k - 1, j, i])) < 1:

                        # Vertical velocity at cell edges (integration).
                        ws[:, k] = ws[:, k - 1] + dws[:, k - 1, j, i]

                # Vertical velocity at cell centers (interpolation).
                for l in range(ws.shape[0]):

                    if np.mean(np.isnan(ws[l, :])) < 1:
                        w[l, :, j, i] = np.interp(depth, zs[l, :, j, i], ws[l, :])

    # Vertical velocity should be NaN where horizontal velocity is.
    w[np.isnan(u)] = np.nan

    return w
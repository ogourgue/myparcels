import math

#################################
# Smagorinsky diffusion kernel. #
#################################

def smagorinsky(particle, fieldset, time):

    # Particle coordinates.
    depth = particle.depth
    lat = particle.lat
    lon = particle.lon

    # Grid resolution to compute gradients.
    dx = 0.01

    # Time step.
    dt = particle.dt

    # Compute gradients.
    updx, vpdx = fieldset.UV[time, depth, lat, lon + dx]
    umdx, vmdx = fieldset.UV[time, depth, lat, lon - dx]
    updy, vpdy = fieldset.UV[time, depth, lat + dx, lon]
    umdy, vmdy = fieldset.UV[time, depth, lat - dx, lon]
    dudx = (updx - umdx) / (2 * dx)
    dudy = (updy - umdy) / (2 * dx)
    dvdx = (vpdx - vmdx) / (2 * dx)
    dvdy = (vpdy - vmdy) / (2 * dx)
    dudx2 = dudx ** 2
    dvdy2 = dvdy ** 2

    # Compute diffusivity.
    A = fieldset.cell_areas[time, 0, lat, lon]
    sq_deg_to_sq_m = (1852 * 60) ** 2 * math.cos(lat * math.pi / 180)
    A = A / sq_deg_to_sq_m
    Kh = fieldset.Cs * A * math.sqrt(dudx2 + .5 * (dudy + dvdx) ** 2 + dvdy2)

    # Compute diffusion.
    dlat = ParcelsRandom.normalvariate(0, 1) * math.sqrt(2 * math.fabs(dt) * Kh)
    dlon = ParcelsRandom.normalvariate(0, 1) * math.sqrt(2 * math.fabs(dt) * Kh)
    particle_dlat += dlat
    particle_dlon += dlon

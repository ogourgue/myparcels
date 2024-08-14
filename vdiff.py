import math

###########################
# Euler diffusion kernel. #
###########################

def euler(particle, fieldset, time):

    # depth(n+1) = depth(n) + b * dW(n)
    b = math.sqrt(2 * fieldset.Kv)
    dW = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    particle_ddepth += b * dW

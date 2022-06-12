import numpy as np


def gravitationalConstant():
    return 4.302e-3  # pc(M_solar)^-1 (km/s)^2


def particleRing(nb, radius, gravityCst, mass):
    particles = []
    velocities = []
    theta = 0
    arcLength = (2 * np.pi) / nb
    v = np.sqrt(gravityCst * mass / radius)
    while len(particles) < nb:
        angle = theta * arcLength
        beta = angle + np.pi / 2
        theta += 1
        particles.append([radius * np.cos(angle), radius * np.sin(angle)])
        velocities.append([v * np.cos(beta), v * np.sin(beta)])
    return np.array(particles), np.array(velocities)





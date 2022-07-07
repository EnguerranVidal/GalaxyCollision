import numpy as np


def gravitationalConstant():
    return 1  # pc(M_solar)^-1 (km/s)^2


def initial_trajectory(periapsis, eccentricity, trueAnomaly, mass):
    theta = (2 * np.pi * trueAnomaly) / 360
    a = periapsis / (1 - eccentricity)
    mu = gravitationalConstant() * mass
    p = a * (1 - eccentricity ** 2)
    h = np.sqrt(p * mu)
    r = p / (1 + eccentricity * np.cos(theta))
    X = np.array([r * np.cos(theta), r * np.sin(theta)])
    V = np.array([-mu * np.sin(theta) / h, mu * (eccentricity + np.cos(theta)) / h])
    return X, V

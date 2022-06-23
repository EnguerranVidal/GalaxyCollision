import matplotlib.pyplot as plt
import numpy as np

from sources.common.gravitation import gravitationalConstant
from sources.common.generation import generateArms2D, generateDisk3D, generateDisk2D


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def plotGalaxyDisk2D(nbStars=1000, percentage=35):
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')  # Changing figure to black
    ax = fig.add_subplot(111)
    ax.set_facecolor('xkcd:black')  # Changing background to black
    G = gravitationalConstant()
    positions, velocities, masses = generateDisk2D(nbStars, 1, 3, G)
    X_bary = positions[:int(nbStars * percentage / 100), 0]
    Y_bary = positions[:int(nbStars * percentage / 100), 1]
    X_dark = positions[int(nbStars * percentage / 100):, 0]
    Y_dark = positions[int(nbStars * percentage / 100):, 1]
    ax.scatter(X_dark, Y_dark, s=10, c='b', alpha=0.3)
    ax.scatter(X_bary, Y_bary, s=3, c='r')
    plt.show()


def plotGalaxyArms2D(nbStars=1000):
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')  # Changing figure to black
    ax = fig.add_subplot(111)
    ax.set_facecolor('xkcd:black')  # Changing background to black
    G = gravitationalConstant()
    positions, velocities, masses = generateArms2D(nbStars, 5, 1, 5, 3, 1, G)
    X = positions[:, 0]
    Y = positions[:, 1]
    ax.scatter(X, Y, s=5)
    plt.show()


def plotGalaxyDisk3D(nbStars=1000, percentage=35):
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')  # Changing figure to black
    ax = fig.add_subplot(projection='3d')
    ax.set_facecolor('xkcd:black')  # Changing background to black
    G = gravitationalConstant()

    positions, velocities, masses = generateDisk3D(nbStars, 1, 3, 0.1, G)
    X_bary = positions[:int(nbStars * percentage / 100), 0]
    Y_bary = positions[:int(nbStars * percentage / 100), 1]
    Z_bary = positions[:int(nbStars * percentage / 100), 2]
    X_dark = positions[int(nbStars * percentage / 100):, 0]
    Y_dark = positions[int(nbStars * percentage / 100):, 1]
    Z_dark = positions[int(nbStars * percentage / 100):, 2]

    # ax.scatter(X_dark, Y_dark, Z_dark, s=10, c='b', alpha=0.3)
    ax.scatter(X_bary, Y_bary, Z_bary, s=3, c='r')
    set_axes_equal(ax)
    plt.axis('off')
    plt.show()

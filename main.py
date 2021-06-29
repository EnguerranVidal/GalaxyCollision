from __galaxycollision import *


def main():
    """ Main program """
    # GALAXY 1 : Main Galaxy
    rings1 = np.linspace(0.5, 3, 20)
    particles1 = np.full_like(rings1, 40)
    gal1 = GalaxyRing2D(rings1, particles1, 10000, 5)
    gal1.initial_state([0, 0], [0, 0])

    # GALAXY 2 : Dwarf Galaxy Incoming
    X, V = initial_trajectory(3, 0.6, -135,
                              gal1.M)  # We send the dwarf galaxy on a collision course along an elliptic conic
    rings2 = np.linspace(0.1, 1, 15)
    particles2 = np.full_like(rings2, 20)
    gal2 = GalaxyRing2D(rings2, particles2, 500, 2.5)
    gal2.initial_state(X, V)

    # SIMULATION
    merger = Galaxy_Collision([gal1, gal2])
    merger.RUN(0.01, 2, method='Runge_Kutta')  # We run some calculations
    merger.display(gif_fps=25, gif_duration=4)  # We create a GIF from the showcased Matplotlib animation


if __name__ == "__main__":
    main()

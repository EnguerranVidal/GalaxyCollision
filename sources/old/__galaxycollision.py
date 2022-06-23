import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
import time


################################################################################
# --------- FUNCTIONS ---------#

def str_to_float_list(string):
    L = string.split()
    n = len(L)
    for i in range(n):
        L[i] = float(L[i])
    return L


def session_name():
    t0 = time.time()
    struct = time.localtime(t0)
    string = str(struct.tm_year) + '-'
    n_months = str(struct.tm_mon)  # MONTHS
    if len(n_months) == 1:
        n_months = '0' + n_months
    string = string + n_months + '-'
    n_days = str(struct.tm_mday)  # DAYS
    if len(n_months) == 1:
        n_days = '0' + n_days
    string = string + n_days + '-'
    n_hours = str(struct.tm_hour)  # HOURS
    if len(n_hours) == 1:
        n_hours = '0' + n_hours
    string = string + n_hours + '-'
    n_mins = str(struct.tm_min)  # MINUTES
    if len(n_mins) == 1:
        n_mins = '0' + n_mins
    string = string + n_mins + '-'
    n_secs = str(struct.tm_sec)  # SECONDS
    if len(n_secs) == 1:
        n_secs = '0' + n_secs
    string = string + n_secs + '.txt'
    return string


def gravitational_constant():
    return 4.302e-3  # pc(M_solar)^-1 (km/s)^2


def particle_ring(N, radius, G, M):
    particles = []
    velocities = []
    theta = 0
    arclen = (2 * np.pi) / N
    v = np.sqrt(G * M / radius)
    while len(particles) < N:
        angle = theta * arclen
        beta = angle + np.pi / 2
        theta += 1
        particles.append([radius * np.cos(angle), radius * np.sin(angle)])
        velocities.append([v * np.cos(beta), v * np.sin(beta)])
    return np.array(particles), np.array(velocities)


def velocity(X, G, M):
    V = np.zeros_like(X)
    R = np.sqrt(X[:, 0] ** 2 + X[:, 1] ** 2)
    v = np.sqrt(G * M / R)
    theta = np.arctan(X[:, 1] / X[:, 0])
    beta = theta + np.pi / 2
    V[:, 0] = v * np.cos(beta)
    V[:, 1] = v * np.sin(beta)
    return V


def select_list(l, n):
    m = len(l)
    skip = int(m / n)
    new_l = []
    for i in range(m):
        if i % skip == 0:
            new_l.append(l[i])
    return new_l


def initial_trajectory(periapsis, eccentricity, true_anomaly, M):
    theta = (2 * np.pi * true_anomaly) / 360
    a = periapsis / (1 - eccentricity)
    mu = gravitational_constant() * M
    p = a * (1 - eccentricity ** 2)
    h = np.sqrt(p * mu)
    r = p / (1 + eccentricity * np.cos(theta))
    X = np.array([r * np.cos(theta), r * np.sin(theta)])
    V = np.array([-mu * np.sin(theta) / h, mu * (eccentricity + np.cos(theta)) / h])
    return X, V


################################################################################
# --------- GALAXY RING 2D CLASS ---------#

class GalaxyRing2D:
    def __init__(self, radii, particles, central_mass, halo_r):
        self.center_velocity = None
        self.center_position = None
        assert len(radii) == len(particles), "len(radii): " + str(len(radii)) + " not equal to len(particles): " + str(
            len(particles))
        G = gravitational_constant()
        positions = []
        velocities = []
        for i in range(len(particles)):
            X, V = particle_ring(particles[i], radii[i], G, central_mass)
            positions.append(X)
            velocities.append(V)
        self.X = np.concatenate(positions)
        self.V = np.concatenate(velocities)
        self.G = G
        self.M = central_mass
        self.rings = len(particles)
        self.particles = self.X.shape[0]
        self.initialized = False
        # Dynamic Friction parameters : Galactic Halo
        self.halo_r = halo_r
        self.halo_v = 2 * np.sqrt(self.G * self.M / self.halo_r)

    def initial_state(self, X, V):
        assert self.initialized == False, "ObjectCluster2D's position and velocity are already initialized."
        X = np.array([X])
        V = np.array([V])
        self.X = self.X + X
        self.V = self.V + V
        self.initialized = True
        self.center_position = X
        self.center_velocity = V

    def interior_mass(self, r):
        indices = r < self.halo_r
        masses = np.zeros_like(r)
        if masses[indices].shape != (0,):
            masses[indices] = self.halo_v ** 2 * r[indices] ** 3 / (self.G * (r[indices] + self.halo_r) ** 2)
        if masses[~indices].shape != (0,):
            masses[~indices] = self.M
        return masses

    def density(self, r):
        r_in = 0.99 * r
        r_out = 1.01 * r
        M_in = self.interior_mass(r_in)
        M_out = self.interior_mass(r_out)
        dM = M_out - M_in
        dV = (4 / 3) * np.pi * (r_out ** 3 - r_in ** 3)
        return dM / dV

    def Dyn_friction(self, r, v, vc, M):
        rho = self.density(r)
        v_vec = v - vc
        v = np.linalg.norm(v_vec)
        return -4 * np.pi * self.G ** 2 * 3 * M * rho * v_vec / (1 + v) ** 3


################################################################################
# --------- BASIC GALAXY CLASS ---------#

class GalaxyBasic:
    def __init__(self, X, V, Xc, Vc, central_mass, halo_r):
        assert X.shape == V.shape, "X and V are not of same shape"
        self.initialized = True
        self.G = gravitational_constant()
        self.M = central_mass
        self.X = X
        self.V = V
        self.center_position = Xc
        self.center_velocity = Vc
        # Dynamic Friction parameters : Galactic Halo
        self.halo_r = halo_r
        self.halo_v = 2 * np.sqrt(self.G * self.M / self.halo_r)

    def interior_mass(self, r):
        indices = r < self.halo_r
        masses = np.zeros_like(r)
        if masses[indices].shape != (0,):
            masses[indices] = self.halo_v ** 2 * r[indices] ** 3 / (self.G * (r[indices] + self.halo_r) ** 2)
        if masses[~indices].shape != (0,):
            masses[~indices] = self.M
        return masses

    def density(self, r):
        r_in = 0.99 * r
        r_out = 1.01 * r
        M_in = self.interior_mass(r_in)
        M_out = self.interior_mass(r_out)
        dM = M_out - M_in
        dV = (4 / 3) * np.pi * (r_out ** 3 - r_in ** 3)
        return dM / dV

    def Dyn_friction(self, r, v, vc, M):
        rho = self.density(r)
        v_vec = v - vc
        v = np.linalg.norm(v_vec)
        return -4 * np.pi * self.G ** 2 * 3 * M * rho * v_vec / (1 + v) ** 3


################################################################################
# --------- MERGER ENGINE CLASS ---------#

class MergerEngine:
    def __init__(self, galaxies=None):
        if galaxies is None:
            galaxies = []
        self.galaxies = galaxies
        self.gravitationCst = gravitational_constant()
        # Massless Particles
        galaxies_X = [galaxy.positions for galaxy in self.galaxies]
        galaxies_V = [galaxy.velocities for galaxy in self.galaxies]
        self.massless_X = np.concatenate(galaxies_X)
        self.massless_V = np.concatenate(galaxies_V)
        self.n_massless = self.massless_X.shape[0]
        # Center Masses
        galaxies_Xc = [galaxy.centerPosition for galaxy in self.galaxies]
        galaxies_Vc = [galaxy.centerVelocity for galaxy in self.galaxies]
        self.center_X = np.array(galaxies_Xc)
        self.center_V = np.array(galaxies_Vc)
        self.n_centers = self.center_X.shape[0]

    def change_galaxies(self, galaxies):
        self.__init__(galaxies)

    def massless_acceleration(self, objects_X, center_X):
        a = np.zeros(objects_X.shape)
        n = objects_X.shape[0]
        m = len(self.galaxies)
        for j in range(n):
            for k in range(m):
                v = center_X[k, :] - objects_X[j, :]
                d = np.linalg.norm(v)
                a[j, :] = a[j, :] + self.gravitationCst * (self.galaxies[k].centralMass * v) / (d ** 3)
        return a

    def centers_acceleration(self, center_X, center_V):
        a = np.zeros(center_X.shape)
        n = center_X.shape[0]
        for j in range(n):
            for k in range(n):
                if j != k:
                    v = center_X[k, :] - center_X[j, :]
                    d = np.linalg.norm(v)
                    f = self.galaxies[k].dynamicFriction(d, center_V[j, :], center_V[k, :],
                                                         self.galaxies[j].centralMass)
                    a[j, :] = a[j, :] + self.gravitationCst * (self.galaxies[k].centralMass * v) / (d ** 3) + f
        return a

    def compute(self, dt, method='Euler_explicit'):
        if method == 'Euler_explicit':
            self.compute_euler_explicit(dt)
        if method == 'Euler_semi_implicit':
            self.compute_euler_semi_implicit(dt)
        # if method=='Euler_symplectic':
        # self.compute_euler_symplectic(dt)
        # if method=='Heun':
        # self.compute_Heun(dt)
        if method == 'Runge_Kutta':
            self.compute_Runge_Kutta(dt)

    def compute_euler_explicit(self, dt):
        # Calculating new positions
        new_V = self.massless_V + dt * self.massless_acceleration(self.massless_X, self.center_X)
        new_Vc = self.center_V + dt * self.centers_acceleration(self.center_X, self.center_V)
        new_X = self.massless_X + dt * self.massless_V
        new_Xc = self.center_X + dt * self.center_V
        # Updating positions
        self.massless_X = new_X
        self.massless_V = new_V
        self.center_X = new_Xc
        self.center_V = new_Vc

    def compute_euler_semi_implicit(self, dt):
        # Calculating new positions
        new_V = self.massless_V + dt * self.massless_acceleration(self.massless_X, self.center_X)
        new_Vc = self.center_V + dt * self.centers_acceleration(self.center_X, self.center_V)
        new_X = self.massless_X + dt * new_V
        new_Xc = self.center_X + dt * new_Vc
        # Updating positions
        self.massless_X = new_X
        self.massless_V = new_V
        self.center_X = new_Xc
        self.center_V = new_Vc

    def compute_Runge_Kutta(self, dt):
        # Calculating new positions
        k1_X = self.massless_V * dt  # K1
        k1_Xc = self.center_V * dt
        k1_V = self.massless_acceleration(self.massless_X, self.center_X) * dt
        k1_Vc = self.centers_acceleration(self.center_X, self.center_V) * dt
        k2_X = (self.massless_V + k1_V / 2) * dt  # K2
        k2_Xc = (self.center_V + k1_Vc / 2) * dt
        k2_V = self.massless_acceleration(self.massless_X + k1_X / 2, self.center_X + k1_Xc / 2) * dt
        k2_Vc = self.centers_acceleration(self.center_X + k1_Xc / 2, self.center_V + k1_Vc / 2) * dt
        k3_X = (self.massless_V + k2_V / 2) * dt  # K3
        k3_Xc = (self.center_V + k2_Vc / 2) * dt
        k3_V = self.massless_acceleration(self.massless_X + k2_X / 2, self.center_X + k2_Xc / 2) * dt
        k3_Vc = self.centers_acceleration(self.center_X + k2_Xc / 2, self.center_V + k2_Vc / 2) * dt
        k4_X = (self.massless_V + k3_V) * dt  # K4
        k4_Xc = (self.center_V + k3_Vc) * dt
        k4_V = self.massless_acceleration(self.massless_X + k3_X, self.center_X + k3_Xc) * dt
        k4_Vc = self.centers_acceleration(self.center_X + k3_Xc, self.center_V + k3_Vc) * dt
        # Updating positions
        new_X = self.massless_X + (k1_X + 2 * k2_X + 2 * k3_X + k4_X) / 6
        new_Xc = self.center_X + (k1_Xc + 2 * k2_Xc + 2 * k3_Xc + k4_Xc) / 6
        new_V = self.massless_V + (k1_V + 2 * k2_V + 2 * k3_V + k4_V) / 6
        new_Vc = self.center_V + (k1_Vc + 2 * k2_Vc + 2 * k3_Vc + k4_Vc) / 6
        self.massless_X = new_X
        self.massless_V = new_V
        self.center_X = new_Xc
        self.center_V = new_Vc


################################################################################
# --------- SIMULATION CLASS ---------#

class Galaxy_Collision:
    def __init__(self, galaxies, path):
        self.current_dir = path
        self.engine = MergerEngine(galaxies)
        self.galaxies = galaxies
        self.n_galaxies = len(galaxies)
        self.tags = np.concatenate([np.full((self.galaxies[i].nbParticles,), i) for i in range(len(self.galaxies))])
        self.time = 0
        self.is_new = True
        self.saves_file = None
        self.n_saves = 0

    def new_session(self, galaxies):
        self.__init__(galaxies)

    def load_session(self, filename):
        new_path = self.current_dir + "\\logs\\" + filename
        assert os.path.exists(new_path), "ERROR : File does not exists, try a different name."
        self.saves_file = filename
        self.is_new = False
        n_saves, n_galaxies, galaxies_info = self.load_header_info()
        self.n_saves = n_saves
        self.n_galaxies = n_galaxies
        self.tags = np.concatenate([np.full((galaxies_info[i][1]), i) for i in range(len(self.galaxies))])
        with open(filename, "r") as file:
            file.readline()
            for i in range(self.n_galaxies):  # Skips the header
                file.readline()
            for i in range(self.n_saves):
                data = self.load_state(file)
        galaxies = []
        for i in range(self.n_galaxies):
            indices = self.tags == i
            X = np.zeros(shape=(galaxies_info[i][1], 2))
            X[:, 0], X[:, 1] = data[0][indices], data[1][indices]
            V = np.zeros(shape=(galaxies_info[i][1], 2))
            V[:, 0], V[:, 1] = data[2][indices], data[3][indices]
            Xc = np.array([data[4][i], data[5][i]])
            Vc = np.array([data[6][i], data[7][i]])
            galaxies.append(GalaxyBasic(X, V, Xc, Vc, galaxies_info[i][0], galaxies_info[i][2]))
        self.engine = MergerEngine(galaxies)

    def save_state(self):
        self.n_saves = self.n_saves + 1
        new_path = self.current_dir + "\\logs\\" + self.saves_file
        X, Y, Xp, Yp, Xc, Yc, Xpc, Ypc = [], [], [], [], [], [], [], []
        for j in range(self.engine.n_massless):
            X.append(str(self.engine.massless_X[j, 0]))
            Y.append(str(self.engine.massless_X[j, 1]))
            Xp.append(str(self.engine.massless_V[j, 0]))
            Yp.append(str(self.engine.massless_V[j, 1]))
        for j in range(self.engine.n_centers):
            Xc.append(str(self.engine.center_X[j, 0]))
            Yc.append(str(self.engine.center_X[j, 1]))
            Xpc.append(str(self.engine.center_V[j, 0]))
            Ypc.append(str(self.engine.center_V[j, 1]))
        with open(new_path, 'a') as file:
            file.write(str(self.time) + '\n')
            file.write(' '.join(X) + '\n')
            file.write(' '.join(Y) + '\n')
            file.write(' '.join(Xp) + '\n')
            file.write(' '.join(Yp) + '\n')
            file.write(' '.join(Xc) + '\n')
            file.write(' '.join(Yc) + '\n')
            file.write(' '.join(Xpc) + '\n')
            file.write(' '.join(Ypc) + '\n')

    def load_state(self, file):
        assert not self.is_new, "ERROR : The save file has yet to be rendered, try doing calculations first."
        data = []
        line = file.readline()
        data.append(float(line))
        for i in range(1, 9):
            line = file.readline()
            data.append(np.array(str_to_float_list(line)))
        return data

    def load_header_info(self):
        assert not self.is_new, "ERROR : The save file has yet to be rendered, try doing calculations first."
        new_path = self.current_dir + "\\logs\\" + self.saves_file
        with open(new_path, 'r') as file:
            initial_line = file.readline()
            initial_line = initial_line.split()
            k, l = initial_line.index('NFRAMES'), initial_line.index('NGAL')
            n_saves = int(initial_line[k + 2])
            n_galaxies = int(initial_line[l + 2])
            galaxies = []
            for i in range(n_galaxies):
                line = file.readline()
                line = line.split()
                gal = "GAL" + str(i)
                k = initial_line.index(gal + 'MASS')
                l = initial_line.index('N' + gal)
                m = initial_line.index(gal + 'HALOR')
                mass = float(initial_line[k + 2])
                particles = int(initial_line[l + 2])
                halo_r = float(initial_line[m + 2])
                galaxies.append((mass, particles, halo_r))
        return n_saves, n_galaxies, galaxies

    def RUN(self, dt, T, method='Euler_explicit'):
        if self.is_new:
            logs_path = self.current_dir + "\\logs"
            if not os.path.exists(logs_path):
                os.mkdir(logs_path)
            self.saves_file = session_name()
            new_path = self.current_dir + "\\logs\\" + self.saves_file
            with open(new_path, 'w') as file:
                file.write('NFRAMES : ' + str(self.n_saves) + ' NGAL : ' + str(self.n_galaxies) + '\n')
                for i in range(self.n_galaxies):
                    gal = "GAL" + str(i)
                    file.write(gal + "MASS : " + str(self.galaxies[i].centralMass) +
                               " N" + gal + " : " + str(self.galaxies[i].nbParticles) + " " +
                               gal + "HALOR : " + str(self.galaxies[i].haloRadius) + ' \n')
            self.save_state()
        new_path = self.current_dir + "\\logs\\" + self.saves_file
        initial_time = self.time
        print("########## Beginning Calculations ##########")
        while self.time < initial_time + T:
            self.engine.compute(dt, method=method)
            self.time = self.time + dt
            self.save_state()
        print("########## Calculations Finished ##########")
        with open(new_path, "r") as file:
            lines = file.readlines()
            lines[0] = 'NFRAMES : ' + str(self.n_saves) + ' NGAL : ' + str(self.n_galaxies) + '\n'
        with open(new_path, "w") as file:
            file.writelines(lines)
        self.is_new = False

    def display(self, create_gif=True, gif_name=None, gif_fps=25, gif_duration=15):
        assert not self.is_new, "ERROR : Cannot display animation since no calculations have taken place."
        if create_gif:
            assert gif_fps * gif_duration < self.n_saves, "ERROR : gif_duration or gif_fps is to high for the " \
                                                          "available number of frames "
            if gif_name is None:
                gif_name = self.saves_file.replace('.txt', '.gif')
            gifs_dir = self.current_dir + "\\gifs"
            if not os.path.exists(gifs_dir):
                os.mkdir(gifs_dir)
        Images = []
        print("########## DISPLAYING ##########")
        fig = plt.figure(figsize=(10, 10))
        fig.patch.set_facecolor('xkcd:black')  # Changing figure to black
        ax = fig.add_subplot(111)
        ax.set_facecolor('xkcd:black')  # Changing background to black
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        new_path = self.current_dir + "\\logs\\" + self.saves_file
        gif_path = self.current_dir + "\\gifs\\" + gif_name
        with open(new_path, "r") as file:
            file.readline()
            for i in range(self.n_galaxies):  # Skips the header
                file.readline()
            # PLOTTING FIRST FRAME
            data = self.load_state(file)
            galaxy_particles = []
            for j in range(self.n_galaxies):
                indices = self.tags == j
                X = data[1][indices]
                Y = data[2][indices]
                galaxy_particles.append(ax.scatter(X, Y, s=5))
            galaxy_centers = ax.scatter(data[5], data[6], c="white", s=10)
            fig.show()
            plt.pause(3)
            # MAIN LOOP
            for i in range(1, self.n_saves):
                plt.pause(0.04)
                data = self.load_state(file)
                Offset_c = np.array([data[5], data[6]]).T
                galaxy_centers.set_offsets(Offset_c)
                for j in range(self.n_galaxies):
                    indices = self.tags == j
                    X = data[1][indices]
                    Y = data[2][indices]
                    Offset = np.array([X, Y]).T
                    galaxy_particles[j].set_offsets(Offset)
                plt.draw()
                if create_gif:
                    fig.canvas.draw()
                    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
                    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
                    Images.append(image)
        if create_gif:
            print("########## GIF CREATION ##########")
            n_images = gif_fps * gif_duration
            Images = select_list(Images, n_images)
            imageio.mimsave(gif_path, Images, fps=gif_fps)
            print("########## GIF FINISHED ##########")

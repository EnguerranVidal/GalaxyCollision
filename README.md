# **GALAXY MERGER** 

<p align="center">
  <img src="https://github.com/EnguerranVidal/Galaxy-Collision/blob/main/docs/showcase_gifs/Galaxy_Collision.gif" width="600" height="600">
</p>


This Python program's goal is to simulate the collision of multiple 
galaxies in a 2D plane and create a GIF animation of their interaction.
For this simulation, we will use a modified N-Body Engine we used for the 
**[Planetary-Orbits-Solar-System](https://github.com/EnguerranVidal/Planetary-Orbits-Solar-System)** Python project back in 2020.


## The Galaxies :

The program's galaxies are made of a center massive black hole (white dot). It also contains big amount of massless particles
representing stars that orbit around the center (blue dots). The particles are placed along circles/rings which share the same speed from
Newtonian physics. We then give the galaxies an initial position and speed. An animation of the different speed rings can be found right above before the collision from the second galaxy..

## The Particles :

As we mentioned earlier, the particles are separated into two groups.

- **Central Black Holes :** 

They are attractors and are also slowed down by galaxies' surrounding halos through dynamic friction, more can be found on Thijs Verkade's thesis **["Simulating Galaxy Collisions in Python for Astronomy Education "](https://fse.studenttheses.ub.rug.nl/22594/1/bAST_2020_VerkadeT.pdf)**. However, we will try to explain some of its principles down below. Galaxies are usually surrounded by a "halo" of dark matter and old stars which are distributed in a spherical cloud. When a massive body penetrates it, some dark matter and stars are put behind the massive body and slow it down.
<p align="center">
  <img src="https://github.com/EnguerranVidal/Galaxy-Collision/blob/main/docs/showcase_images/halos.PNG">
</p>


The resulting acceleration is similar to air friction. It is usually called "dynamical friction" and is stronger the highest the massive body speed is. The used equation for this acceleration can be found right below, it direclty affects the massive body's velocity and ease the "merging" process for galaxies since without this force, colliding galaxies would just flung out to infinity from energy conservation. The density of the halo cloud is calculated using the method from the source thesis.
<p align="center">
  <img src="https://github.com/EnguerranVidal/Galaxy-Collision/blob/main/docs/showcase_images/friction.PNG">
</p>


- **Massless Particles :** 

They represent stars which masses are insignificant when compared to the center black holes we discussed earlier. To ease the calculations the N-Body problem poses, we make them massless, which means they are not taken into account in the overall gravity calculations and do not affect the central black holes.

## The N-Body Engine :

The calculations engine has been directly taken from the  **[Planetary-Orbits-Solar-System](https://github.com/EnguerranVidal/Planetary-Orbits-Solar-System)** Python Numerical Physics project I did with Jonathan OERS back in 2020. It was of course modified to fit the needs of this project, especially to separate the calculations of next-steps' positions and velocities for central massive bodies (black holes) and massless particles.
The use of integration schemes to solve the N-Body problem differential equations is implemented through the Euler method and the Runge-Kutta 4 method which are the two most famous. We will try to implement more integratiosn schemes later and also try and compare those to the use of "odeint" from the **scipy** scientific library.

<p align="center">
  <img src="https://github.com/EnguerranVidal/Galaxy-Collision/blob/main/docs/showcase_images/integration.PNG">
</p>

## The Main Simulation Class :

The main class commands and runs calculation done by the above engine. It then stores the resulting data in a txt file put in the **[\logs directory](https://github.com/EnguerranVidal/Galaxy-Collision/tree/main/logs)** with a header meant to portray the simulation parameters for easy replicability and further calculations through a system of "loadable sessions". Such header can be seen right below, it hosts some info about the diferent galaxies.
<p align="center">
  <img src="https://github.com/EnguerranVidal/Galaxy-Collision/blob/main/docs/showcase_images/header.PNG">
</p>
The data can then displayed in a Matplotlib animation. The animation can be stored too as a GIF file through the use of the imageio library. However, as of now, the GIF is not compressed and can therefore be quite heavy, we will try to fix this annoyance in the next update if possible.

# CONTENTS
This repository contains the following files :
- **[main.py](https://github.com/EnguerranVidal/Galaxy-Collision/blob/main/main.py)** : contains an example runnable program which result in the creation of the top page animation of two galaxies colliding.
- **[__galaxycollision.py](https://github.com/EnguerranVidal/Galaxy-Collision/blob/main/__galaxycollision.py)** : The main file containing the classes and fucntions used throughout the project. 
 



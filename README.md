# **GALAXY MERGER** 



This Python program's goal is to simulate the collision of multiple 
galaxies in a 2D plane and create a GIF animation of their interaction.
For this simulation, we will use a modified N-Body Engine we used for the 
**[Planetary-Orbits-Solar-System](https://github.com/EnguerranVidal/Planetary-Orbits-Solar-System)** Python project back in 2020.


## The Galaxies :

The program's galaxies are made of a center massive black hole. It also contains big amount of massless particles
representing stars that orbit around the center. The particles are defined in circles/rings which share the same speed from
Newtonian physics. We build the galaxies like so and then give them an initial position and speed.



## The Engine :

As we mentioned earlier, the particles are separated into two groups.

- Central Black Holes : They are attractors and are also slowed down by galaxies' surrounding halos through dynamic friction.

- Massless Particles : They represent stars which masses are insignificant when compared to the center black holes. 
To ease the calculations the N-Body problem poses, we make them massless, which means they are not taken into account in the overall gravity calculations.


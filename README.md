# moleculardynamics
This is a particle simulation program written in C. It generates a .xyz file for visualization using VMD. 
Verlet method is used for simulation while Berendsen Thermostat is used for the heat bath implementation.
The default simulation is a simulation of 512 Argon atom in a cube with side length 45 angstrom. The initial temperature is 600 K and the heat bath is at 85 K. It runs for 5000 steps with step size 0.01 ps.

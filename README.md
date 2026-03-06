
# Pendulum Impact Solver and Simulation (Python)

## Overview
This project implements a numerical solver and simulation for a rigid-body pendulum system undergoing an impact event.

The solver computes the system state before and after impact, while the simulation visualizes the resulting motion.

## Objectives
- Model the dynamics of a pendulum system
- Compute the system response during an impact event
- Implement a numerical solver for the impact conditions
- Visualize the resulting motion

## Engineering Approach
The system is modeled using rigid body dynamics and conservation laws.

The solver computes:
- angular velocities before and after impact
- impulse forces
- energy transfer during collision

The simulation then integrates the system motion and visualizes the results.

## Tools
- Python
- NumPy
- Matplotlib

## Repository Structure
code/
- impact_solver.py
- simulation.py
- animation.py

images/
- simulation results
- plots

## Results
The simulation provides:
- motion visualization
- system trajectories
- analysis of impact dynamics

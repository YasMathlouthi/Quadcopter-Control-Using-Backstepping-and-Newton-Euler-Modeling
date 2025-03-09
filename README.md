# Quadcopter Control Using Backstepping and Newton-Euler Modeling

## Introduction
This repository contains the MATLAB code and documentation for the control of a quadcopter using the Newton-Euler dynamics modeling and backstepping control method. The project aims to develop a robust control system to manage the complex, non-linear dynamics of a quadcopter.

## Project Structure
- `ARS5_vf.m`: Main MATLAB script for the quadcopter simulation.
- `Guide.txt`: Instructions for configuring the simulation parameters.
- `ars5_report.pdf`: Comprehensive project report detailing the methodology, mathematical modeling, and simulation results.

## Getting Started
To run the simulation:
1. **Clone the repository:**
   ```bash
   git clone https://github.com/YasMathlouthi/Quadcopter-Control-Using-Backstepping-and-Newton-Euler-Modeling
2. **Open MATLAB and navigate to the project directory.**
3. **Run the ARS5_vf.m script to start the simulation:** Modify the trajectory variable in the script to switch between different flight trajectories such as 'circulaire', 'spirale', and 'reference'.
4. **Configuration:**
Adjust the simulation parameters in the ARS5_vf.m script to explore different aspects of quadcopter control:

- g: Gravitational constant.
- m: Mass of the quadcopter.
- I_xx, I_yy, I_zz: Moments of inertia.
- Control gains for different axes and motions.
5. **Simulation:**
The script is configured to simulate various trajectories:
- Circular
- Spiral
- Reference (straight line) Each trajectory tests the quadcopter's ability to follow a set path under controlled conditions
   

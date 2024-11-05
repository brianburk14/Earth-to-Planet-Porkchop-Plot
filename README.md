# About the Earth-to-Planet Porkchop Plot Project on GitHub

## Overview
This project generates porkchop plots to determine the optimal launch engergy (C3) for interplanetary transfers from Earth to a target planet (e.g., Mars or Venus). Porkchop plots are widely used in mission planning to visually represent potential trajectories and assess their feasibility based on launch and arrival dates.

The program utilizes NASA Horizons System website for trajectory data to produce visualizations that guide optimal transfer windows based on delta-v and energy requirements. Two files are read into MATLAB, which are the only inputs needed for the user to generate the plots. The first file contains the data for the state vector for Earth between two given dates. The second file contains the data for the state vector for the target celestial body, assuming to orbit around the sun like Earth, between its two given dates.

When running this program, please allow around five seconds for the plots to load due to its computational intensity. The plots should show up as contours, in the shapes of porkchops where the center of the contours should return the lowest C3 value, resulting in minimal delta_v for departure. You should be able to see the numbers associated with the contours in the plots. Time of flight information will also be provided, via the colored lines. This can help determine the optimal date for launch in case of time requirements.

## Requirements
1. MATLAB (tested on version R2023a)
2. Access to NASA Horizons System for trajectory data: https://ssd.jpl.nasa.gov/horizons/app.html#/

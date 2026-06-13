# Wind Speed Prediction Using H-Infinity Kalman Filter

## Overview

This repository contains the MATLAB implementation and simulation data for the conference paper:

**"Robust Wind Speed Prediction Using H-Infinity Kalman Filter"**

The project focuses on estimating the effective wind speed in wind turbines using three different estimation methods:

* Linear Kalman Filter (KF)
* Extended Kalman Filter (EKF)
* H-Infinity Kalman Filter (HKF)

The results show that the H-Infinity Kalman Filter provides higher estimation accuracy under uncertain wind conditions.

---

## Repository Structure

```
data/
    Wind speed datasets (4.5 m/s and 10.5 m/s)

figures/
    Simulation results and estimation error plots

main/
    MATLAB source codes

paper/
    Conference paper

README.md
```

---

## Simulation Environment

* MATLAB
* Wind turbine model: NREL 5 MW
* Wind profiles:

  * Average wind speed: 4.5 m/s
  * Average wind speed: 10.5 m/s

---

## Results

The performance of three estimation algorithms is compared using Mean Squared Error (MSE).

The H-Infinity Kalman Filter achieves the lowest estimation error among the tested methods.

---

## Conference Paper

This work has been accepted as a conference paper.

---

## Author

Armin Sarkoobi

M.Sc. Student in Electrical Engineering (Control)

University of Birjand

---

## Citation

If you use this work, please cite the corresponding conference paper.

---

## Disclaimer

This repository is intended for research and educational purposes.
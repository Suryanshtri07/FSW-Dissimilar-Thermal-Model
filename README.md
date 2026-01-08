# Analytical Modeling of Heat Generation in Dissimilar FSW (Cu-Al)

## 1. Project Overview
This project implements a numerical approach to estimate the heat generation and temperature distribution during **Friction Stir Welding (FSW)** of dissimilar materials (Copper and Aluminum). Unlike conventional models that assume uniform material properties, this model accounts for the distinct shear stress characteristics at the interface of the Cu-Al system.

## 2. Mathematical Framework
The model utilizes an **inverse estimation method** based on the contact shear stress condition:

$$dQ = \omega \cdot \tau_{interface} \cdot r \cdot dA$$

Where:
* $\omega$ = Angular velocity (rad/s)
* $\tau_{interface}$ = Shear stress (dependent on material: $\tau_{Cu}$ or $\tau_{Al}$)
* $r$ = Radial distance from tool center

The tool geometry is modeled with an **Inverted Parabolic Profile** to simulate the concave shoulder design used in high-quality FSW tools to facilitate material flow.

## 3. Key Features
* **Dissimilar Interface Handling:** The domain is discretized into angular sectors, assigning distinct shear stress values ($\tau_{Cu} = 25 MPa$, $\tau_{Al} = 15 MPa$) to the respective advancing and retreating sides.
* **Geometric Accuracy:** Implements a parabolic Z-profile ($Z = -k(r-R_p)(r-R_s)$) rather than a simplified flat shoulder.
* **Thermal Visualization:** Generates a 3D thermal map showing the asymmetric heat input characteristic of dissimilar joining.

## 4. Results
The simulation highlights the heat concentration on the Copper side due to higher flow stress.
![Thermal Distribution](./results/thermal_plot.png)

## 5. Usage
To run the simulation:
```bash
pip install -r requirements.txt
python src/fsw_model.py

# Ponchon-Savarit Calculator for Binary Distillation

The **Ponchon-Savarit Calculator** is a web-based tool designed for chemical engineers to analyze binary distillation columns using the Ponchon-Savarit method. This method integrates material and energy balances on an enthalpy-composition (H-x-y) diagram, offering greater accuracy for non-ideal mixtures compared to the McCabe-Thiele method. The calculator automates computations for theoretical stages, feed stage location, condenser and reboiler duties, and minimum reflux ratio, with interactive visualizations powered by Plotly.js.

This repository contains:
- **HTML/JavaScript implementation**: A user-friendly web interface for distillation calculations.
- **Python implementation**: A reference version translated from the HTML calculator, using SciPy for numerical computations.

Inspired by the calculator at [https://chemenggcalc.com/ponchon-savarit-diagram-calculator-distillation/](https://chemenggcalc.com/ponchon-savarit-calculator-distillation/), this project is open-source for educational and engineering use.

## Features
- **Inputs**:
  - Process parameters: Feed composition (\( z_F \)), flow rate (\( F \)), distillate/bottoms compositions (\( x_D, x_W \)), feed condition (\( q \)), reflux ratio (\( R \)).
  - VLE/Enthalpy data: Liquid/vapor mole fractions (xData, yData), liquid/vapor enthalpies (\( H_L, H_V \)).
- **Outputs**:
  - Number of theoretical stages and feed stage.
  - Condenser (\( Q_C \)) and reboiler (\( Q_R \)) duties in kW.
  - Minimum reflux ratio (\( R_{\text{min}} \)).
  - Interactive H-x-y diagram with tie lines, operating lines, and difference points (\( \Delta_R, \Delta_S \)).
- **Validation**: Ensures monotonicity of VLE data and valid input ranges.
- **Visualization**: Plotly.js for zoomable, downloadable plots.
- **Python Backend**: Replicates calculations with NumPy, SciPy, and Plotly.

## Demo
Access the calculator online at [https://chemenggcalc.github.io/Ponchon-savarit-calculator/]([https://chemenggcalc.com/ponchon-savarit-calculator-distillation/](https://chemenggcalc.github.io/Ponchon-savarit-calculator/)


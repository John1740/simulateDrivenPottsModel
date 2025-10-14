# simulateDrivenPottsModel

Julia code to Simulate a Potts Model driven by a constant force pointing in clockwise direction. 

The different Simulation types simulate different interaction types
- MF_Gillespie_Sim simulates an exact stochstic trajectory of a finite system with Mean-Field interactions using the gillespie algorithm.
- MF_RungeKutta_Sim simulates the trajectory of an infinite system with Mean-Field interactions at the Thermodynamic Limit using Runge Kutta 4.
- NN_Gilespie_Sim simulates an exact stochstic trajectory of a finite system with Nearest-Neighbor interactions using a unique implementation of the Gillespie algorithm. 
- NN_MonteCarlo_Sim simulates an approximate stochstic trajectory a finite system with Nearest-Neighbor interactions using a time discretization approximation.

# Requirements
- Julia v1.10.10 or higher
  - Package: StatsBase v0.34.6 or higher
  - Package: Distributions v0.25.122 or higher
  - Package: JSON v1.1.0 or higher
  - Package: JLD2 v0.5.15 or higher

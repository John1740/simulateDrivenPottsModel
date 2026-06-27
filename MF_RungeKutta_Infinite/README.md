# MF_RungeKutta_Infinite

This folder contains a Runge–Kutta integration of the infinite-system mean-field Potts model. It solves the deterministic mean-field equations in the thermodynamic limit and is useful for comparing finite-size stochastic dynamics to the infinite-system behavior.

## Usage

- Config: `config.json` (example provided). Configuration options are explained below.
- Run:

```bash
cd MF_RungeKutta_Infinite
julia execute_simulation.jl
```

- Output: saved observables in JLD2 files and `simulation_parameters.json`.

## Configuration
The simulation is configured via `config.json`. Each key is documented with its default value and a short description.

Key descriptions
- `tau` (float, default 1.0): base timescale for transitions; multiplies or rescales the computed update rates.
- `beta` (float, default 1.0): inverse temperature parameter controlling the effective stochasticity in the mean-field equations.
- `driving_force` (float, default 0.0): strength of the external driving force applied to the Potts variables.
- `interaction_strength` (float, default 1.0): mean-field coupling constant that sets interaction energy between oscillators.
- `number_of_states` (int, default 1): number of Potts states (q).
- `total_number_of_episodes` (int, default 1): number of independent simulations or trajectories.
- `total_simulation_time` (float, default 10.0): simulated physical time to run each episode until.
- `skip_simulation_time` (float, default 0.0): initial transient time to skip when saving observables.
- `dt` (float, default 0.01): time step used by the Runge–Kutta integrator.
- `save_location` (string, default empty): directory or prefix where output files will be written; empty uses current folder.
- `rates` (string, default "arrhenius"): rate law used in the mean-field equations.
- `rates_eta` (number, default 0): additional parameter passed to the rate function.
- `seed` (int, default 0): random number generator seed (not always used by deterministic integration; included for consistency).

Observable save counts
Each observable key below controls how many times the observable is saved (integer). A value of `0` disables saving that observable; positive values request that many equidistant saves.

- `occupation_vector` (int): save the mean occupation counts of the q states over time.
- `fourier_modes` (int): save Fourier-mode coefficients of the occupation vector.

Notes
- Output is written in JLD2 format.
- `simulation_parameters.json` is saved alongside the observables for reproducibility.

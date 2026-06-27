# NN_MonteCarlo_Approx

This folder contains an approximate nearest-neighbour Monte Carlo implementation for driven Potts models. The dynamics use discrete time steps, making this simulation a baseline for comparing approximation error and qualitative behavior against the exact nearest-neighbour Gillespie model.

## Usage

- Config: `config.json` (example provided). Configuration options are explained below.
- Run:

```bash
cd NN_MonteCarlo_Approx
julia execute_simulation.jl
```

- Output: saved observables in JLD2 files and `simulation_parameters.json`.

## Configuration
The simulation is configured via `config.json`. Each key is documented with its default value and a short description.

Key descriptions
- `tau` (float, default 1.0): base timescale for transitions.
- `beta` (float, default 1.0): inverse temperature parameter controlling stochasticity in the transition rates.
- `driving_force` (float, default 0.0): strength of the external clockwise driving force applied to the Potts variables.
- `interaction_strength` (float, default 1.0): nearest-neighbour coupling constant (J) that sets interaction energy between neighbours.
- `number_of_states` (int, default 3): number of Potts states (q).
- `number_of_pots` (int, default 1): number of independent Potts systems or replicas simulated together.
- `dimensions` (array of ints, default [1,1]): lattice size in each spatial dimension, e.g. `[Lx, Ly]` or `[Lx, Ly, Lz]`.
- `total_number_of_episodes` (int, default 1): number of independent simulation runs.
- `total_simulation_time` (float, default 10.0): simulated physical time to run each episode until.
- `skip_simulation_time` (float, default 0.0): initial transient time to skip when saving observables.
- `lattice_type` (string, default "cubic"): lattice geometry (options: `cubic`, `hexagonal`).
- `start_position` (string, default "center"): shorthand for initial condition selection (options: `center`, `random`, `ordered`).
- `boundry_conditions` (string, default "periodic"): boundary condition type (options: `periodic`, `wall`).
- `save_location` (string, default empty): directory or prefix where output files will be written; empty uses current folder.
- `rates` (string, default "arrhenius"): rate law to use for transitions.
- `rates_eta` (number, default 0): additional parameter passed to the rate function.
- `seed` (int, default 0): random number generator seed.

Observable save counts
Each observable key below controls how many times the observable is saved (integer). A value of `0` disables saving that observable; positive values request that many equidistant saves.

- `potts_state_array` (int): save full microscopic state snapshots (large files). Use sparingly.
- `occupation_vector` (int): save the histogram / occupation counts of the q states over time.
- `fourier_modes` (int): save Fourier-mode coefficients of the occupation vector.
- `angle_change` (int): save mean angular velocity (or phase change) of selected Fourier modes.
- `fourier_abs_avg` (int): save mean absolute values of Fourier modes and their raw second and fourth moments.

Notes
- Output is written in JLD2 format.
- `simulation_parameters.json` is saved alongside the observables for reproducibility.

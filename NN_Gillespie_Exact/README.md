# NN_Gillespie_Exact 

This folder contains an exact Gillespie implementation of q-state Potts models on a lattice with nearest-neighbour interactions and a constant force pointing clockwise. The implementation uses a novel algorithm designed to preserve exact stochastic dynamics while scaling efficiently to larger lattices (see algorithm on the main README).


## Usage

- Config: `config.json` (example provided). Configuration options are explained below
- Run:

```bash
cd NN_Gillespie_Exact
julia execute_simulation.jl
```

- Output: saved observables in JLD2 files and `simulation_parameters.json`.

##  Configuration
The simulation is configured via `config.json` (example shown below). Each key is documented with its default value and a short description — fix or extend any descriptions as needed. `simulation_parameters.json` is written alongside outputs for reproducibility.


Key descriptions
- `tau` (float, default 1.0): base timescale for transitions; multiplies or rescales the computed holding times.
- `beta` (float, default 1.0): inverse temperature parameter controlling stochasticity in the transition rates.
- `driving_force` (float, default 0.0): strength of the external clockwise driving force (f) applied to the Potts variables.
- `interaction_strength` (float, default 1.0): nearest-neighbour coupling constant (J) that sets interaction energy between neighbours.
- `number_of_states` (int, default 3): number of Potts states (q).
- `number_of_pots` (int, default 1): number of independent Potts systems or replicas simulated together (default 1).
- `dimensions` (array of ints, default [1,1]): lattice size in each spatial dimension, e.g. `[Lx, Ly]` or `[Lx, Ly, Lz]`.
- `total_number_of_episodes` (int, default 1): number of independent simulation runs (ensemble size).
- `total_simulation_time` (float, default 10.0): simulated physical time to run each episode until.
- `skip_simulation_time` (float, default 0.0): initial transient time to skip when saving observables (useful to discard transient).
- `lattice_type` (string, default "cubic"): lattice geometry (options: `cubic`, `hexagonal`).
- `start_position` (string, default "center"): shorthand for initial condition (options: `center`, `random`, `ordered`).
- `boundary_conditions` (string, default "periodic"): boundary condition type (options: `periodic`, `wall`).
- `save_location` (string, default empty): directory or prefix where output files will be written; empty uses current folder.
- `rates` (string, default "arrhenius"): rate law to use for transitions (options: `arrhenius`, `glauber`, `one_direction`).
- `rates_eta` (number, default 0): additional parameter passed to the rate function (model-specific; e.g., asymmetry).
- `seed` (int, default 0): random number generator seed (use 0 or negative to allow non-deterministic runs).

Observable save counts
Each observable key below controls how many times the observable is saved (integer). A value of `0` disables saving that observable; positive integers request that many equidistant saves over the recorded simulation window. Time series observables are saved as JLD2 files.

- `potts_state_array` (int): save full microscopic state snapshots (large files). Use sparingly.
- `occupation_vector` (int): save the histogram / occupation counts of the q states over time.
- `fourier_modes` (int): save Fourier-mode coefficients of the occupation vector.
- `angle_change` (int): save mean angular velocity (or phase change) of selected Fourier modes.
- `fourier_abs_avg` (int): save mean absolute values of Fourier modes and their raw second and fourth moments (for Binder-like analysis).
- `entropy` (int): save entropy-like quantities such as dissipated work, inflow, and activity.
- `magnetization` (int): save global magnetization time series.

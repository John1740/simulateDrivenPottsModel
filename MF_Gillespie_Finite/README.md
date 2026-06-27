# MF_Gillespie_Finite

This folder contains an exact Gillespie implementation for finite mean-field Potts systems. Each oscillator interacts equally with every other oscillator, producing stochastic trajectories that capture finite-size fluctuations in the mean-field limit.

## Usage

- Config: `config.json` (example provided). Configuration options are explained below.
- Run:

```bash
cd MF_Gillespie_Finite
julia execute_simulation.jl
```

- Output: saved observables in JLD2 files and `simulation_parameters.json`.

## Configuration
The simulation is configured via `config.json`. Each key is documented with its default value and a short description.

Key descriptions
- `tau` (float, default 1.0): base timescale for transitions; multiplies or rescales the computed holding times.
- `beta` (float, default 1.0): inverse temperature parameter controlling stochasticity in the transition rates.
- `driving_force` (float, default 0.0): strength of the external driving force applied to the Potts variables.
- `interaction_strength` (float, default 1.0): mean-field coupling constant that sets interaction energy between oscillators.
- `number_of_states` (int, default 3): number of Potts states (q).
- `number_of_pots` (int, default 100): number of oscillators in the finite mean-field system.
- `total_number_of_episodes` (int, default 1): number of independent simulation runs (ensemble size).
- `total_simulation_time` (float, default 10.0): simulated physical time to run each episode until.
- `skip_simulation_time` (float, default 0.0): initial transient time to skip when saving observables.
- `start_position` (string, default "center"): shorthand for initial condition selection (e.g. `center`, `random`, `ordered`).
- `save_location` (string, default empty): directory or prefix where output files will be written; empty uses current folder.
- `rates` (string, default "arrhenius"): rate law to use for transitions.
- `rates_eta` (number, default 0): additional parameter passed to the rate function (model-specific).
- `seed` (int, default 0): random number generator seed (use 0 or negative to allow non-deterministic runs).

Observable save counts
Each observable key below controls how many times the observable is saved (integer). A value of `0` disables saving that observable; positive values request that many equidistant saves over the recorded simulation window.

- `occupation_vector` (int): save the histogram / occupation counts of the q states over time.
- `fourier_modes` (int): save Fourier-mode coefficients of the occupation vector.
- `angle_change` (int): save mean angular velocity (or phase change) of selected Fourier modes.
- `fourier_abs_avg` (int): save mean absolute values of Fourier modes and their raw second and fourth moments.
- `entropy` (int): save entropy-like quantities such as dissipated work, inflow, and activity.

Notes
- Output is written in JLD2 format.
- `simulation_parameters.json` is saved alongside the observables for reproducibility.

module sim_data
    export sim_parameters, sim_state, observable_data

    struct sim_parameters
        tau::Float64
        beta::Float64
        driving_force::Float64
        interaction_strength::Float64
        number_of_states ::Int64
        total_number_of_episodes::Int64
        total_simulation_time::Float64
        skip_simulation_time::Float64
        dt:: Float64
        save_location :: String
        rates :: String
        rates_eta :: Float64
        seed :: Int64
    end

    mutable struct sim_state
        number_of_current_episode::Int64
        episode_time::Float64
        occupation_vector::Vector{Float64}
        rates_function::Function
        current_rates_up::Vector{Float64}
        current_rates_down::Vector{Float64}

    end

    struct observable_number_of_saves
        occupation_vector::Int64
        fourier_modes::Int64
    end

    mutable struct observable_data
        occupation_vector::Array{Float64}
        fourier_modes::Array{Complex{Float64}}
    end

end
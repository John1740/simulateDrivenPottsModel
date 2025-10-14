module sim_data_gb
    export sim_parameters, sim_state, observable_data

    struct sim_parameters
        tau::Float64
        beta::Float64
        driving_force::Float64
        interaction_strength::Float64
        number_of_states ::Int64
        number_of_pots ::Int64
        total_number_of_episodes::Int64
        total_simulation_time::Float64
        skip_simulation_time::Float64
        start_position:: String
        save_location :: String
        rates :: String
        rates_eta :: Float64
        seed :: Int64
        version :: String
    end

    mutable struct sim_state
        number_of_current_episode::Int64
        number_of_current_step::Int64
        episode_time::Float64
        prev_episode_time::Float64
        occupation_vector::Vector{Int64}
        fourier_modes::Vector{Complex{Float64}}
        angle_change::Vector{Float64}
        fourier_abs_avg::Float64
        fourier_abs_avg_square::Float64
        fourier_abs_avg_quartic::Float64
        energy::Float64
        work::Float64
        inflow::Float64
        rates_function::Function
        current_rates_up::Vector{Float64}
        current_rates_down::Vector{Float64}
        k_in_up::Vector{Float64}
        k_in_down::Vector{Float64}

    end

    struct observable_number_of_saves
        occupation_vector::Int64
        fourier_modes::Int64
        angle_change::Int64
        fourier_abs_avg::Int64
        entropy::Int64
    end

    mutable struct observable_data
        occupation_vector::Array{Int64}
        fourier_modes::Array{Complex{Float64}}
        angle_change::Array{Float64}
        fourier_abs_avg::Array{Float64}
        fourier_abs_avg_square::Array{Float64}
        fourier_abs_avg_quartic::Array{Float64}
        energy::Array{Float64}
        work::Array{Float64}
        inflow::Array{Float64}
    end

end
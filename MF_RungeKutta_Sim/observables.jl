module observables

    include("sim_data.jl")
    import ..sim_data: sim_parameters, sim_state, observable_data, observable_number_of_saves

    export initialize_observable_arrays, write_observables_to_array
    
    function initialize_observable_arrays(parameters::sim_parameters,observable_saves::observable_number_of_saves) :: observable_data
        occupation_vector_d=parameters.number_of_states
        fourier_modes_d=Int(div(parameters.number_of_states,2))

        return observable_data(
            zeros(Float64,parameters.total_number_of_episodes,observable_saves.occupation_vector+1,occupation_vector_d),
            zeros(Complex{Float64},parameters.total_number_of_episodes,observable_saves.fourier_modes+1,fourier_modes_d))
    end

    function write_observables_to_array(parameters::sim_parameters,states::sim_state,obserables_saves::observable_number_of_saves,observables::observable_data)

        if obserables_saves.occupation_vector>0
            index=floor((states.episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.occupation_vector)
            if index+1>0
                if index>floor((states.episode_time-parameters.dt-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.occupation_vector)
                    observables.occupation_vector[states.number_of_current_episode,Int(index)+1,:]=calculate_occupation_vector(states)
                end
            end
        end

        if obserables_saves.fourier_modes>0
            index=floor((states.episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.fourier_modes)
            if index+1>0
                if index>floor((states.episode_time-parameters.dt-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.fourier_modes)
                    observables.fourier_modes[states.number_of_current_episode,Int(index)+1,:]=calculate_fourier_modes(states)
                end
            end
        end
        
        return Nothing
    end

    #######

    function calculate_occupation_vector(states::sim_state) :: Vector{Float64}
        return states.occupation_vector
    end

    function calculate_fourier_modes(states::sim_state) :: Vector{Complex{Float64}}
        len=length(states.occupation_vector)
        modes=Int(floor(len/2))
        fourier_modes=zeros(Complex{Float64},modes)

        for k in 1:modes
            for n in 1:len
                fourier_modes[k] += exp(1im * 2 * π * k * n /len) * states.occupation_vector[n]
            end
        end
        return fourier_modes
    end



end # module end
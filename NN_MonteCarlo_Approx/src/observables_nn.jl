module observables_nn

    include("sim_data_nn.jl")
    import ..sim_data_nn: sim_parameters, sim_state, observable_data, observable_number_of_saves

    export initialize_observable_arrays, write_observables_to_array, write_initial_observables_to_array
    
    function initialize_observable_arrays(parameters::sim_parameters,observable_saves::observable_number_of_saves)::observable_data
        potts_state_array_d=parameters.number_of_pots
        occupation_vector_d=parameters.number_of_states
        fourier_modes_d=Int(div(parameters.number_of_states,2))
        
        return observable_data(
            zeros(Int64,parameters.total_number_of_episodes,observable_saves.potts_state_array+1,potts_state_array_d),
            zeros(Int64,parameters.total_number_of_episodes,observable_saves.occupation_vector+1,occupation_vector_d),
            zeros(Complex{Float64},parameters.total_number_of_episodes,observable_saves.fourier_modes+1,fourier_modes_d),
            zeros(Float64,parameters.total_number_of_episodes,observable_saves.angle_change+1,1),
            zeros(Float64,parameters.total_number_of_episodes,observable_saves.fourier_abs_avg+1,1),
            zeros(Float64,parameters.total_number_of_episodes,observable_saves.fourier_abs_avg+1,1),
            zeros(Float64,parameters.total_number_of_episodes,observable_saves.fourier_abs_avg+1,1)
            )
    end

    function write_observables_to_array(parameters::sim_parameters,states::sim_state,obserables_saves::observable_number_of_saves,observables::observable_data)#::Nothing

        if obserables_saves.potts_state_array>0
            index=floor(Int64,(states.episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.potts_state_array)+1
            prev_index=floor(Int64,(states.prev_episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time) * obserables_saves.potts_state_array)+1

            
            if index>0
                if index>prev_index

                    
                    for i in prev_index+1:index
                        if i<=obserables_saves.potts_state_array+1 && i>0
                            observables.potts_state_array[states.number_of_current_episode,i,:]=calculate_potts_state_array(states)
                        end
                    end
                end
            end

        end


        if obserables_saves.occupation_vector>0
            index=floor(Int64,(states.episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.occupation_vector)+1
            prev_index=floor(Int64,(states.prev_episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time) * obserables_saves.occupation_vector)+1


            if index>0
                if index>prev_index
                    for i in prev_index+1:index
                        if i<=obserables_saves.occupation_vector+1 && i>0
                            observables.occupation_vector[states.number_of_current_episode,i,:]=calculate_occupation_vector(states)
                        end
                    end
                end
            end

        end


        if obserables_saves.fourier_modes>0
            index=floor(Int64,(states.episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.fourier_modes)+1
            prev_index=floor(Int64,(states.prev_episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.fourier_modes)+1

            if index>0 
                if index>prev_index
                    for i in prev_index+1:index
                        if i<=obserables_saves.fourier_modes+1 && i>0
                            observables.fourier_modes[states.number_of_current_episode,i,:]=calculate_fourier_modes(states)
                        end
                    end
                end
            end

        end

        if obserables_saves.angle_change>0
            index=floor(Int64,(states.episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.angle_change)+1
            prev_index=floor(Int64,(states.prev_episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.angle_change)+1

            if index>0 
                if index>prev_index
                    for i in prev_index+1:index
                        if i<=obserables_saves.angle_change+1 && i>0
                            observables.angle_change[states.number_of_current_episode,i,1]=calculate_angle_change(parameters,states)
                        end
                    end
                end
            end

        end

        if obserables_saves.fourier_abs_avg>0
            index=floor(Int64,(states.episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.fourier_abs_avg)+1
            prev_index=floor(Int64,(states.prev_episode_time-parameters.skip_simulation_time)/(parameters.total_simulation_time-parameters.skip_simulation_time)* obserables_saves.fourier_abs_avg)+1

            if index>0 
                if index>prev_index
                    for i in prev_index+1:index
                        if i<=obserables_saves.fourier_abs_avg+1 && i>0
                            observables.fourier_abs_avg[states.number_of_current_episode,i,1]=calculate_fourier_abs_avg(parameters,states)
                            observables.fourier_abs_avg_square[states.number_of_current_episode,i,1]=calculate_fourier_abs_avg_square(parameters,states)
                            observables.fourier_abs_avg_quartic[states.number_of_current_episode,i,1]=calculate_fourier_abs_avg_quartic(parameters,states)
                        end
                    end
                end
            end

        end


        return nothing
    end

    function write_initial_observables_to_array(parameters::sim_parameters,states::sim_state,obserables_saves::observable_number_of_saves,observables::observable_data)#::Nothing

        if obserables_saves.potts_state_array>0
            observables.potts_state_array[states.number_of_current_episode,1,:]=calculate_potts_state_array(states)
        end

        if obserables_saves.occupation_vector>0
            observables.occupation_vector[states.number_of_current_episode,1,:]=calculate_occupation_vector(states)
        end


        if obserables_saves.fourier_modes>0
            observables.fourier_modes[states.number_of_current_episode,1,:]=calculate_fourier_modes(states)
        end

        if obserables_saves.angle_change>0
            observables.angle_change[states.number_of_current_episode,1,1]=calculate_angle_change(parameters,states)
        end

        if obserables_saves.fourier_abs_avg>0
            observables.fourier_abs_avg[states.number_of_current_episode,1,1]=calculate_fourier_abs_avg(parameters,states)
            observables.fourier_abs_avg_square[states.number_of_current_episode,1,1]=calculate_fourier_abs_avg_square(parameters,states)
            observables.fourier_abs_avg_quartic[states.number_of_current_episode,1,1]=calculate_fourier_abs_avg_quartic(parameters,states)
        end


        return nothing
    end

    #######

    function calculate_potts_state_array(states::sim_state)::Vector{Int64}
        return states.potts_state_array
    end

    function calculate_occupation_vector(states::sim_state)::Vector{Int64}
        return states.occupation_vector
    end

    function calculate_fourier_modes(states::sim_state)::Vector{Complex{Float64}}
        return states.fourier_modes
    end

    function calculate_angle_change(parameters::sim_parameters,states::sim_state)::Float64
        if states.episode_time-parameters.skip_simulation_time<=0.0
            return 0
        end

        return states.angle_change/(states.episode_time-parameters.skip_simulation_time)
    end

    function calculate_fourier_abs_avg(parameters::sim_parameters,states::sim_state)::Float64
        #return states.fourier_abs_avg/states.number_of_current_step

        if states.episode_time-parameters.skip_simulation_time<=0.0
            return 0
        end
        return states.fourier_abs_avg/(states.episode_time-parameters.skip_simulation_time)
    end

    function calculate_fourier_abs_avg_square(parameters::sim_parameters,states::sim_state)::Float64

        if states.episode_time-parameters.skip_simulation_time<=0.0
            return 0
        end
        return states.fourier_abs_avg_square/(states.episode_time-parameters.skip_simulation_time)
    end

    function calculate_fourier_abs_avg_quartic(parameters::sim_parameters,states::sim_state)::Float64

        if states.episode_time-parameters.skip_simulation_time<=0.0
            return 0
        end
        return states.fourier_abs_avg_quartic/(states.episode_time-parameters.skip_simulation_time)
    end


end # module end
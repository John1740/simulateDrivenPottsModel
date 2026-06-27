module sim_nn
__precompile__()
include("sim_data_nn.jl")

using Random
using Dates
using StatsBase
using Distributions

import .sim_data_nn

include("load_save_nn.jl")
include("rates_nn.jl")
include("observables_nn.jl")
include("lattice_nn.jl")



function initialize_simulation(parameters::sim_data_nn.sim_parameters)::sim_data_nn.sim_state
    if parameters.seed != 0
        Random.seed!(parameters.seed)
    end

    nn_list = zeros(Int64, 1, 1)
    modes = Int(floor(parameters.number_of_states / 2))

    state = sim_data_nn.sim_state(
        1, # current episode
        1, # current step
        0.0, # episode time
        0.0, # prev episode time
        zeros(Int64, parameters.number_of_pots), # potts state array
        nn_list, # neighbor list
        zeros(Float64, 1, 1), #rates function
        zeros(Float64, parameters.number_of_pots), # current rates up
        zeros(Float64, parameters.number_of_pots), # current rates down
        zeros(Int64, parameters.number_of_pots), # d Energy up
        zeros(Int64, parameters.number_of_pots), # d Energy down
        1, # max neighbors
        zeros(Int32, 1, 1), # d_state_array
        zeros(Int32, 1, 1), # d_state_occupation
        zeros(Int32, 1, 1), # d_state_position
        zeros(Int64, parameters.number_of_states), # occupation_vector
        zeros(Complex{Float64}, parameters.number_of_states, modes), # fourier_exponents
        zeros(Complex{Float64}, modes), # fourier_modes
        zeros(Float64, modes), # angle_change
        0.0, # fourier_abs_avg
        0.0, # fourier_abs_avg_square
        0.0, # fourier_abs_avg_quartic
        0.0 # work
    )


    return state

end

function initialize_episode(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)#::Nothing
    #states.occupation_vector=zeros(Float64, parameters.number_of_states)
    states.episode_time = 0.0
    states.prev_episode_time = 0.0
    states.number_of_current_step = 0

    if parameters.start_position == "center"
        states.potts_state_array = rand(1:parameters.number_of_states, parameters.number_of_pots)
    elseif parameters.start_position == "edge"
        states.potts_state_array = ones(Int64, parameters.number_of_pots)
    elseif parameters.start_position == "circle"
        p0 = 1 / 6
        pn = (1 - p0) / (parameters.number_of_states - 1)
        weights = pn * ones(Float64, parameters.number_of_states)
        weights[1] = p0
        states.potts_state_array = sample(1:parameters.number_of_states, Weights(weights), parameters.number_of_pots)
    elseif parameters.start_position == "wave"
        states.potts_state_array = ones(Int64, parameters.number_of_pots)
        for i in 1:parameters.number_of_pots
            states.potts_state_array[i] = min(floor(Int, i * parameters.number_of_states / parameters.number_of_pots) + 1, parameters.number_of_states)
        end
    else
        error("start position not found")
    end

    states.angle_change = zeros(Float64, length(states.fourier_modes))
    states.fourier_abs_avg = 0.0
    states.fourier_abs_avg_square = 0.0
    states.fourier_abs_avg_quartic = 0.0
    states.work = 0.0
    states.fourier_modes = zeros(Complex{Float64}, length(states.fourier_modes))

    calculate_occupation_vector!(parameters, states)
    calculate_all_d_energy!(parameters, states)
    initialize_d_states!(parameters, states)


    return nothing
end



function run_simulation()#:: Nothing
    parameters = load_save_nn.load_simulation_parameters("config.json")
    println("Initializing Simulation")
    states = initialize_simulation(parameters)
    lattice_nn.initialize_lattice(parameters, states)

    rates_nn.get_rates_function(parameters, states)

    observable_saves = load_save_nn.load_observable_number_of_saves("config.json")
    observable_data = observables_nn.initialize_observable_arrays(parameters, observable_saves)


    # counter = 0
    # base_name = "started.txt"
    # file_name = base_name
    # while isfile(file_name)
    #     counter += 1
    #     file_name = "$(base_name[1:end-4])($counter)$(base_name[end-3:end])"
    # end
    # touch(file_name)


    while states.number_of_current_episode <= parameters.total_number_of_episodes
        run_episode(parameters, states, observable_saves, observable_data)
        states.number_of_current_episode += 1
    end
    load_save_nn.save_all(parameters, observable_saves, observable_data)

    # counter = 0
    # base_name = "finished.txt"
    # file_name = base_name
    # while isfile(file_name)
    #     counter += 1
    #     file_name = "$(base_name[1:end-4])($counter)$(base_name[end-3:end])"
    # end
    # touch(file_name)

    return nothing
end

function run_episode(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state, observable_saves::sim_data_nn.observable_number_of_saves, observable_data::sim_data_nn.observable_data)#::Nothing
    println("Initializing Episode ", states.number_of_current_episode, " out of ", parameters.total_number_of_episodes)
    initialize_episode(parameters, states)
    observables_nn.write_initial_observables_to_array(parameters, states, observable_saves, observable_data)

    ##
    modes = Int(floor(parameters.number_of_states / 2))
    #exp_values=zeros(Complex{Float64},parameters.number_of_states,modes)
    states.fourier_exponents = zeros(Complex{Float64}, parameters.number_of_states, modes)
    for k in 1:modes
        for n in 1:parameters.number_of_states
            states.fourier_exponents[n, k] += exp(1im * 2 * π * k * (n - 1) / parameters.number_of_states)
        end
    end

    println("Starting Episode")
    start_time = Dates.now()
    prev_progress = -1
    while states.episode_time < parameters.total_simulation_time

        states.number_of_current_step += 1


        run_step_g(parameters, states, observable_saves, observable_data)


        if mod(states.number_of_current_step, 5000) == 0
            progress = min(floor(Int64, states.episode_time / parameters.total_simulation_time * 100), 100)
            if progress > prev_progress
                print("\r$progress%")
            end
            prev_progress = progress
        end

    end
    progress = min(floor(Int64, states.episode_time / parameters.total_simulation_time * 100), 100)
    print("\r$progress%")
    println("")
    end_time = Dates.now()
    elapsed_time = (end_time - start_time).value / 1000
    println("Elapsed time: ", elapsed_time, "seconds")

    #println(states.episode_time)
    #println(states.occupation_vector)
    #println(states.number_of_current_step)
    println(states.fourier_abs_avg / (states.episode_time - parameters.skip_simulation_time))
    #println(states.fourier_modes)
    #println(states.number_of_current_step)
    #println("Angle Change: ",states.angle_change)


    return nothing
end



function run_step_g(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state, observable_saves::sim_data_nn.observable_number_of_saves, observable_data::sim_data_nn.observable_data)#::Nothing
    states.prev_episode_time = states.episode_time
    transition_up = states.d_state_occupation[:, 1] .* states.rates_function[:, 1]
    transition_down = states.d_state_occupation[:, 2] .* states.rates_function[:, 2]
    k_up = sum(transition_up)
    k_down = sum(transition_down)

    k_out = k_up + k_down

    m = length(transition_up)

    dt = rand(Exponential(1 / k_out))

    states.episode_time += dt

    if states.episode_time + dt > parameters.skip_simulation_time
        if observable_saves.angle_change > 0
            calculate_fourier(parameters, states)
        end
    end



    observables_nn.write_observables_to_array(parameters, states, observable_saves, observable_data)


    up_or_down = sample([1, 0], Weights([k_up, k_down]))

    if up_or_down == 1
        state_of_transition = sample(1:m, Weights(transition_up))
        state_occupation = states.d_state_occupation[state_of_transition, 1]


        index_of_transition = states.d_state_array[state_of_transition, 1, rand(1:state_occupation)]
        #println(index_of_transition)
        #println(states.neighbor_list[index_of_transition,:])
        state_i_old = states.potts_state_array[index_of_transition]
        states.occupation_vector[state_i_old] -= 1
        states.potts_state_array[index_of_transition] = mod(state_i_old, parameters.number_of_states) + 1
        state_i_new = states.potts_state_array[index_of_transition]
        states.occupation_vector[state_i_new] += 1

        states.work += parameters.driving_force

        calculate_one_d_energy!(index_of_transition, parameters, states)


        for j in states.neighbor_list[index_of_transition, :]
            if j == 0
                continue

            end
            #calculate_one_d_energy!(Int32(j),parameters,states)
            state_j = states.potts_state_array[j]
            same = 0
            lower = 0
            upper = 0
            if state_i_new == state_j
                same += 1
            end
            if mod(state_i_new, parameters.number_of_states) + 1 == state_j
                lower += 1
            end
            if mod(state_i_new - 2, parameters.number_of_states) + 1 == state_j
                upper += 1
            end

            if state_i_old == state_j
                same -= 1
            end
            if mod(state_i_old, parameters.number_of_states) + 1 == state_j
                lower -= 1
            end
            if mod(state_i_old - 2, parameters.number_of_states) + 1 == state_j
                upper -= 1
            end

            states.energy_change_up[j] += upper - same
            states.energy_change_down[j] += lower - same

        end

    else
        state_of_transition = sample(1:m, Weights(transition_down))
        state_occupation = states.d_state_occupation[state_of_transition, 2]

        index_of_transition = states.d_state_array[state_of_transition, 2, rand(1:state_occupation)]

        state_i_old = states.potts_state_array[index_of_transition]
        states.occupation_vector[state_i_old] -= 1
        states.potts_state_array[index_of_transition] = mod(state_i_old - 2, parameters.number_of_states) + 1
        state_i_new = states.potts_state_array[index_of_transition]
        states.occupation_vector[state_i_new] += 1

        states.work -= parameters.driving_force

        calculate_one_d_energy!(index_of_transition, parameters, states)

        for j in states.neighbor_list[index_of_transition, :]
            if j == 0
                continue
            end
            #calculate_one_d_energy!(Int32(j),parameters,states)
            state_j = states.potts_state_array[j]
            same = 0
            lower = 0
            upper = 0
            if state_i_new == state_j
                same += 1
            end
            if mod(state_i_new, parameters.number_of_states) + 1 == state_j
                lower += 1
            end
            if mod(state_i_new - 2, parameters.number_of_states) + 1 == state_j
                upper += 1
            end

            if state_i_old == state_j
                same -= 1
            end
            if mod(state_i_old, parameters.number_of_states) + 1 == state_j
                lower -= 1
            end
            if mod(state_i_old - 2, parameters.number_of_states) + 1 == state_j
                upper -= 1
            end

            states.energy_change_up[j] += upper - same
            states.energy_change_down[j] += lower - same
        end
    end

    move_pot_state(index_of_transition, parameters, states)
    for nn in states.neighbor_list[index_of_transition, :]
        move_pot_state(Int32(nn), parameters, states)
    end

    return nothing
end

function debug_check_for_position_error(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)
    for i in 1:parameters.number_of_pots
        if states.energy_change_up[i] + states.max_neighbors + 1 != states.d_state_position[i, 1, 1]
            println("error: Position 1")
            println(i)
        end

        if states.energy_change_down[i] + states.max_neighbors + 1 != states.d_state_position[i, 1, 2]
            println("error: Position 2")
            println(i)
        end


        if states.d_state_array[states.d_state_position[i, 1, 1], 1, states.d_state_position[i, 2, 1]] != i
            println("error: Position 3")
            println(i)
        end
        if states.d_state_array[states.d_state_position[i, 1, 2], 2, states.d_state_position[i, 2, 2]] != i
            println("error: Position 5")
            println(i)
        end

    end

end

function calculate_fourier(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)
    modes = length(states.fourier_modes)
    fourier_modes_new = zeros(Complex{Float64}, modes)

    for k in 1:modes
        for n in 1:parameters.number_of_states
            fourier_modes_new[k] += states.fourier_exponents[n, k] * states.occupation_vector[n] / parameters.number_of_pots
        end
    end


    angle_change = zeros(Float64, modes)
    abs_fourier = 0
    abs_fourier_old = 0
    for i in 1:modes
        abs_fourier = abs(fourier_modes_new[i])
        abs_fourier_old = abs(states.fourier_modes[i])
        if abs_fourier == 0.0 || abs_fourier_old == 0.0
            angle_change[i] = 0.0
        else
            angle_change[i] = (real(fourier_modes_new[i]) * (imag(fourier_modes_new[i] - states.fourier_modes[i])) - imag(fourier_modes_new[i]) * (real(fourier_modes_new[i] - states.fourier_modes[i]))) / (abs_fourier^2)
            # angle_old=angle(states.fourier_modes[end])
            # angle_new=angle(fourier_modes_new[end])
            # angle_change=(angle_new-angle_old)
            # if angle_change > π
            #    angle_change -= 2π
            # elseif angle_change < -π
            #     angle_change += 2π
            # end

        end

        if abs(angle_change[i]) > π / 36
            angle_change[i] = 0
        end
    end

    #states.fourier_abs_avg=states.fourier_abs_avg*(states.number_of_current_step-1)/states.number_of_current_step+abs_fourier*1/states.number_of_current_step


    states.fourier_abs_avg += abs_fourier_old * (states.episode_time - states.prev_episode_time)
    states.fourier_abs_avg_square += (abs_fourier_old^2) * (states.episode_time - states.prev_episode_time)
    states.fourier_abs_avg_quartic += (abs_fourier_old^4) * (states.episode_time - states.prev_episode_time)



    states.fourier_modes = fourier_modes_new

    states.angle_change .+= angle_change

    return Nothing

end

function initialize_d_states!(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)#::Nothing
    M = (2 * states.max_neighbors + 1)
    states.d_state_array = zeros(Int32, M, 2, parameters.number_of_pots)
    states.d_state_occupation = zeros(Int32, M, 2)
    states.d_state_position = zeros(Int32, parameters.number_of_pots, 2, 2)

    for pot in 1:parameters.number_of_pots
        #up
        m = states.energy_change_up[pot] + states.max_neighbors + 1
        states.d_state_occupation[m, 1] += 1
        d_pos = states.d_state_occupation[m, 1]
        states.d_state_array[m, 1, d_pos] = pot
        states.d_state_position[pot, 1, 1] = m
        states.d_state_position[pot, 2, 1] = d_pos
        #down
        m = states.energy_change_down[pot] + states.max_neighbors + 1
        states.d_state_occupation[m, 2] += 1
        d_pos = states.d_state_occupation[m, 2]
        states.d_state_array[m, 2, d_pos] = pot
        states.d_state_position[pot, 1, 2] = m
        states.d_state_position[pot, 2, 2] = d_pos
    end

end

function calculate_occupation_vector!(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)#::Nothing
    states.occupation_vector = zeros(Int64, parameters.number_of_states)
    for pot in 1:parameters.number_of_pots
        states.occupation_vector[states.potts_state_array[pot]] += 1
    end

    return nothing

end


function calculate_all_d_energy!(parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)
    for pot in 1:parameters.number_of_pots
        same_state = 0
        upper_state = 0
        lower_state = 0

        for nn in 1:states.max_neighbors
            if states.neighbor_list[pot, nn] == 0
                continue
            end
            if states.potts_state_array[states.neighbor_list[pot, nn]] == states.potts_state_array[pot]
                same_state += 1
            end
            if states.potts_state_array[states.neighbor_list[pot, nn]] == mod(states.potts_state_array[pot] - 2, parameters.number_of_states) + 1
                lower_state += 1
            end
            if states.potts_state_array[states.neighbor_list[pot, nn]] == mod(states.potts_state_array[pot], parameters.number_of_states) + 1
                upper_state += 1
            end

        end
        states.energy_change_up[pot] = upper_state - same_state
        states.energy_change_down[pot] = lower_state - same_state
    end

    return nothing
end

function calculate_one_d_energy!(pot::Int32, parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)

    same_state = 0
    upper_state = 0
    lower_state = 0
    for nn in 1:states.max_neighbors

        if states.neighbor_list[pot, nn] == 0
            continue
        end
        if states.potts_state_array[states.neighbor_list[pot, nn]] == states.potts_state_array[pot]
            same_state += 1
        end
        if states.potts_state_array[states.neighbor_list[pot, nn]] == mod(states.potts_state_array[pot] - 2, parameters.number_of_states) + 1
            lower_state += 1
        end
        if states.potts_state_array[states.neighbor_list[pot, nn]] == mod(states.potts_state_array[pot], parameters.number_of_states) + 1
            upper_state += 1
        end

    end
    states.energy_change_up[pot] = upper_state - same_state
    states.energy_change_down[pot] = lower_state - same_state


    return nothing
end

function move_pot_state(pot::Int32, parameters::sim_data_nn.sim_parameters, states::sim_data_nn.sim_state)
    #up
    if pot == 0
        return nothing
    end
    for i in 1:2
        prev_pot_d_state = states.d_state_position[pot, 1, i]
        if i == 1
            new_pot_d_state = states.energy_change_up[pot] + states.max_neighbors + 1
        else
            new_pot_d_state = states.energy_change_down[pot] + states.max_neighbors + 1
        end
        if new_pot_d_state != prev_pot_d_state
            if states.d_state_position[pot, 2, i] == states.d_state_occupation[prev_pot_d_state, i]
                states.d_state_array[prev_pot_d_state, i, states.d_state_position[pot, 2, i]] = 0
                states.d_state_occupation[prev_pot_d_state, i] -= 1
            else
                last_p_in_m = states.d_state_array[prev_pot_d_state, i, states.d_state_occupation[prev_pot_d_state, i]]
                states.d_state_position[last_p_in_m, 2, i] = states.d_state_position[pot, 2, i]

                states.d_state_array[prev_pot_d_state, i, states.d_state_position[pot, 2, i]] = last_p_in_m

                states.d_state_array[prev_pot_d_state, i, states.d_state_occupation[prev_pot_d_state, i]] = 0
                states.d_state_occupation[prev_pot_d_state, i] -= 1
            end
            states.d_state_occupation[new_pot_d_state, i] += 1
            states.d_state_array[new_pot_d_state, i, states.d_state_occupation[new_pot_d_state, i]] = pot

            states.d_state_position[pot, 1, i] = new_pot_d_state
            states.d_state_position[pot, 2, i] = states.d_state_occupation[new_pot_d_state, i]
        end
    end

    #down
    # prev_pot_d_state=states.d_state_position[pot,1,2]
    # new_pot_d_state=states.energy_change_down[pot]+states.max_neighbors+1
    # if new_pot_d_state!=prev_pot_d_state
    #     if states.d_state_position[pot,2,2]==states.d_state_occupation[prev_pot_d_state,2]
    #         states.d_state_array[prev_pot_d_state,2,states.d_state_position[pot,2,2]]=0
    #         states.d_state_occupation[prev_pot_d_state,2]-=1
    #     else
    #         states.d_state_array[prev_pot_d_state,2,states.d_state_position[pot,2,2]]= states.d_state_array[prev_pot_d_state,2,states.d_state_occupation[prev_pot_d_state,2]]
    #         states.d_state_array[prev_pot_d_state,2,states.d_state_occupation[prev_pot_d_state,2]]=0
    #         states.d_state_occupation[prev_pot_d_state,2]-=1
    #     end
    #     states.d_state_occupation[new_pot_d_state,2]+=1
    #     states.d_state_array[new_pot_d_state,2,states.d_state_occupation[new_pot_d_state,2]]=pot

    #     states.d_state_position[pot,1,2]=new_pot_d_state
    #     states.d_state_position[pot,2,2]=states.d_state_occupation[new_pot_d_state,2]
    # end

end
end

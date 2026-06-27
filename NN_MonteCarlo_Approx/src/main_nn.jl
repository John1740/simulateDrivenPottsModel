#__precompile__()
module sim_nn

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
       
        nn_list=zeros(Int64,1,1)
        modes=Int(floor(parameters.number_of_states/2))

        state=sim_data_nn.sim_state(
            1,
            1,
            0.0,
            0.0,
            zeros(Int64, parameters.number_of_pots),
            nn_list,
            zeros(Float64,1,1),#rates
            0.0,
            1,#max neigbours
            0,
            0,
            zeros(Int64,parameters.number_of_states),#occupation_vector
            zeros(Complex{Float64},modes),#fourier_modes
            0.0,#angle_change
            0.0,#fourier_abs_avg
            0.0,#fourier_abs_avg_square
            0.0#fourier_abs_avg_quartic
            )
        

        return state

    end

    function initialize_episode(parameters::sim_data_nn.sim_parameters,states::sim_data_nn.sim_state)#::Nothing
        #states.occupation_vector=zeros(Float64, parameters.number_of_states)
        states.episode_time=0
        states.prev_episode_time=0
        states.number_of_current_step=0

        if parameters.start_position=="center"
            states.potts_state_array=rand(1:parameters.number_of_states, parameters.number_of_pots)
        elseif parameters.start_position=="edge"
            states.potts_state_array=ones(Int64, parameters.number_of_pots)
        end

        states.angle_change=0
        states.fourier_abs_avg=0
        states.fourier_abs_avg_square=0
        states.fourier_abs_avg_quartic=0
        states.fourier_modes=zeros(Complex{Float64},length(states.fourier_modes))

        calculate_occupation_vector!(parameters,states)
        
  
        return  nothing
    end



    function run_simulation()#:: Nothing
        parameters=load_save_nn.load_simulation_parameters("config.json")
        println("Initializing Simulation")
        states=initialize_simulation(parameters)
        lattice=lattice_nn.initialize_lattice(parameters,states)

        rates_nn.get_rates_function(parameters,states)

        observable_saves=load_save_nn.load_observable_number_of_saves("config.json")
        observable_data=observables_nn.initialize_observable_arrays(parameters,observable_saves)

        while states.number_of_current_episode<=parameters.total_number_of_episodes
            run_episode(parameters,states,observable_saves,observable_data)
            states.number_of_current_episode+=1
        end
        load_save_nn.save_all(parameters,observable_saves,observable_data)
        return nothing
    end

    function run_episode(parameters::sim_data_nn.sim_parameters,states::sim_data_nn.sim_state,observable_saves::sim_data_nn.observable_number_of_saves,observable_data::sim_data_nn.observable_data)#::Nothing
        println("Initializing Episode ", states.number_of_current_episode ," out of " ,parameters.total_number_of_episodes)
        initialize_episode(parameters,states)
        observables_nn.write_initial_observables_to_array(parameters,states,observable_saves,observable_data)

        ##
        modes=Int(floor(parameters.number_of_states/2))
        exp_values=zeros(Complex{Float64},parameters.number_of_states,modes)

        for k in 1:modes
            for n in 1:parameters.number_of_states
                exp_values[n,k] += exp(1im * 2 * π * k * (n-1) /parameters.number_of_states)
            end
        end

        ##
        println("Starting Episode")
        start_time = Dates.now()
        while states.episode_time<parameters.total_simulation_time
            states.number_of_current_step+=1
            if states.episode_time>=parameters.skip_simulation_time
                calculate_fourier(parameters,states,exp_values)
                observables_nn.write_observables_to_array(parameters,states,observable_saves,observable_data)
            end
            
            run_step_g(parameters,states,observable_saves,observable_data)
            if mod(states.number_of_current_step,5000)==0
                progress=min(floor(Int64,states.episode_time/parameters.total_simulation_time*100),100)
                print("\r$progress%")
            end

        end

        observables_nn.write_observables_to_array(parameters,states,observable_saves,observable_data)

        progress=min(floor(Int64,states.episode_time/parameters.total_simulation_time*100),100)
        print("\r$progress%")
        println("")
        end_time = Dates.now()
        elapsed_time = (end_time - start_time).value / 1000
        println("Elapsed time: ", elapsed_time," seconds")

        #println(states.episode_time)
        #println(states.occupation_vector)
        #println(states.number_of_current_step)
        #println(states.fourier_abs_avg/(states.episode_time-parameters.skip_simulation_time))
        #println(states.number_of_current_step)
        println("Angle Change; ",states.angle_change)
        println("Acceptance Rate ; ",states.accepted/(states.accepted+states.rejected))

        return nothing
    end

    function run_step_g(parameters::sim_data_nn.sim_parameters,states::sim_data_nn.sim_state,observable_saves::sim_data_nn.observable_number_of_saves,observable_data::sim_data_nn.observable_data)#::Nothing
        states.prev_episode_time=states.episode_time
        pot_to_update=rand(1:parameters.number_of_pots)
        up_or_down=rand(1:2)

        same_state=0
        upper_state=0
        lower_state=0

        for nn in 1:states.max_neigbours
            #um ungerade nn und komplizierte bounddry conditions zu berücksichtigen wird ineffizent(doppel) gerechentet
            if states.neigbour_list[pot_to_update,nn]==0
                continue
            elseif states.potts_state_array[states.neigbour_list[pot_to_update,nn]]==states.potts_state_array[pot_to_update]
                same_state+=1
                
            elseif states.potts_state_array[states.neigbour_list[pot_to_update,nn]]==mod(states.potts_state_array[pot_to_update]-2,parameters.number_of_states)+1
                lower_state+=1
            elseif states.potts_state_array[states.neigbour_list[pot_to_update,nn]]==mod(states.potts_state_array[pot_to_update],parameters.number_of_states)+1
                upper_state+=1

            end

        end
        difference_up=upper_state-same_state
        difference_down=lower_state-same_state
        tmp=2
        if up_or_down==1
            p=states.rates_function[difference_up+states.max_neigbours+1,1]/states.max_rate*tmp
            #p=1-exp(-states.rates_function[difference_up+states.max_neigbours+1,1]/states.max_rate*tmp)
        else
            p=states.rates_function[difference_down+states.max_neigbours+1,2]/states.max_rate*tmp
            #p=1-exp(-states.rates_function[difference_down+states.max_neigbours+1,2]/states.max_rate*tmp)
        end

        if p>rand()
            if up_or_down==1
                states.occupation_vector[ states.potts_state_array[pot_to_update]]+=-1
                states.potts_state_array[pot_to_update]=mod(states.potts_state_array[pot_to_update],parameters.number_of_states)+1
                states.occupation_vector[ states.potts_state_array[pot_to_update]]+=1
            else
                states.occupation_vector[ states.potts_state_array[pot_to_update]]+=-1
                states.potts_state_array[pot_to_update]=mod(states.potts_state_array[pot_to_update]-2,parameters.number_of_states)+1
                states.occupation_vector[ states.potts_state_array[pot_to_update]]+=1
            end
            states.accepted+=1
        else
            states.rejected+=1
        end

        states.episode_time+=1/(2*parameters.number_of_pots*states.max_rate)*tmp

    end


    function calculate_fourier(parameters::sim_data_nn.sim_parameters,states::sim_data_nn.sim_state,exp_values::Array{Complex{Float64}})
        modes=length(states.fourier_modes)
        fourier_modes_new=zeros(Complex{Float64},modes)

        for k in 1:modes
            for n in 1:parameters.number_of_states
                fourier_modes_new[k] += exp_values[n,k] * states.occupation_vector[n]/parameters.number_of_pots
            end
        end

        abs_foruier=abs(fourier_modes_new[end])
        abs_foruier_old=abs(states.fourier_modes[end])
        
        if abs_foruier==0.0 || abs_foruier_old==0.0
            angle_change=0
        else
            angle_change=(real(fourier_modes_new[end])*(imag(fourier_modes_new[end]-states.fourier_modes[end]))-imag(fourier_modes_new[end])*(real(fourier_modes_new[end]-states.fourier_modes[end])))/(abs_foruier^2)
            # angle_old=angle(states.fourier_modes[end])
            # angle_new=angle(fourier_modes_new[end])
            # angle_change=(angle_new-angle_old)
            # if angle_change > π
            #    angle_change -= 2π
            # elseif angle_change < -π
            #     angle_change += 2π
            # end
        
        end

        if abs(angle_change)> π/36
            angle_change=0
        end

        #states.fourier_abs_avg=states.fourier_abs_avg*(states.number_of_current_step-1)/states.number_of_current_step+abs_foruier*1/states.number_of_current_step
    
        
        states.fourier_abs_avg+=abs_foruier_old*(states.episode_time-states.prev_episode_time)
        states.fourier_abs_avg_square+=(abs_foruier_old^2)*(states.episode_time-states.prev_episode_time)
        states.fourier_abs_avg_quartic+=(abs_foruier_old^4)*(states.episode_time-states.prev_episode_time)



        states.fourier_modes=fourier_modes_new

        states.angle_change+=angle_change

        return Nothing

    end

    function calculate_occupation_vector!(parameters::sim_data_nn.sim_parameters,states::sim_data_nn.sim_state)#::Nothing
        for pot in 1:parameters.number_of_pots
            states.occupation_vector[states.potts_state_array[pot]]+=1
        end

        return nothing

    end

end


module sim_gb
    __precompile__()
    include("sim_data_gb.jl")

    using Random
    using Dates
    using StatsBase
    using Distributions

    import .sim_data_gb

    include("load_save_gb.jl")
    include("rates_gb.jl")
    include("observables_gb.jl")


    
    function initialize_simulation(parameters::sim_data_gb.sim_parameters)::sim_data_gb.sim_state
        if parameters.seed !=0
            Random.seed!(parameters.seed)
        end
        modes=Int(floor(parameters.number_of_states/2))
        return sim_data_gb.sim_state(
            1,
            1,
            0.0,
            0.0,
            zeros(Int64, parameters.number_of_states),
            zeros(Complex{Float64},modes),
            zeros(Float64,modes),#angle change
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            rates_gb.get_rates_function(parameters),
            zeros(Float64, parameters.number_of_states),
            zeros(Float64, parameters.number_of_states),
            zeros(Float64, parameters.number_of_states),
            zeros(Float64, parameters.number_of_states))
    end

    function initialize_episode(parameters::sim_data_gb.sim_parameters,state::sim_data_gb.sim_state)#::Nothing
        #state.occupation_vector=zeros(Float64, parameters.number_of_states)
        state.episode_time=0
        state.prev_episode_time=0
        state.number_of_current_step=0
        
        if parameters.start_position=="center"
            v = div(parameters.number_of_states, parameters.number_of_pots)
            state.occupation_vector=fill(v,parameters.number_of_states)
            for _ in 1:parameters.number_of_pots-v*parameters.number_of_states
                state_index = rand(1:parameters.number_of_states)
                state.occupation_vector[state_index] += 1
            end
        elseif parameters.start_position=="edge"
            state.occupation_vector=zeros(Int64,parameters.number_of_states)#[parameters.number_of_pots,0,0]
            state.occupation_vector[1]=parameters.number_of_pots
        end
        state.fourier_modes=zeros(Complex{Float64},length(state.fourier_modes))

        #for k in 1:length(state.fourier_modes)
        #    for n in 1:parameters.number_of_states
        #        state.fourier_modes[k] += exp(1im * 2 * π * k * (n-1) /parameters.number_of_states) * state.occupation_vector[n]/parameters.number_of_pots
        #    end
        #end

        state.angle_change=zeros(Float64,length(state.fourier_modes))
        state.fourier_abs_avg=0.0
        state.fourier_abs_avg_square=0.0
        state.fourier_abs_avg_quartic=0.0
        state.work=0.0
        state.inflow=0.0
        state.energy=0.0

        return  nothing
    end



    function run_simulation()#:: Nothing
        parameters=load_save_gb.load_simulation_parameters("config.json")
        println("Initializing Simulation")
        states=initialize_simulation(parameters)
        observable_saves=load_save_gb.load_observable_number_of_saves("config.json")
        observable_data=observables_gb.initialize_observable_arrays(parameters,observable_saves)

        counter = 0
        base_name = "started.txt"
        file_name = base_name
        while isfile(file_name)
            counter += 1
            file_name = "$(base_name[1:end-4])($counter)$(base_name[end-3:end])"
        end
        touch(file_name)


        while states.number_of_current_episode<=parameters.total_number_of_episodes
            run_episode(parameters,states,observable_saves,observable_data)
            states.number_of_current_episode+=1
        end
        print(observable_data.work)
        load_save_gb.save_all(parameters,observable_saves,observable_data)


        counter = 0
        base_name = "finished.txt"
        file_name = base_name
        while isfile(file_name)
            counter += 1
            file_name = "$(base_name[1:end-4])($counter)$(base_name[end-3:end])"
        end
        touch(file_name)


        return nothing
    end

    function run_episode(parameters::sim_data_gb.sim_parameters,states::sim_data_gb.sim_state,observable_saves::sim_data_gb.observable_number_of_saves,observable_data::sim_data_gb.observable_data)#::Nothing
        println("Initializing Episode ", states.number_of_current_episode ," out of " ,parameters.total_number_of_episodes)
        initialize_episode(parameters,states)
        observables_gb.write_initial_observables_to_array(parameters,states,observable_saves,observable_data)

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

        calculate_rates!(states)
        #calculate_k_in!(states)
        prev_progress=-1
        while states.episode_time<parameters.total_simulation_time
            states.number_of_current_step+=1
            if states.episode_time>parameters.skip_simulation_time
                calculate_fourier(parameters,states,exp_values)
            end
            run_step(parameters,states,observable_saves,observable_data)
            ##
            if mod(states.number_of_current_step,5000)==0
                progress=min(floor(Int64,states.episode_time/parameters.total_simulation_time*100),100)
                if progress>prev_progress
                    print("\r$progress%")
                end
                prev_progress=progress
            end
            
            ##
        end
        progress=min(floor(Int64,states.episode_time/parameters.total_simulation_time*100),100)
        print("\r$progress%")
        println("")
        end_time = Dates.now()
        elapsed_time = (end_time - start_time).value / 1000
        println("Elapsed time: ", elapsed_time," seconds")

        println("Angle Change; ",states.angle_change)

        #println(observable_data.fourier_abs_avg)
        #println(states.angle_change)
        #println(states.k_in_up)
        #println(states.k_in_down)

        return nothing
    end

    function calculate_fourier(parameters::sim_data_gb.sim_parameters,states::sim_data_gb.sim_state,exp_values::Array{Complex{Float64}})
        modes=length(states.fourier_modes)
        fourier_modes_new=zeros(Complex{Float64},modes)

        for k in 1:modes
            for n in 1:parameters.number_of_states
                fourier_modes_new[k] += exp_values[n,k] * states.occupation_vector[n]/parameters.number_of_pots
            end
        end

        abs_foruier=abs(fourier_modes_new[end])
        
        angle_change=zeros(Float64,modes)
        abs_foruier=0
        abs_foruier_old=0
        for i in 1:modes
            abs_foruier=abs(fourier_modes_new[i])
            abs_foruier_old=abs(states.fourier_modes[i])
            if abs_foruier==0.0 || abs_foruier_old==0.0
                angle_change[i]=0.0
            else
                    angle_change[i]=(real(fourier_modes_new[i])*(imag(fourier_modes_new[i]-states.fourier_modes[i]))-imag(fourier_modes_new[i])*(real(fourier_modes_new[i]-states.fourier_modes[i])))/(abs_foruier^2)
                # angle_old=angle(states.fourier_modes[end])
                # angle_new=angle(fourier_modes_new[end])
                # angle_change=(angle_new-angle_old)
                # if angle_change > π
                #    angle_change -= 2π
                # elseif angle_change < -π
                #     angle_change += 2π
                # end
            
            end

            if abs(angle_change[i])> π/36
                angle_change[i]=0
            end
        end
        #states.fourier_abs_avg=states.fourier_abs_avg*(states.number_of_current_step-1)/states.number_of_current_step+abs_foruier*1/states.number_of_current_step
 
        abs_foruier_old=abs(states.fourier_modes[end])

        dt=(states.episode_time-states.prev_episode_time)
        states.fourier_abs_avg+=abs_foruier_old*dt
        states.fourier_abs_avg_square+=(abs_foruier_old^2)*dt
        states.fourier_abs_avg_quartic+=(abs_foruier_old^4)*dt
        states.energy-=(sum(states.occupation_vector.^2)-parameters.number_of_pots)*parameters.interaction_strength/2/parameters.number_of_pots*dt



        states.fourier_modes=fourier_modes_new

        states.angle_change+=angle_change

        return nothing

    end

    function run_step(parameters::sim_data_gb.sim_parameters,states::sim_data_gb.sim_state,observable_saves::sim_data_gb.observable_number_of_saves,observable_data::sim_data_gb.observable_data)#::Nothing
        
        #calculate_rates!(states) # updated with faster version for q>3
        #calculate_k_in!(states)
        states.prev_episode_time=states.episode_time
        k_up = sum(states.current_rates_up)
        k_down = sum(states.current_rates_down)
        
        k_out=k_up+k_down
        states.episode_time += rand(Exponential(1 / k_out))
        states.inflow=-k_out

        observables_gb.write_observables_to_array(parameters,states,observable_saves,observable_data)

        up_or_down = sample([1, 0], Weights([k_up, k_down]))
        
        if up_or_down==1
            index_of_transition=sample(1:parameters.number_of_states,Weights(states.current_rates_up))
            states.occupation_vector[index_of_transition]+=-1
            states.occupation_vector[mod(index_of_transition,parameters.number_of_states)+1]+=1

            states.work+=parameters.driving_force
           
            update_rates!(states,parameters,index_of_transition,1)
            #update_k_in!(states,parameters,index_of_transition,1)
        else
            index_of_transition=sample(1:parameters.number_of_states,Weights(states.current_rates_down))
            states.occupation_vector[index_of_transition]+=-1
            states.occupation_vector[mod(index_of_transition-2,parameters.number_of_states)+1]+=1

            states.work-=parameters.driving_force
            update_rates!(states,parameters,index_of_transition,-1)
            #update_k_in!(states,parameters,index_of_transition,-1)

        end


        return nothing
    end

    function calculate_rates!(states::sim_data_gb.sim_state)#::Nothing
        states.current_rates_up[1]=states.rates_function(states.occupation_vector[2]-states.occupation_vector[1],1)*states.occupation_vector[1]
        states.current_rates_down[1]=states.rates_function(states.occupation_vector[end]-states.occupation_vector[1],-1)*states.occupation_vector[1]
        for i in 2:length(states.occupation_vector)-1
            states.current_rates_up[i]=states.rates_function(states.occupation_vector[i+1]-states.occupation_vector[i],1)*states.occupation_vector[i]
            states.current_rates_down[i]=states.rates_function(states.occupation_vector[i-1]-states.occupation_vector[i],-1)*states.occupation_vector[i]
        end
        states.current_rates_up[end]=states.rates_function(states.occupation_vector[1]-states.occupation_vector[end],1)*states.occupation_vector[end]
        states.current_rates_down[end]=states.rates_function(states.occupation_vector[end-1]-states.occupation_vector[end],-1)*states.occupation_vector[end]

        return nothing
    end

    function calculate_k_in_old!(states::sim_data_gb.sim_state)#::Nothing
        if states.occupation_vector[1]==0
            states.k_in_up[1]=0
            states.k_in_down[1]=0
        else
            states.k_in_up[1]=states.rates_function(states.occupation_vector[1]-states.occupation_vector[end]-2,1)*(states.occupation_vector[end]+1)
            states.k_in_down[1]=states.rates_function(states.occupation_vector[1]-states.occupation_vector[2]-2,-1)*(states.occupation_vector[2]+1)
        end

        for i in 2:length(states.occupation_vector)-1
            if states.occupation_vector[i]==0
                states.k_in_up[i]=0
                states.k_in_down[i]=0
            else
                states.k_in_up[i]=states.rates_function(states.occupation_vector[i]-states.occupation_vector[i-1]-2,1)*(states.occupation_vector[i-1]+1)
                states.k_in_down[i]=states.rates_function(states.occupation_vector[i]-states.occupation_vector[i+1]-2,-1)*(states.occupation_vector[i+1]+1)
            end
        end

        if states.occupation_vector[end]==0
            states.k_in_up[end]=0
            states.k_in_down[end]=0
        else
            states.k_in_up[end]=states.rates_function(states.occupation_vector[end]-states.occupation_vector[end-1]-2,1)*(states.occupation_vector[end-1]+1)
            states.k_in_down[end]=states.rates_function(states.occupation_vector[end]-states.occupation_vector[1]-2,-1)*(states.occupation_vector[1]+1)
        end

        return nothing
    end

    function update_k_in_old!(states::sim_data_gb.sim_state,parameters::sim_data_gb.sim_parameters,ind::Int64,up_or_down::Int64)#::Nothing

        if up_or_down==1
            ind_l1=mod(ind-2,parameters.number_of_states)+1
            ind_u1=mod(ind,parameters.number_of_states)+1
            ind_u2=mod(ind+1,parameters.number_of_states)+1

            if states.occupation_vector[ind]==0
                states.k_in_up[ind]=0
            else
                states.k_in_up[ind]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_l1]-2,1)*(states.occupation_vector[ind_l1]+1)
                states.k_in_down[ind]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_u1]-2,-1)*(states.occupation_vector[ind_u1]+1)
            end
            if states.occupation_vector[ind_u1]==0
                states.k_in_up[ind_u1]=0
            else
                states.k_in_up[ind_u1]=states.rates_function(states.occupation_vector[ind_u1]-states.occupation_vector[ind]-2,1)*(states.occupation_vector[ind]+1)
                states.k_in_down[ind_u1]=states.rates_function(states.occupation_vector[ind_u1]-states.occupation_vector[ind_u2]-2,-1)*(states.occupation_vector[ind_u2]+1)
            end
          
            if states.occupation_vector[ind_u2]==0
                states.k_in_up[ind_u2]=0
            else
                states.k_in_up[ind_u2]=states.rates_function(states.occupation_vector[ind_u2]-states.occupation_vector[ind_u1]-2,1)*(states.occupation_vector[ind_u1]+1)
            end
                
            if states.occupation_vector[ind_l1]==0
                states.k_in_up[ind_l1]=0
            else
                states.k_in_down[ind_l1]=states.rates_function(states.occupation_vector[ind_l1]-states.occupation_vector[ind]-2,-1)*(states.occupation_vector[ind]+1)
            end

        else
            ind_l1=mod(ind-2,parameters.number_of_states)+1
            ind_l2=mod(ind-3,parameters.number_of_states)+1
            ind_u1=mod(ind,parameters.number_of_states)+1
            

            if states.occupation_vector[ind]==0
                states.k_in_up[ind]=0
            else
                states.k_in_up[ind]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_l1]-2,1)*(states.occupation_vector[ind_l1]+1)
                states.k_in_down[ind]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_u1]-2,-1)*(states.occupation_vector[ind_u1]+1)
            end

            if states.occupation_vector[ind_u1]==0
                states.k_in_up[ind_u1]=0
            else
                states.k_in_up[ind_u1]=states.rates_function(states.occupation_vector[ind_u1]-states.occupation_vector[ind]-2,1)*(states.occupation_vector[ind]+1)
            end
          
            if states.occupation_vector[ind_l2]==0
                states.k_in_up[ind_l2]=0
            else
                states.k_in_down[ind_l2]=states.rates_function(states.occupation_vector[ind_l2]-states.occupation_vector[ind_l1]-2,-1)*(states.occupation_vector[ind_l1]+1)

            end
                
            if states.occupation_vector[ind_l1]==0
                states.k_in_up[ind_l1]=0
            else
                states.k_in_up[ind_l1]=states.rates_function(states.occupation_vector[ind_l1]-states.occupation_vector[ind_l2]-2,1)*(states.occupation_vector[ind_l2]+1)
                states.k_in_down[ind_l1]=states.rates_function(states.occupation_vector[ind_l1]-states.occupation_vector[ind]-2,-1)*(states.occupation_vector[ind]+1)
            end
            
        end
        return nothing
    end

    function update_rates!(states::sim_data_gb.sim_state,parameters::sim_data_gb.sim_parameters,ind::Int64,up_or_down::Int64)#::Nothing

        if up_or_down==1
            ind_l1=mod(ind-2,parameters.number_of_states)+1
            ind_u1=mod(ind,parameters.number_of_states)+1
            ind_u2=mod(ind+1,parameters.number_of_states)+1

            states.current_rates_up[ind]=states.rates_function(states.occupation_vector[ind_u1]-states.occupation_vector[ind],1)*states.occupation_vector[ind]
            states.current_rates_up[ind_l1]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_l1],1)*states.occupation_vector[ind_l1]
            states.current_rates_up[ind_u1]=states.rates_function(states.occupation_vector[ind_u2]-states.occupation_vector[ind_u1],1)*states.occupation_vector[ind_u1]

            states.current_rates_down[ind]=states.rates_function(states.occupation_vector[ind_l1]-states.occupation_vector[ind],-1)*states.occupation_vector[ind]
            states.current_rates_down[ind_u1]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_u1],-1)*states.occupation_vector[ind_u1]
            states.current_rates_down[ind_u2]=states.rates_function(states.occupation_vector[ind_u1]-states.occupation_vector[ind_u2],-1)*states.occupation_vector[ind_u2]
        else
            ind_l1=mod(ind-2,parameters.number_of_states)+1
            ind_l2=mod(ind-3,parameters.number_of_states)+1
            ind_u1=mod(ind,parameters.number_of_states)+1
            

            states.current_rates_up[ind]=states.rates_function(states.occupation_vector[ind_u1]-states.occupation_vector[ind],1)*states.occupation_vector[ind]
            states.current_rates_up[ind_l1]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_l1],1)*states.occupation_vector[ind_l1]
            states.current_rates_up[ind_l2]=states.rates_function(states.occupation_vector[ind_l1]-states.occupation_vector[ind_l2],1)*states.occupation_vector[ind_l2]

            states.current_rates_down[ind]=states.rates_function(states.occupation_vector[ind_l1]-states.occupation_vector[ind],-1)*states.occupation_vector[ind]
            states.current_rates_down[ind_u1]=states.rates_function(states.occupation_vector[ind]-states.occupation_vector[ind_u1],-1)*states.occupation_vector[ind_u1]
            states.current_rates_down[ind_l1]=states.rates_function(states.occupation_vector[ind_l2]-states.occupation_vector[ind_l1],-1)*states.occupation_vector[ind_l1]

        end


        return nothing
    end
    return nothing
end

module main
    using Random
    import .sim_data

    include("sim_data.jl")
    include("load_save.jl")
    include("rates.jl")
    include("observables.jl")


    
    function initialize_simulation(parameters::sim_data.sim_parameters) :: sim_data.sim_state

        if parameters.seed !=0
            Random.seed!(parameters.seed)
        end

        return sim_data.sim_state(1,0,zeros(Float64, parameters.number_of_states),rates.get_rates_function(parameters),zeros(Float64, parameters.number_of_states),zeros(Float64, parameters.number_of_states))
        
    end

    function initialize_episode(parameters::sim_data.sim_parameters,state::sim_data.sim_state)
        #state.occupation_vector=zeros(Float64, parameters.number_of_states)
        state.episode_time=0
        #state.occupation_vector=fill(1/parameters.number_of_states, parameters.number_of_states)
        #state.occupation_vector+=0.3/parameters.number_of_states*rand(parameters.number_of_states)
        #state.occupation_vector./=sum(state.occupation_vector)
        state.occupation_vector=[0.01, 0.01, 0.98]

        return Nothing
    end



    function run_simulation()
        parameters=load_save.load_simulation_parameters("config.json")
        states=initialize_simulation(parameters)
        observable_saves=load_save.load_observable_number_of_saves("config.json")
        observable_data=observables.initialize_observable_arrays(parameters,observable_saves)
        while states.number_of_current_episode<=parameters.total_number_of_episodes
            run_episode(parameters,states,observable_saves,observable_data)
            states.number_of_current_episode+=1
        end
        #println(observable_data.occupation_vector)
        load_save.save_all(parameters,observable_saves,observable_data)
        return Nothing
    end

    function run_episode(parameters::sim_data.sim_parameters,states::sim_data.sim_state,observable_saves::sim_data.observable_number_of_saves,observable_data::sim_data.observable_data)
        initialize_episode(parameters,states)
        while states.episode_time<parameters.total_simulation_time
            observables.write_observables_to_array(parameters,states,observable_saves,observable_data)
            run_step(parameters,states)
        end
        observables.write_observables_to_array(parameters,states,observable_saves,observable_data)
        return Nothing
    end

    function run_step(parameters::sim_data.sim_parameters,states::sim_data.sim_state)
        RK4(parameters,states)
        #println(states.occupation_vector)
        states.episode_time+=parameters.dt
        return Nothing
    end

    function RK1(parameters::sim_data.sim_parameters,states::sim_data.sim_state)
        calculate_rates!(states)
        dp=calculate_dp(states)
        states.occupation_vector+=dp*parameters.dt
        #states.occupation_vector./=sum(states.occupation_vector)

        return Nothing

    end


    function RK4(parameters::sim_data.sim_parameters,states::sim_data.sim_state)
        calculate_transition_rates!(states,states.occupation_vector)
        k1=calculate_k_i(states)*parameters.dt

        calculate_transition_rates!(states,states.occupation_vector+k1/2)
        k2=calculate_k_i(states)*parameters.dt

        calculate_transition_rates!(states,states.occupation_vector+k2/2)
        k3=calculate_k_i(states)*parameters.dt

        calculate_transition_rates!(states,states.occupation_vector+k3)
        k4=calculate_k_i(states)*parameters.dt

        states.occupation_vector+=1/6*(k1+2*k2+2*k3+k4)
        #states.occupation_vector./=sum(states.occupation_vector)

        return Nothing

    end


    function calculate_rates!(states::sim_data.sim_state)
        states.current_rates_up[1]=states.rates_function(states.occupation_vector[2]-states.occupation_vector[1],1)
        states.current_rates_down[1]=states.rates_function(states.occupation_vector[end]-states.occupation_vector[1],-1)
        for i in 2:length(states.occupation_vector)-1
            states.current_rates_up[i]=states.rates_function(states.occupation_vector[i+1]-states.occupation_vector[i],1)
            states.current_rates_down[i]=states.rates_function(states.occupation_vector[i-1]-states.occupation_vector[i],-1)
        end
        states.current_rates_up[end]=states.rates_function(states.occupation_vector[1]-states.occupation_vector[end],1)
        states.current_rates_down[end]=states.rates_function(states.occupation_vector[end-1]-states.occupation_vector[end],-1)

        return Nothing
    end


    function calculate_transition_rates!(states::sim_data.sim_state,occupation_vector::Vector{Float64})
        states.current_rates_up[1]=states.rates_function(occupation_vector[2]-occupation_vector[1],1)
        states.current_rates_down[1]=states.rates_function(occupation_vector[end]-occupation_vector[1],-1)
        for i in 2:length(occupation_vector)-1
            states.current_rates_up[i]=states.rates_function(occupation_vector[i+1]-occupation_vector[i],1)
            states.current_rates_down[i]=states.rates_function(occupation_vector[i-1]-occupation_vector[i],-1)
        end
        states.current_rates_up[end]=states.rates_function(occupation_vector[1]-occupation_vector[end],1)
        states.current_rates_down[end]=states.rates_function(occupation_vector[end-1]-occupation_vector[end],-1)
        
        states.current_rates_up=states.current_rates_up.*occupation_vector
        states.current_rates_down=states.current_rates_down.*occupation_vector
        return Nothing
    end
        
    function calculate_k_i(states::sim_data.sim_state) :: Vector{Float64}
        dp=-states.current_rates_up-states.current_rates_down
        dp[1]+=states.current_rates_up[end]+states.current_rates_down[2]
        dp[end]+=states.current_rates_up[end-1]+states.current_rates_down[1]
        for i in 2:length(states.occupation_vector)-1
            dp[i]+=states.current_rates_up[i-1]+states.current_rates_down[i+1]
        end
        return dp
    end

    function calculate_dp(states::sim_data.sim_state) :: Vector{Float64}
        transition_rates_up=states.current_rates_up.*states.occupation_vector
        transition_rates_down=states.current_rates_down.*states.occupation_vector
        dp=-transition_rates_up-transition_rates_down
        dp[1]+=transition_rates_up[end]+transition_rates_down[2]
        dp[end]+=transition_rates_up[end-1]+transition_rates_down[1]
        for i in 2:length(states.occupation_vector)-1
            dp[i]+=transition_rates_up[i-1]+transition_rates_down[i+1]
        end
        return dp
    end

end

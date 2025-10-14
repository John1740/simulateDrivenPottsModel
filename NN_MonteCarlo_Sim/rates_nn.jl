
module rates_nn
   
    include("sim_data_nn.jl")

    import ..sim_data_nn: sim_parameters, sim_state

    export get_rates_function

    function get_rates_function(parameters::sim_parameters,states::sim_state)#::Nothing

        
        eval(Symbol(parameters.rates*"_rates"))(parameters,states)
        return nothing
        
    end


    function arrhenius_rates(parameters::sim_parameters,states::sim_state)#::Nothing
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        number_of_pots=parameters.number_of_pots
        driving_force=parameters.driving_force
        eta=parameters.rates_eta

        max_neigbours=states.max_neigbours
        
        function rates_f(d::Int64,s::Int64)::Float64
            return 1 / tau * exp(beta / 2 * (interaction_strength / max_neigbours *
            (1 - s * eta) * (d) + s * driving_force))
        end

        rates=zeros(Float64,2*max_neigbours+1,2)
        for i in -max_neigbours:max_neigbours
            rates[i+max_neigbours+1,1]=rates_f(i,1)
            rates[i+max_neigbours+1,2]=rates_f(i,-1)
        end

        states.max_rate=rates_f(max_neigbours,1)
        states.rates_function=rates
        return nothing

    end 

    function one_direction_rates(parameters::sim_parameters,states::sim_state)#::Nothing
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        number_of_pots=parameters.number_of_pots
        driving_force=parameters.driving_force


        max_neigbours=states.max_neigbours

        function rates_f(d::Int64,s::Int64)::Float64
            return exp(interaction_strength / max_neigbours * d)
        end

        rates=zeros(Float64,2*max_neigbours+1,2)
        for i in -max_neigbours:max_neigbours
            rates[i+max_neigbours+1,1]=rates_f(i,1)
            rates[i+max_neigbours+1,2]=0
        end

        states.max_rate=rates_f(max_neigbours,1)
        states.rates_function=rates
        return nothing

    end 

    function glauber_rates(parameters::sim_parameters) :: Function
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        driving_force=parameters.driving_force
        eta=parameters.rates_eta

        max_neigbours=states.max_neigbours

        function rates_f(d::Int64,s::Int64) :: Float64
            return 2 / tau * exp(s * beta / 2 * driving_force)/(exp(-interaction_strength/parameters.max_neigbours * d)+1)
        end
        
        rates=zeros(Float64,2*max_neigbours+1,2)
        for i in -max_neigbours:max_neigbours
            rates[i+max_neigbours+1,1]=rates_f(i,1)
            rates[i+max_neigbours+1,2]=rates_f(i,-1)
        end

        states.max_rate=rates_f(max_neigbours,1)
        states.rates_function=rates
        return nothing

    end 


end
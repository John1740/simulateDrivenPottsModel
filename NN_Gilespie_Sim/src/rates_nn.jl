module rates_nn

    include("sim_data_nn.jl")

    import ..sim_data_nn: sim_parameters, sim_state
    using JLD2

    export get_rates_function

    function get_rates_function(parameters::sim_parameters,states::sim_state)#::Nothing

        
        eval(Symbol(parameters.rates*"_rates"))(parameters,states)
        calculate_w0(parameters,states)
        return nothing
        
    end

    function calculate_w0(parameters::sim_parameters,states::sim_state)
        function prob(i)
            k::Int64=abs(i)
            p=0
            for r::Int64 in k:floor(states.max_neighbors/2+k/2)#range(k,int(np.floor(nn/2+k/2)+1)):
                g::Int64=r-k
                p+=binomial(states.max_neighbors, r+g)*(2/parameters.number_of_states)^(r+g)*(1-2/parameters.number_of_states)^(states.max_neighbors-r-g)*binomial(r+g, r)*(1/2)^(r+g)
            end
            return p

        end

        w0=0
        for i in -states.max_neighbors:states.max_neighbors
            w0+=prob(i)*(states.rates_function[i+states.max_neighbors+1,1]-states.rates_function[i+states.max_neighbors+1,2])

        end


        @save "decoh_work.jld2" w0
    end


    function arrhenius_rates(parameters::sim_parameters,states::sim_state)#::Nothing
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        number_of_pots=parameters.number_of_pots
        driving_force=parameters.driving_force
        eta=parameters.rates_eta

        max_neighbors=states.max_neighbors
        
        function rates_f(d::Int64,s::Int64)::Float64
            return 1 / tau * exp(beta / 2 * (interaction_strength / max_neighbors *
            (1 - s * eta) * (d) + s * driving_force))
        end

        rates=zeros(Float64,2*max_neighbors+1,2)
        for i in -max_neighbors:max_neighbors
            rates[i+max_neighbors+1,1]=rates_f(i,1)
            rates[i+max_neighbors+1,2]=rates_f(i,-1)
        end

        states.rates_function=rates
        return nothing

    end 

    function arrhenius_inf_rates(parameters::sim_parameters,states::sim_state)#::Nothing
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        number_of_pots=parameters.number_of_pots
        driving_force=parameters.driving_force
        eta=parameters.rates_eta

        max_neighbors=states.max_neighbors
        
        function rates_f(d::Int64,s::Int64)::Float64
            return 1 / tau * exp(beta / 2 * (interaction_strength / max_neighbors *
            (1 - s * eta) * (d) ))
        end

        rates=zeros(Float64,2*max_neighbors+1,2)
        for i in -max_neighbors:max_neighbors
            rates[i+max_neighbors+1,1]=rates_f(i,1)
            rates[i+max_neighbors+1,2]=0
        end

        states.rates_function=rates
        return nothing

    end 

    function one_direction_rates(parameters::sim_parameters,states::sim_state)#::Nothing
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        number_of_pots=parameters.number_of_pots
        driving_force=parameters.driving_force


        max_neighbors=states.max_neighbors

        function rates_f(d::Int64,s::Int64)::Float64
            return exp(interaction_strength / max_neighbors * d)
        end

        rates=zeros(Float64,2*max_neighbors+1,2)
        for i in -max_neighbors:max_neighbors
            rates[i+max_neighbors+1,1]=rates_f(i,1)
            rates[i+max_neighbors+1,2]=0
        end

        states.rates_function=rates
        return nothing

    end 

    function glauber_rates(parameters::sim_parameters,states::sim_state) #:: Function
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        driving_force=parameters.driving_force
        eta=parameters.rates_eta

        max_neighbors=states.max_neighbors

        function rates_f(d::Int64,s::Int64) :: Float64
            return 2 / tau * exp(s * beta / 2 * driving_force)/(exp(-interaction_strength/max_neighbors * d)+1)
        end
        
        rates=zeros(Float64,2*max_neighbors+1,2)
        for i in -max_neighbors:max_neighbors
            rates[i+max_neighbors+1,1]=rates_f(i,1)
            rates[i+max_neighbors+1,2]=rates_f(i,-1)
        end


        states.rates_function=rates
        return nothing

    end 


end
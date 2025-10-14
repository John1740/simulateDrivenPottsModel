
module rates_gb
   
    include("sim_data_gb.jl")

    import ..sim_data_gb: sim_parameters

    export get_rates_function

    function get_rates_function(parameters::sim_parameters)::Function

        
        rates =eval(Symbol(parameters.rates*"_rates"))(parameters)
        return rates
    end


    function arrhenius_rates(parameters::sim_parameters)::Function
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        number_of_pots=parameters.number_of_pots
        driving_force=parameters.driving_force
        eta=parameters.rates_eta
        function rates(p::Int64,s::Int64)::Float64
            return 1 / tau * exp(beta / 2 * (interaction_strength / number_of_pots *
            (1 - s * eta) * (p+1) + s * driving_force))
        end
        return rates

    end 

    function glauber_rates(parameters::sim_parameters) :: Function
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        driving_force=parameters.driving_force
        eta=parameters.rates_eta
        function rates(p::Int64,s::Int64) :: Float64
            return 2 / tau * exp(s * beta / 2 * driving_force)/(exp(-interaction_strength/parameters.number_of_pots *
             p+1)+1)
        end
        return rates

    end 



end
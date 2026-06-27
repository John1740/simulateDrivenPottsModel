
module rates
   
    include("sim_data.jl")

    import ..sim_data: sim_parameters

    export get_rates_function

    function get_rates_function(parameters::sim_parameters):: Function

        
        rates =eval(Symbol(parameters.rates*"_rates"))(parameters)
        return rates
    end


    function arrhenius_rates(parameters::sim_parameters) :: Function
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        driving_force=parameters.driving_force
        eta=parameters.rates_eta
        function rates(p::Float64,s::Int64) :: Float64
            return 1 / tau * exp(beta / 2 * (interaction_strength *
            (1 - s * eta) * (p) + s * driving_force))
        end
        return rates

    end 

    function glauber_rates(parameters::sim_parameters) :: Function
        tau=parameters.tau
        beta=parameters.beta
        interaction_strength=parameters.interaction_strength
        driving_force=parameters.driving_force
        eta=parameters.rates_eta
        function rates(p::Float64,s::Int64) :: Float64
            return 2 / tau * exp(s * beta / 2 * driving_force)/(exp(-interaction_strength *
             p)+1)
        end
        return rates

    end 



end
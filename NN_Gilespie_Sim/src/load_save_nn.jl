module load_save_nn

    using JSON
    using JLD2

    include("sim_data_nn.jl")
    import ..sim_data_nn: sim_parameters, observable_number_of_saves, observable_data


    export load_simulation_parameters,load_observable_number_of_saves

    function load_simulation_parameters(json_file::String)::sim_parameters
        # Load JSON data from file
        version="NN_2.0"
        json_data = JSON.parse(read(json_file, String))
    
        # Default values
        default_values = Dict(
            "tau" => 1.0,
            "beta" => 1.0,
            "driving_force" => 0.0,
            "interaction_strength" => 1.0,
            "number_of_states" => 3,
            "number_of_pots" => 1,
            "dimensions" => [1,1],
            "total_number_of_episodes" => 1,
            "total_simulation_time" => 10.0,
            "skip_simulation_time" => 0.0,
            "lattice_type" => "cubic",
            "start_position" => "center",
            "boundary_conditions" => "periodic",
            "save_location" => "",
            "rates" => "arrhenius",
            "rates_eta" => 0,
            "seed" => 0,
            "version" => "0"
        )
    
        # Merge default values with JSON data
        merged_data = merge(default_values, json_data[1])

        allowed_keys = Set(fieldnames(sim_parameters))
        extra_keys = setdiff(Set(Symbol.(keys(merged_data))), allowed_keys)
    
        if !isempty(extra_keys)
            error("Error: Unexpected key(s) in JSON data: $(join(extra_keys, ", "))")
        end
    
        # Create the struct with the merged values
        return sim_parameters(
            merged_data["tau"],
            merged_data["beta"],
            merged_data["driving_force"],
            merged_data["interaction_strength"],
            merged_data["number_of_states"],
            prod(merged_data["dimensions"]),
            merged_data["dimensions"],
            merged_data["total_number_of_episodes"],
            merged_data["total_simulation_time"],
            merged_data["skip_simulation_time"],
            merged_data["lattice_type"],
            merged_data["start_position"],
            merged_data["boundary_conditions"],
            merged_data["save_location"],
            merged_data["rates"],
            merged_data["rates_eta"],
            merged_data["seed"],
            version
        )
    end

    function load_observable_number_of_saves(json_file::String) ::observable_number_of_saves
        # Load JSON data from file
        json_data = JSON.parse(read(json_file, String))
    
        # Default values
        default_values = Dict(
            "potts_state_array" => 0,
            "occupation_vector" => 0,
            "fourier_modes" => 0,
            "angle_change" => 0,
            "fourier_abs_avg" => 0,
            "entropy" => 0,
            "magnetization" => 0
        )
    
        # Merge default values with JSON data
        merged_data = merge(default_values, json_data[2])

        allowed_keys = Set(fieldnames(observable_number_of_saves))
        extra_keys = setdiff(Set(Symbol.(keys(merged_data))), allowed_keys)
    
        if !isempty(extra_keys)
            error("Error: Unexpected key(s) in JSON data: $(join(extra_keys, ", "))")
        end
    
        # Create the struct with the merged values
        return observable_number_of_saves(
            merged_data["potts_state_array"],
            merged_data["occupation_vector"],
            merged_data["fourier_modes"],
            merged_data["angle_change"],
            merged_data["fourier_abs_avg"],
            merged_data["entropy"],
            merged_data["magnetization"]
        )
    end

    function save_all(parameters::sim_parameters,observable_saves::observable_number_of_saves,observable_array::observable_data)#::Nothing
        if observable_saves.potts_state_array>0
            JLD2.save("potts_state_array.jld2","potts_state_array",observable_array.potts_state_array)

        end

        if observable_saves.occupation_vector>0
            JLD2.save("occupation_vector.jld2","occupation_vector",observable_array.occupation_vector)

        end

        if observable_saves.fourier_modes>0
            JLD2.save("fourier_modes.jld2","fourier_modes",observable_array.fourier_modes)

        end
        
        if observable_saves.angle_change>0
            JLD2.save("angle_change.jld2","angle_change",observable_array.angle_change)

        end

        if observable_saves.magnetization>0
            JLD2.save("magnetization.jld2","magnetization",observable_array.magnetization)

        end

        if observable_saves.fourier_abs_avg>0
            JLD2.save("fourier_abs_avg.jld2","fourier_abs_avg",observable_array.fourier_abs_avg)
            JLD2.save("fourier_abs_avg_square.jld2","fourier_abs_avg_square",observable_array.fourier_abs_avg_square)
            JLD2.save("fourier_abs_avg_quartic.jld2","fourier_abs_avg_quartic",observable_array.fourier_abs_avg_quartic)
        end

        if observable_saves.entropy>0
            JLD2.save("work.jld2","work",observable_array.work)
            JLD2.save("inflow.jld2","inflow",observable_array.inflow)
            JLD2.save("activity.jld2","activity",observable_array.activity)
        end

        json_data = JSON.json(parameters)

        open("simulation_parameters.json", "w") do file
            write(file, json_data)

        end
        
        return nothing
    end

end
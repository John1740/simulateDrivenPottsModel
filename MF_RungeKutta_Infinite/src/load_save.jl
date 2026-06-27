module load_save

    using JSON
    using JLD2

    include("sim_data.jl")
    import ..sim_data: sim_parameters, observable_number_of_saves, observable_data


    export load_simulation_parameters,load_observable_number_of_saves

    function load_simulation_parameters(json_file::String) :: sim_parameters
        # Load JSON data from file
        json_data = JSON.parse(read(json_file, String))
    
        # Default values
        default_values = Dict(
            "tau" => 1.0,
            "beta" => 1.0,
            "driving_force" => 0.0,
            "interaction_strength" => 1.0,
            "number_of_states" => 1,
            "total_number_of_episodes" => 1,
            "total_simulation_time" => 10.0,
            "skip_simulation_time" => 0.0,
            "dt" => 0.01,
            "save_location" => "",
            "rates" => "arrhenius",
            "rates_eta" => 0,
            "seed" => 0
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
            merged_data["total_number_of_episodes"],
            merged_data["total_simulation_time"],
            merged_data["skip_simulation_time"],
            merged_data["dt"],
            merged_data["save_location"],
            merged_data["rates"],
            merged_data["rates_eta"],
            merged_data["seed"]
        )
    end

    function load_observable_number_of_saves(json_file::String) :: observable_number_of_saves
        # Load JSON data from file
        json_data = JSON.parse(read(json_file, String))
    
        # Default values
        default_values = Dict(
            "occupation_vector" => 0,
            "fourier_modes" => 0
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
            merged_data["occupation_vector"],
            merged_data["fourier_modes"]
        )
    end

    function save_all(parameters::sim_parameters,observable_saves::observable_number_of_saves,observable_array::observable_data)
        if observable_saves.occupation_vector>0
            JLD2.save("occupation_vector.jld2","occupation_vector",observable_array.occupation_vector)

        end

        if observable_saves.fourier_modes>0
            JLD2.save("fourier_modes.jld2","fourier_modes",observable_array.fourier_modes)

        end
        
        json_data = JSON.json(parameters)

        open("simulation_parameters.json", "w") do file
            write(file, json_data)

        end
        
        return Nothing
    end

end
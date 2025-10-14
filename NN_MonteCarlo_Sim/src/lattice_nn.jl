module lattice_nn

    using JLD2
   
    include("sim_data_nn.jl")

    import ..sim_data_nn: sim_parameters, sim_state

    export initialize_lattice

    function initialize_lattice(parameters::sim_parameters,states::sim_state)#::Nothing

        allowed_boundry_conditions=["wall","periodic"]
        if parameters.boundry_conditions ∉ allowed_boundry_conditions
            error("Boundry conditon not known/implemented.")
        end
        
        allowed_lattice_types=["cubic","mean_field","hexagonal"]
        if parameters.lattice_type ∉ allowed_lattice_types
            error("Lattice type not known/implemented.")
        end

        eval(Symbol(parameters.lattice_type*"_lattice"))(parameters,states)
        return nothing
        
    end


    function mean_field_lattice(parameters::sim_parameters,states::sim_state)
        println("mean field is only for testing!(very Slow)")
        neigbourlist=zeros(Int64,parameters.number_of_pots,parameters.number_of_pots-1)

        for j in 1:parameters.number_of_pots

            neigbourlist[j,:]=[i for i in 1:parameters.number_of_pots if i != j]
        end
        states.neigbour_list=neigbourlist
        states.max_neigbours=parameters.number_of_pots-1

        return nothing
    end



    function cubic_lattice(parameters::sim_parameters,states::sim_state)#::Nothing
        
        dimension=length(parameters.dimensions)
        max_nn=dimension*2
        neigbourlist=zeros(Int64,parameters.number_of_pots,max_nn)
        coordinates=zeros(Int64,parameters.number_of_pots,dimension)
        dr=1


        function get_index(pos::Vector{Int64})
            pos_periodic=pos
            for i in 1:length(pos)
                if pos[i] <1 || pos[i] >parameters.dimensions[i] 
                    if parameters.boundry_conditions=="wall"
                        return 0
                    elseif parameters.boundry_conditions=="periodic"
                        pos_periodic[i]= mod(pos[i]-1,parameters.dimensions[i])+1
                    end
                end
            end

            index=pos[1]
            for i in 2:length(pos)
                index+=(pos[i]-1)*prod(parameters.dimensions[1:i-1])
            end
            return index
        end

        function get_coordinates(pos::Vector{Int64})

            return dr*pos
        end



        position=ones(Int64,length(parameters.dimensions))
        for i in 1:parameters.number_of_pots
            pot_index=get_index(position)
            for dims in 1:dimension
                dv=zeros(Int64,dimension)
                dv[dims]=1
                neigbourlist[pot_index,dims*2-1]=get_index(position+dv)
                neigbourlist[pot_index,dims*2]=get_index(position-dv)
            end
            
            coordinates[pot_index,:]=get_coordinates(position)

            dim=1
            while dim<=dimension
                position[dim]+=1
                if position[dim]<=parameters.dimensions[dim]
                    break
                else
                    position[dim]=1
                    dim+=1
                end
                
            end
            
        end

        @save "lattice.jld2" coordinates neigbourlist
        println("latice saved")


        states.neigbour_list=neigbourlist
        states.max_neigbours=max_nn

        return nothing
        
    end

    function hexagonal_lattice(parameters::sim_parameters,states::sim_state)#::Nothing
        
        dimension=length(parameters.dimensions)
        
        if any(x -> x % 2 != 0, parameters.dimensions) 
            error("Hexagonal lattice needs even amout in every dimension.")
        end

        if dimension ∉ [2,3] 
            error("This dimension of hcp is not implemented.")
        end

        max_nn=(dimension-1)*6
        neigbourlist=zeros(Int64,parameters.number_of_pots,max_nn)
        coordinates=zeros(Float64,parameters.number_of_pots,dimension)
        dx=1
        dy=dx*sqrt(3)/2
        dz=dx*sqrt(3)/2

        if dimension==2

            function get_index(pos::Vector{Int64})
                pos_periodic=pos
                for i in 1:length(pos)
                    if pos[i] <1 || pos[i] >parameters.dimensions[i] 
                        if parameters.boundry_conditions=="wall"
                            return 0
                        elseif parameters.boundry_conditions=="periodic"
                            pos_periodic[i]= mod(pos[i]-1,parameters.dimensions[i])+1
                        end
                    end
                end

                index=pos[1]
                for i in 2:length(pos)
                    index+=(pos[i]-1)*prod(parameters.dimensions[1:i-1])
                end
                return index
            end

            function get_coordinates(pos::Vector{Int64})
                mult=sqrt(3)/2*ones(dimension)
                mult[1]=dx
                coords=pos.*mult
                if iseven(pos[2])
                    coords[1]+=dx/2 
                    #coords[2]-=dy/2 
                end
                return coords
            end



            position=zeros(Int64,length(parameters.dimensions))
            position[1]=1 
            for i in 1:parameters.number_of_pots
                pot_index=get_index(position)


                neigbourlist[pot_index,1]=get_index(position+[1,0])
                neigbourlist[pot_index,2]=get_index(position+[-1,0])

                neigbourlist[pot_index,3]=get_index(position+[0,1])
                neigbourlist[pot_index,4]=get_index(position+[0,-1])

                if iseven(position[2])
                    di=+1
                else
                    di=-1
                end
                neigbourlist[pot_index,5]=get_index(position+[di,1])
                neigbourlist[pot_index,6]=get_index(position+[di,-1])

                coordinates[pot_index,:]=get_coordinates(position)

                dim=1
                while dim<=dimension
                    position[dim]+=1
                    if position[dim]<=parameters.dimensions[dim]
                        break
                    else
                        position[dim]=1
                        dim+=1
                    end
                    
                end

                
            end

        elseif dimension==3
            function get_index3(pos::Vector{Int64})
                pos_periodic=pos
                for i in 1:length(pos)
                    if pos[i] <1 || pos[i] >parameters.dimensions[i] 
                        if parameters.boundry_conditions=="wall"
                            return 0
                        elseif parameters.boundry_conditions=="periodic"
                            pos_periodic[i]= mod(pos[i]-1,parameters.dimensions[i])+1
                        end
                    end
                end

                index=pos[1]
                for i in 2:length(pos)
                    index+=(pos[i]-1)*prod(parameters.dimensions[1:i-1])
                end
                return index
            end
            function get_coordinates3(pos::Vector{Int64})
                mult=sqrt(3)/2*ones(dimension)
                mult[1]=dx
                coords=pos.*mult
                if iseven(pos[2])
                    coords[1]+=dx/2 
                end
                if iseven(pos[3])
                    coords[1]+=dx/2
                    coords[2]+=dx/2/sqrt(3)
                end
                return coords
            end



            position=zeros(Int64,length(parameters.dimensions))
            position[1]=1 
            #position=ones(Int64,length(parameters.dimensions))
            for i in 1:parameters.number_of_pots

                pot_index=get_index3(position)


                neigbourlist[pot_index,1]=get_index3(position+[1,0,0])
                neigbourlist[pot_index,2]=get_index3(position+[-1,0,0])

                neigbourlist[pot_index,3]=get_index3(position+[0,1,0])
                neigbourlist[pot_index,4]=get_index3(position+[0,-1,0])

                if iseven(position[2])
                    di2=1
                else
                    di2=0
                end
                neigbourlist[pot_index,5]=get_index3(position+[2*di2-1,1,0])
                neigbourlist[pot_index,6]=get_index3(position+[2*di2-1,-1,0])

                if iseven(position[3])
                    di3=1
                else
                    di3=0
                end

                neigbourlist[pot_index,7]=get_index3(position+[0+di3,0,1])
                neigbourlist[pot_index,8]=get_index3(position+[di2-1+di3,-1+2*di3,1])
                neigbourlist[pot_index,9]=get_index3(position+[-1+di3,0,1])

                neigbourlist[pot_index,10]=get_index3(position+[0+di3,0,-1])
                neigbourlist[pot_index,11]=get_index3(position+[di2-1+di3,-1+2*di3,-1])
                neigbourlist[pot_index,12]=get_index3(position+[-1+di3,0,-1])



                coordinates[pot_index,:]=get_coordinates3(position)

                dim=1
                while dim<=dimension
                    position[dim]+=1
                    if position[dim]<=parameters.dimensions[dim]
                        break
                    else
                        position[dim]=1
                        dim+=1
                    end
                    
                end

                
            end
        end

        @save "lattice.jld2" coordinates neigbourlist
        #println("latice saved")

        #println(neigbourlist)
        states.neigbour_list=neigbourlist
        states.max_neigbours=max_nn

        return nothing
        
    end

end
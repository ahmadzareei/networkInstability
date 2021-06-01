using MATLAB
using Notifier
using Plots

data = []
#number of runs
Ns = 4
N_list = collect(1:0.5:5)
M_list = collect(0:0.25:2)
T = 6001
#record step
record_step = 30#record step
num_record_per_sim = round(Int,T/record_step)+1

# seed_list = Array(1:12)
seed_list  = collect(1:10)

#maximum erosion rate
alpha = 0.01*20
normalized_Q = 0.
eigRecord = false
recordPoints = collect(1.5:0.5:2.5)
# using ProfileView
Nx,Ny,Nz = 50,50,12
working_mode = "erosion"
graphDict = Dict()
for i in 1:length(seed_list)
    graphDict[i] = initializer_jl(seed_list[i],Nx,Ny,R_min=1,R_max=14,randomness = 0.2)
end
for n in N_list
    for m in M_list
        data_folder = string("E:/TempCode/MatlabFlow/matData/",working_mode,"/DelaunayNet/",
        Nx,"by",Ny,"T",T,"d0.2","/","N",round(n,digits=2),"/","M",round(m,digits=2),"/","a",round(alpha,digits=2),"/")
        # data_folder = string("/matData/erosion/GridNet/sameLength/",Nx,"by",Ny,"T",T,"/")
        # data_folder = string("/matData/erosion/DelaunayNet/",Nx,"by",Ny,"T",T,"/")
        # if !ispath(data_folder)
        mkpath(data_folder)
    end
    for m in M_list
        #Simulation iterations
        #data folder
        data_folder = string("E:/TempCode/MatlabFlow/matData/",working_mode,"/DelaunayNet/",
        Nx,"by",Ny,"T",T,"d0.2","/","N",round(n,digits=2),"/","M",round(m,digits=2),"/","a",round(alpha,digits=2),"/")
        # data_folder = string("/matData/erosion/GridNet/sameLength/",Nx,"by",Ny,"T",T,"/")
        # data_folder = string("/matData/erosion/DelaunayNet/",Nx,"by",Ny,"T",T,"/")
        # if !ispath(data_folder)
        # mkpath(data_folder)
        # end
        #order parameters
        Data = zeros(length(seed_list),4,num_record_per_sim)

        Threads.@threads for i in 1:length(seed_list)
            # for i in 1:length(seed_list)
            rnd_seed = seed_list[i]
            # if n == 1
            #     T1 = 3001
            #     data_folder1 = string("E:/TempCode/MatlabFlow/matData/",working_mode,"/DelaunayNet/",
            #     Nx,"by",Ny,"T",T1,"/","N",round(n,digits=2),"/")
            #     G,locDict,lengthDict,radDict,knownPs,groundNodes,edgeIndex = initializer_jl(rnd_seed,Nx,Ny,R_min=1,R_max=14)
            #     d = read_matfile(string(data_folder1,"matLargeDataS",rnd_seed,".mat"))
            #     L_array = [lengthDict[e] for e in collect(edges(G))]
            #     R_array = jarray(d["R_t"])[end,:]
            #     C_array = R_array.^4 ./L_array
            #     P,eigMax = U_solver(G,groundNodes,knownPs,C_array,eigRecord=eigRecord) #Pressure
            #     Q = I_from_U(G,P,C_array) # Volume flow
            #     if normalized_Q != 0.
            #         outflow = sum(abs.(Q[ground_edges]))
            #         Q = normalized_Q/outflow.*Q
            #         P = normalized_Q/outflow.*P
            #     end
            #     num_record_per_sim1 =round(Int,T1/record_step)+1
            #     R_array_t = zeros(num_record_per_sim,ne(G))
            #     outflow_t = zeros(num_record_per_sim)
            #     weightP_t = zeros(num_record_per_sim,ne(G))
            #     R_array_t[1:num_record_per_sim1,:] = jarray(d["R_t"])
            #     outflow_t[1:num_record_per_sim1] = jarray(d["Outflow_t"])
            #     weightP_t[1:num_record_per_sim1,:] = jarray(d["WeightP_t"])
            #
            #     matFileLoc = string(data_folder,"matLargeDataS",rnd_seed,".mat")
            #     write_matfile(matFileLoc;
            #     R_t = R_array_t,
            #     Outflow_t = outflow_t,
            #     WeightP_t = weightP_t)
            #
            #     @time network_erosion!(data_folder,G,groundNodes,knownPs,normalized_Q, P, Q, R_array, T-T1, n , alpha ,
            #     L_array = L_array, matWrite = true,matFile = matFileLoc,
            #     edgeIndex = edgeIndex,record_step = record_step,record_data = nothing,
            #     record_start = num_record_per_sim1,
            #     eigRecord=eigRecord , breakpoint=true)
            # else
            # G,locDict,lengthDict,radDict,knownPs,groundNodes,edgeIndex = initializer_grid_jl(rnd_seed,Nx,Ny,R_min=1,R_max=14,sameLength=true)
            G,locDict,lengthDict,radDict,knownPs,groundNodes,edgeIndex = graphDict[i]
            # mat_G = mxarray(Matrix(adjacency_matrix(G)))
            sources = zeros(ne(G))
            targets = zeros(ne(G))
            k = 1
            for edge in edges(G)
                sources[k] = src(edge)
                targets[k] = dst(edge)
                k += 1
            end
            write_matfile(string(data_folder,"ST",rnd_seed,".mat"); s = sources, t=targets)
            pointPos = zeros(nv(G),2)
            for j in 1:nv(G)
                pointPos[j,1] = locDict[j][1]
                pointPos[j,2] = locDict[j][2]
                # pointPos[j,3] = locDict[j][3]
            end
            pointPos = mxarray(pointPos)
            ground_edges  = right_edges(G,groundNodes)
            #run the simulation and save data
            # num of iterations
            #total flow (if set to 0. then it is constant pressure case)
            C = Dict(e=>radDict[e]^4/lengthDict[e] for e in edges(G)) #Conducctance
            R_array = [radDict[e] for e in collect(edges(G))]
            L_array = [lengthDict[e] for e in collect(edges(G))]
            C_array = [C[edge] for edge in collect(edges(G))]
            P,eigMax = U_solver(G,groundNodes,knownPs,C_array,eigRecord=eigRecord) #Pressure
            Q = I_from_U(G,P,C_array) # Volume flow
            if normalized_Q != 0.
                outflow = sum(abs.(Q[ground_edges]))
                Q = normalized_Q/outflow.*Q
                P = normalized_Q/outflow.*P
            end
            write_matfile(string(data_folder,"configArrayS",rnd_seed,".mat"); posArray = pointPos, lenArray = L_array )
            # path = string(eval(@__DIR__),data_folder)
            # if !isdir(path)
            #     mkdir(path)
            # end

            record_data = zeros(4,num_record_per_sim)
            # mxcall(:initializer_jl,1,rnd_seed)
            #load network from initializer.mat


            # write_matfile(string(eval(@__DIR__),data_folder,"N_",round(n,digits=2),"_0",".mat"),
            # diameterP = R_array,
            # outflow = sum(abs.(Q[ground_edges])),
            # weightP = Q,
            # )

            R_array_t = zeros(num_record_per_sim,ne(G))
            outflow_t = zeros(num_record_per_sim)
            weightP_t = zeros(num_record_per_sim,ne(G))
            R_array_t[1,:] = copy(R_array)
            outflow_t[1] = sum(abs.(Q[ground_edges]))
            weightP_t[1,:] = Q
            matFileLoc = string(data_folder,"matLargeDataS",rnd_seed,".mat")
            write_matfile(matFileLoc;
            R_t = R_array_t,
            Outflow_t = outflow_t,
            WeightP_t = weightP_t)

            Ne = ne(G)
            record_data[1,1] = 1/(Ne-1)*(Ne- sum(Q.^2)^2/sum(Q.^4))
            record_data[2,1] = sum(R_array)/Ne
            record_data[3,1] = sum(abs.(Q[ground_edges]))/maximum(P)
            record_data[4,1] = eigMax
            #@time network_erosion!(G,groundNodes,knownPs,Dict(),normalized_Q, P, Q, radDict, T, n , c ,L_dict = lengthDict, matWrite = true,matWritePara =string("N_",n),edgeIndex = edgeIndex,record_step = record_step)
            # @time network_clog!(data_folder,G,groundNodes,knownPs,normalized_Q, P, Q, R_array, T, n , alpha ,
            # L_array = L_array, matWrite = true,matWritePara =string("N_",round(n,digits=2)),
            # edgeIndex = edgeIndex,record_step = record_step,record_data = record_data)
            if working_mode == "clog"
                @time network_clog!(data_folder,G,groundNodes,knownPs,normalized_Q, P, Q, R_array, T, n , alpha ,
                L_array = L_array, matWrite = true,matFile = matFileLoc,
                edgeIndex = edgeIndex,record_step = record_step,record_data = record_data,eigRecord=eigRecord)
            else
                @time network_erosion!(data_folder,G,groundNodes,knownPs,normalized_Q, P, Q, R_array, T, n , alpha ,
                m = m, L_array = L_array, matWrite = true,matFile = matFileLoc,
                edgeIndex = edgeIndex,record_step = record_step,record_data = record_data,
                eigRecord=eigRecord ,  recordPoints = recordPoints)
            end
            # d1, d2, d3, R1, R2, R3 = mxcall(:shortestPathGuess,6,"N_0_100.mat",string("RS_",i,"_N_",0))
            # append!(data,[d1,d2,d3,R1,R2,R3,sum(1 ./collect(values(C)))/length(C)]');
            Data[i,:,:] = record_data
        end
        write_matfile(string(data_folder,"A_Matlab_Data",".mat"),Data = Data,normQ = normalized_Q, seed = seed_list)
    end
end
# plot(N_list,Data[:,1,1:2:end],st = :scatter,marker_zs = Data[:,2,1:2:end],label=nothing,
# m= (:lighttest,),colorbar_title="<R>",xaxis = "N",yaxis="L",title=working_mode)
notify("Task completed",sound = true)

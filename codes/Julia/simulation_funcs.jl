using SparseArrays
using LightGraphs
using LinearAlgebra
using MAT
using Random
using VoronoiDelaunay
using Plots
using QHull
include("QDelaunay.jl")
include("VoronoiPy.jl")

### plotting functions for weighted graph w or w/o pressure distribution given
struct GraphWithWeights
    x
    y
    lws
end

@recipe function f(g::GraphWithWeights)
    seriestype  :=  :path
    legend --> false
    grid --> false
    axis --> false
    linewidth --> g.lws
    linealpha --> g.lws/3.5
    linecolor   --> :black
#     linecolor --> map(lw ->RGB(cgrad(:heat)[lw]),g.lws./3)
    html_output_format --> :png
    x,y = g.x,g.y
end

struct GraphWithPressure
    x
    y
    Us
end

struct Graph3dWeights
    x
    y
    z
    lws
end
@recipe function f(g::Graph3dWeights)
    seriestype  :=  :path
    legend --> false
    grid --> true
    axis --> true
    label --> true
    camera --> (45,45)
    linewidth --> g.lws
    linealpha --> g.lws/3.5
    linecolor   --> :black
#     linecolor --> map(lw ->RGB(cgrad(:heat)[lw]),g.lws./3)
    html_output_format --> :png
    x,y,z = g.x,g.y,g.z
end

@recipe function f(F::GraphWithPressure,NBins = nothing)
    legend --> false
    axis --> false
    grid --> false
    html_output_format --> :png
    seriestype --> :heatmap
    seriesalpha --> 0.5
    seriescolor -->[:blue,:white,:red]
    colorbar --> nothing
    if NBins != nothing
        N_bins = NBins
    else
        N_bins = floor(Int,sqrt(length(F.x))/1.5)
    end
    h1 = fit(Histogram, (F.x,F.y), nbins=N_bins)
    h2 = fit(Histogram, (F.x,F.y), StatsBase.weights(F.Us),nbins=N_bins)
    # x --> h1.edges[1]
    # y --> h1.edges[2]
    dx = h1.edges[1][2] - h1.edges[1][1]
    dy = h1.edges[2][2] - h1.edges[2][1]
    data --> (h1.edges[1], h1.edges[2],(h2.weights./h1.weights)')
    # print(h2.edges)
end


@inbounds function PlotGraphWithLocPressure(G::Graph, Loc::Dict, I=nothing, U= nothing; projection = nothing)
    #return x,y,lw that can reproduce the network
    ne = LightGraphs.ne(G)
    d = incidence_matrix(G,oriented = true)
    x = Array{Float64,2}(undef,2,ne)
    y = Array{Float64,2}(undef,2,ne)
    lws = ones(Float64,1,ne)
    if projection != nothing
        k,j = projection[1],projection[1]
    else
        k,j = 1,2
    end
    i = 1
    if I != nothing
        width_max = maximum(abs.(I))
        for i in 1:ne
            u,v = findnz(d[:,i])[1]
            x[1,i]=Loc[u][k]
            x[2,i]=Loc[v][k]
            y[1,i]=Loc[u][j]
            y[2,i]=Loc[v][j]
            lws[1,i] = max(0,abs(I[i])/width_max*3-1e-4)
        end
        # for e in edges(G)
        #     u,v = src(e), dst(e)
        #
        #     i += 1
        # end
    else
        for e in edges(G)
            u,v = src(e), dst(e)
            x[1,i]=Loc[u][k]
            x[2,i]=Loc[v][k]
            y[1,i]=Loc[u][j]
            y[2,i]=Loc[v][j]
            i += 1
        end
    end
    Nv = nv(G)
    xv = ones(Float64,Nv)
    yv = ones(Float64,Nv)
    Us = ones(Float64,Nv)
    i = 1
    if U != nothing
        U_max = maximum(values(U))
        for v in vertices(G)
            xv[i]=Loc[v][k]
            yv[i]=Loc[v][j]
            Us[i] = max(0,abs(U[v])/U_max-1e-6)
            i += 1
        end
    else
        for v in vertices(G)
            xv[i]=Loc[v][k]
            yv[i]=Loc[v][j]
            i += 1
        end
    end
    GraphWithWeights(x, y, lws),GraphWithPressure(xv,yv,Us)
end

@inbounds function PlotGraph3dWithWeights(G::Graph, Loc::Dict, I=nothing)
    #return x,y,lw that can reproduce the network
    ne = LightGraphs.ne(G)
    d = incidence_matrix(G,oriented = true)
    x = Array{Float64,2}(undef,2,ne)
    y = Array{Float64,2}(undef,2,ne)
    z = Array{Float64,2}(undef,2,ne)
    lws = ones(Float64,1,ne)
    i = 1
    if I != nothing
        width_max = maximum(abs.(I))
        for i in 1:ne
            u,v = findnz(d[:,i])[1]
            x[1,i]=Loc[u][1]
            x[2,i]=Loc[v][1]
            y[1,i]=Loc[u][2]
            y[2,i]=Loc[v][2]
            z[1,i]=Loc[u][3]
            z[2,i]=Loc[v][3]
            lws[1,i] = max(0,abs(I[i])/width_max*3-1e-4)
        end
        # for e in edges(G)
        #     u,v = src(e), dst(e)
        #
        #     i += 1
        # end
    else
        for e in edges(G)
            u,v = src(e), dst(e)
            x[1,i]=Loc[u][1]
            x[2,i]=Loc[v][1]
            y[1,i]=Loc[u][2]
            y[2,i]=Loc[v][2]
            z[1,i]=Loc[u][3]
            z[2,i]=Loc[v][3]
            i += 1
        end
    end
    Graph3dWeights(x, y, z, lws)
end

### random 2D delaunay triangulation
function random_delaunay(N ; pos=false,random_seed=0)

    Random.seed!(random_seed)
    width = (max_coord - min_coord)*0.99
    tess = DelaunayTessellation()
    i_to_pos = Dict([(i,(min_coord+rand()*width+0.01,min_coord+rand()*width+0.01) ) for i in 1:N])
    n_left = round(Int,sqrt(N))
    n_right,n_top,n_bottom = n_left,n_left,n_left
    n_total = N + 4*n_left
    m = n_left+n_right
    for i in 1:n_left
        i_to_pos[N+i] = (min_coord+0.01, min_coord+width*(i-0.5)/n_left)
    end
    for i in 1:n_right
        i_to_pos[N+n_left+i] = (max_coord-0.01, min_coord+width*(i-0.5)/n_right)
    end
    for i in 1:n_top
        i_to_pos[N+m+i] = (min_coord+width*(i-0.5)/n_top, max_coord-0.01)
    end
    for i in 1:n_bottom
        i_to_pos[N+m+n_top+i] = (min_coord+width*(i-0.5)/n_bottom, min_coord+0.01)
    end
    pos_to_i = Dict(value=>key for (key,value) in i_to_pos)
    a = Point2D[Point(i_to_pos[i][1], i_to_pos[i][2]) for i in 1:n_total]
    push!(tess,a)

    A = spzeros(n_total,n_total)
    for edge in delaunayedges(tess)
        pos_1 = (getx(geta(edge)),gety(geta(edge)))
        pos_2 = (getx(getb(edge)),gety(getb(edge)))
        if abs2(pos_1[1]-pos_2[1])+abs2(pos_1[2]-pos_2[2]) > 0.002^2
            A[pos_to_i[pos_1],pos_to_i[pos_2]] = 1
            A[pos_to_i[pos_2],pos_to_i[pos_1]] = 1
        end
    end
    if pos
        return Graph(A),copy(i_to_pos)
    else
        return Graph(A)
    end
end

function random_2d_graph(Nx,Ny;randomSeed=0, randomness = 0.2) # nx is used to refer networkx
  #first create a standard grid of nx*ny points in Dx*Dy region
    nodes = collect(1:Nx*Ny)
    n_total = Nx*Ny
    Dx = max_coord - min_coord
    Dy = max_coord - min_coord
    xs = min_coord .+ fld.(nodes .- 1,Ny) .*Dx/(Nx-1)
    ys = min_coord .+ rem.(nodes .- 1,Ny) *Dy/(Ny-1)
    #then the boundary nodes are clear
    leftNodes = collect(Ny÷4:3*Ny÷4)
    rightNodes = collect(Ny*(Nx-1)+Ny÷4:Nx*Ny-Ny÷4)
    # now add some random x and y displacement to the internal nodes (keep the boundary ones fixed)
    Random.seed!(randomSeed)
    dxs = randomness*Random.randn(Ny*(Nx-2))*Dx/(Nx-1) #randomness defines how random the nodes are moved
    dys = randomness*Random.randn(Ny*(Nx-2))*Dy/(Ny-1)
    xs[Ny+1:Ny*(Nx-1)] = xs[Ny+1:Ny*(Nx-1)] + dxs
    ys[Ny+1:Ny*(Nx-1)] = ys[Ny+1:Ny*(Nx-1)] + dys
    xs = max.( min.( xs,max_coord),min_coord)
    ys = max.( min.( ys,max_coord),min_coord)

    #make a delaunay triangulation of the point data

    tess = DelaunayTessellation()
    a = Point2D[Point(xs[i], ys[i]) for i in 1:n_total]
    pos_to_i = Dict(Point(xs[i], ys[i])=>i for i in 1:n_total)


    # create a set for edges that are indexes of the points
    # edges = set()
    # for each Delaunay triangle
    push!(tess,a)

    A = spzeros(n_total,n_total)
    for edge in delaunayedges(tess)
        pos_1 = (getx(geta(edge)),gety(geta(edge)))
        pos_2 = (getx(getb(edge)),gety(getb(edge)))
        if abs2(pos_1[1]-pos_2[1])+abs2(pos_1[2]-pos_2[2]) > 0.002^2
          A[pos_to_i[geta(edge)],pos_to_i[getb(edge)]] = 1
          A[pos_to_i[getb(edge)],pos_to_i[geta(edge)]] = 1
        end
    end
    xs = (xs .- min_coord)*Nx
    ys = (ys .- min_coord)*Ny
    # xs, ys = ys, xs
    points = Dict(i=>(xs[i],ys[i]) for i in 1:n_total)
    return Graph(A),points,leftNodes,rightNodes
end


function random_3d_graph(Nx,Ny,Nz;randomSeed=0, randomness = 0.2)
    #generate random 3d delaunay networks
    nodes = collect(1:Nx*Ny*Nz)
    n_total = Nx*Ny*Nz
    coords = zeros(n_total,3)
    i = 1
    for x in 1:Nx
        for y in 1:Ny
            for z in 1:Nz
                coords[i,:] = [x y z]
                i += 1
            end
        end
    end
    botNodes = collect(1:Nz*Ny)
    topNodes = collect(n_total-Nz*Ny+1:n_total)
    # now add some random x and y displacement to the internal nodes (keep the boundary ones fixed)
    Random.seed!(randomSeed)
    dxs = randomness*Random.randn(Nz*Ny*(Nx-2)) #randomness defines how random the nodes are moved
    dys = randomness*Random.randn(Nz*Ny*(Nx-2))
    dzs = randomness*Random.randn(Nz*Ny*(Nx-2))
    coords[Nz*Ny+1:n_total-Nz*Ny,1] = coords[Nz*Ny+1:n_total-Nz*Ny,1] + dxs
    coords[Nz*Ny+1:n_total-Nz*Ny,2] = coords[Nz*Ny+1:n_total-Nz*Ny,2] + dys
    coords[Nz*Ny+1:n_total-Nz*Ny,3] = coords[Nz*Ny+1:n_total-Nz*Ny,3] + dzs
    points = Dict(i=>coords[i,:] for i in 1:n_total)
    ch = QDelaunay.delaunay3d(coords)
    A = spzeros(n_total,n_total)

    #make a delaunay triangulation of the point data
    # points =[(xs[i],ys[i]) for i in 1:Nx*Ny]
    # posArray = (np.vstack((xs,ys)))
    Edges = Set()
    for simplex in ch.simplices
        push!(Edges,(simplex[1],simplex[2]))
        push!(Edges,(simplex[2],simplex[3]))
        push!(Edges,(simplex[3],simplex[4]))
        push!(Edges,(simplex[4],simplex[1]))
    end
    removed_Edges = Set()
    for edge in Edges
        coord1 = coords[edge[1],:]
        coord2 = coords[edge[2],:]
        if sum((coord1-coord2).^2)>5
            push!(removed_Edges,edge)
        end
    end
    setdiff!(Edges,removed_Edges)
    for edge in Edges
        A[edge[1],edge[2]] = 1
        A[edge[2],edge[1]] = 1
    end
    G = Graph(A)

    return G,points,botNodes,topNodes
    # end
end

function initializer_jl(rs,Nx,Ny;R_min=1,R_max=10,randomness = 0.2)
    #initialize 2d networks with uniformly random distribution of radii
    Random.seed!(rs)
    G, locDict,leftNodes,rightNodes = random_2d_graph(Nx,Ny,randomSeed=rs,randomness=randomness)
    edge_array = collect(edges(G))
    edgeIndex = Dict(i=>edge_array[i] for i in 1:ne(G))
    radDict = Dict(edge_array[i]=>Array(R_min.+(R_max-R_min)*Random.rand(ne(G)))[i] for i in 1:ne(G))
    lengthDict = Dict()
    for e in edge_array
        lengthDict[e] = sqrt((locDict[src(e)][1]-locDict[dst(e)][1])^2+(locDict[src(e)][2]-locDict[dst(e)][2])^2)
    end
    knownPs = Dict(i => 1 for i in leftNodes)
    groundNodes = rightNodes
    print("Network initialized of N = ", string(length(vertices(G))), " M = ",string(length(edges(G))),"\n")
    return G,locDict,lengthDict,radDict,knownPs,groundNodes,edgeIndex
end

function initializer_3d_jl(rs,Nx,Ny,Nz;R_min=1,R_max=10,randomness = 0.2)
    #initialize 3D Delaunay networks with uniformly random distribution of radii
    Random.seed!(rs)
    G, locDict,leftNodes,rightNodes = random_3d_graph(Nx,Ny,Nz,randomSeed=rs,randomness = randomness)
    edge_array = collect(edges(G))
    edgeIndex = Dict(i=>edge_array[i] for i in 1:ne(G))
    radDict = Dict(edge_array[i]=>Array(R_min.+(R_max-R_min)*Random.rand(ne(G)))[i] for i in 1:ne(G))
    lengthDict = Dict()
    for e in edge_array
        lengthDict[e] = sqrt(sum((locDict[src(e)][i]-locDict[dst(e)][i])^2 for i in 1:3))
    end
    knownPs = Dict(i => 1 for i in leftNodes)
    groundNodes = rightNodes
    print("Network initialized","\n")
    return G,locDict,lengthDict,radDict,knownPs,groundNodes,edgeIndex
end

function initializer_3d_vor_jl(rs,Nx,Ny,Nz;R_min=1,R_max=10,randomness = 0.2)
    #initialize 3D Voronoi networks with uniformly random distribution of radii
    Random.seed!(rs)
    adj_G,posArray,leftNodes,rightNodes =  VoronoiPy.Random3DNetwork(Nx,Ny,Nz,randomness,rs)
    locDict = Dict(i=>posArray[i,:] for i in 1:size(posArray)[1])
    G = Graph(adj_G)
    edge_array = collect(edges(G))
    edgeIndex = Dict(i=>edge_array[i] for i in 1:ne(G))
    radDict = Dict(edge_array[i]=>Array(R_min.+(R_max-R_min)*Random.rand(ne(G)))[i] for i in 1:ne(G))
    lengthDict = Dict()
    for e in edge_array
        lengthDict[e] = sqrt(sum((locDict[src(e)][i]-locDict[dst(e)][i])^2 for i in 1:3))
    end
    knownPs = Dict(i => 1 for i in leftNodes)
    groundNodes = rightNodes
    print("Network initialized of N = ", string(length(vertices(G))), " M = ",string(length(edges(G))),"\n")
    return G,locDict,lengthDict,radDict,knownPs,groundNodes,edgeIndex
end

function grid_graph(N_x,N_y;pos=false)
    #grid graph generator
    N = N_x*N_y
    xs = [fld(i,N_y) for i in 0:N-1]
    ys = [rem(i,N_y) for i in 0:N-1]
    A = zeros(N,N)
    for i in 1:N, j in 1:N
        abs(xs[i]-xs[j])+abs(ys[i]-ys[j])== 1 ? A[i,j]=1 : nothing
    end
    G = Graph(A)
    if pos
        return G,Dict(i=>(xs[i]/(N_x-1),ys[i]/(N_y-1)) for i in 1:N)
    else
        return G
    end
end

function initializer_grid_jl(rs,Nx,Ny;R_min=1,R_max=10,sameLength=false,branch=false)
    #initialize 2d grid networks with uniformly random distribution of radii
    Random.seed!(rs)
    # G, locDict = random_delaunay(N,pos=true,random_seed=rs)
    N = Nx*Ny
    G, locDict = grid_graph(Ny,Nx,pos=true)
    edge_array = collect(edges(G))
    edgeIndex = Dict(i=>edge_array[i] for i in 1:ne(G))
    radDict = Dict(edge_array[i]=>Array(R_min.+(R_max-R_min)*Random.rand(ne(G)))[i] for i in 1:ne(G))
    lengthDict = Dict()
    for e in edge_array
        # lengthDict[e] = sqrt((locDict[src(e)][1]-locDict[dst(e)][1])^2+(locDict[src(e)][2]-locDict[dst(e)][2])^2)

        #same length
        if sameLength
            lengthDict[e] = 1.
        else
            lengthDict[e] = sqrt((locDict[src(e)][1]-locDict[dst(e)][1])^2+(locDict[src(e)][2]-locDict[dst(e)][2])^2)
        end
    end
    n_left = Ny
    n_right = n_left
    if branch
        knownPs = Dict(n_left ÷ 2 => 1)
    else
        knownPs = Dict(i => 1 for i in 1:n_left)
    end
    groundNodes =collect(N-n_right+1:N)
    print("Network initialized")
    return G,locDict,lengthDict,radDict,knownPs,groundNodes,edgeIndex
end

function U_solver(G::Graph, groundNodes, knownUs::Dict, C; eigRecord=false)
    #solve Kirchoff's law for given ground nodes and given voltage boundary conditions
    #pure matrix implementation
    Es = collect(edges(G))
    d = incidence_matrix(G,oriented=true)
    if C isa Base.Dict
        #incase C is still given as a dictionary
        c = [C[(src(edge),dst(edge))] for edge in collect(edges(G)) ]
    else
        c = C
    end
    D = d*Diagonal(c)*transpose(d)
    n = nv(G)
    m = length(groundNodes)+length(knownUs)
    z = zeros(n+m)
    E = spzeros(m,n)
    j = 1
    for i in groundNodes[1:end]
        E[j,i] = 1
        j += 1
    end
    for i in keys(knownUs)
        E[j,i] = 1
        z[n+j] = knownUs[i]
        j += 1
    end
    Q = [D E';E spzeros(m,m)]
#     print(cond(Array(Q)))
    U = Q\z
    if eigRecord
        return U[1:n],cond(Array(Q),2)
    else
        return U[1:n],0
    end
end

function I_from_U(G,U,C)
    #solving current from given conductance and voltage
    #pure matrix implementation
    d = incidence_matrix(G,oriented=true)
    return Diagonal(C)*transpose(d)*U
end

function right_edges(G,groundNodes)
    #return the edge list connected to ground nodes
    i = 1
    right_edge_array = []
    for edge in collect(edges(G))
        if src(edge) in groundNodes || dst(edge) in groundNodes
            append!(right_edge_array,i)
        end
        i += 1
    end
    return right_edge_array
end

function network_clog!(data_folder,G,Ground_nodes,Known_Us,normalized_I,U,I,Rad,T,n,c;
    L_array = nothing, matWrite = false, matFile = nothing,
    edgeIndex = nothing, record_step = 20, record_data=nothing, eigRecord =false, expBound = false,Vc=1.0,dt=0.01)
    ground_edges  = right_edges(G,Ground_nodes)
    #clogging iterations
    #pure matrix implementation
    records = 1
    Ne = ne(G)
    left_nodes = collect(keys(Known_Us))
    if record_data != nothing
        K0 = record_data[3,1]
    else
        K0 = 1.
    end

    if matWrite
        d = read_matfile(matFile)
        R_array_t = jarray(d["R_t"])
        outflow_t = jarray(d["Outflow_t"])
        weightP_t = jarray(d["WeightP_t"])
    end
    for t in 1:T
        if expBound
            #c=dt
            dR = c*dt*abs.(I./(Rad.^n)).*exp.(-abs.(I)./(Rad.^2)./Vc)./L_array
        else
            max_delta = maximum(abs.(I./(Rad.^n)))
            dR = c/max_delta*abs.(I./(Rad.^n))
        end
        # Rad = max.(Rad .- dR, 1e-4)
        #break when radius is smaller than 1e-4
        Rad = (Rad .- dR)
        if minimum(Rad) < 1e-4
            print("One tube blocked, break")
            break
        end
        if L_array == nothing
            C = Rad.^4
        else
            C = Rad.^4 ./L_array
        end
        U,eigMax = try
            U_solver(G,Ground_nodes,Known_Us,C,eigRecord=eigRecord)
        catch y
            print("N=",n," breaks due to disconnect")
            break
        end
        I = I_from_U(G,U,C)
        if normalized_I != 0.
            outflow = sum(abs.(I[ground_edges]))
            K = outflow/U[left_nodes[1]]
            I = normalized_I/outflow.*I
            U = normalized_I/outflow.*U
        else
            outflow = sum(abs.(I[ground_edges]))
            K = outflow/U[left_nodes[1]]
        end
        if K < 1e-6*K0
            print("Reach 1e-6 min K")
            break
        end
        if t%record_step == 0
            records += 1
            if matWrite
                R_array_t[records,:] =  copy(Rad)
                outflow_t[records] = sum(abs.(I[ground_edges]))
                weightP_t[records,:] = I
            end
            if record_data != nothing
                #data 1 = order paramertes, data 2 = average R , data 3 = K
                record_data[1,records] = 1/(Ne-1)*(Ne- sum(I.^2)^2/sum(I.^4))
                record_data[2,records] = sum(Rad)/Ne
                record_data[3,records] = K
                record_data[4,records] = eigMax
            end
        end
    end
    if matWrite
        write_matfile(matFile;
        R_t = R_array_t,
        Outflow_t = outflow_t,
        WeightP_t = weightP_t)
    end
end

function network_erosion!(data_folder,G,Ground_nodes,Known_Us,normalized_I,U,I,Rad,T,n,c;
    m = 1,
    L_array = nothing, matWrite = false, matFile = nothing,
    edgeIndex = nothing, record_step = 20, record_data=nothing, record_start = 1,
    eigRecord = false, recordPoints = nothing)
    ground_edges  = right_edges(G,Ground_nodes)
    #erosion iterations
    #pure matrix implementation
    records = record_start
    Ne = ne(G)
    left_nodes = collect(keys(Known_Us))
    j = 2
    if record_data != nothing
        K0 = record_data[3,record_start]
    else
        K0 = 1.
    end

    if matWrite
        d = read_matfile(matFile)
        R_array_t = jarray(d["R_t"])
        outflow_t = jarray(d["Outflow_t"])
        weightP_t = jarray(d["WeightP_t"])
    end
    for t in 1:T
        max_delta = maximum(abs.(I).^m./(Rad.^n))
        dR = c/max_delta*abs.(I).^m./(Rad.^n)
        Rad = Rad .+ dR
        if L_array == nothing
            C = Rad.^4
        else
            C = Rad.^4 ./L_array
        end
        U,eigMax = try
            U_solver(G,Ground_nodes,Known_Us,C,eigRecord=eigRecord)
        catch y
            print("N=",n," breaks due to disconnect")
            break
        end
        I = I_from_U(G,U,C)
        if normalized_I != 0.
            outflow = sum(abs.(I[ground_edges]))
            I = normalized_I/outflow.*I
            U = normalized_I/outflow.*U
            K = 1/U[left_nodes[1]]
        else
            outflow = sum(abs.(I[ground_edges]))
            K = outflow/U[left_nodes[1]]
        end
        if K > 1e8*K0
            print("Reach 1e8 max K")
            break
        end
        if t%record_step == 0
            records += 1
            print("Simulate t = ",string(t),"\n")
            if record_data != nothing
                #data 1 = order paramertes, data 2 = average R , data 3 = K
                record_data[1,records] = 1/(Ne-1)*(Ne- sum(I.^2)^2/sum(I.^4))
                record_data[2,records] = sum(Rad)/Ne
                record_data[3,records] = K
                record_data[4,records] = eigMax
            end
            if recordPoints != nothing
                R_ave = sum(Rad)/Ne
                if R_ave > recordPoints[j]*7.5
                    print("Reach radius = ",string(R_ave))
                    if matWrite
                        R_array_t[j,:] =  copy(Rad)
                        outflow_t[j] = sum(abs.(I[ground_edges]))
                        weightP_t[j,:] = I
                    end
                    j += 1
                    if j>length(recordPoints)
                        break
                    end
                end
            end
        end
    end
    if matWrite
        write_matfile(matFile;
        R_t = R_array_t,
        Outflow_t = outflow_t,
        WeightP_t = weightP_t)
    end
end

function load_graph_from_mat(filename)
    #read MAT file and return Julia graph object
    vars = matread(filename)
    s = vars["s"]
    srcs = round.(Int,s)
    t = vars["t"]
    dsts = round.(Int,t)
    G = SimpleGraph(round(Int,vars["n_tot"]))
    radDict = Dict()
    edgeIndex = Dict()
    for i in 1:length(s)
        add_edge!(G,srcs[i],dsts[i])
        radDict[(srcs[i],dsts[i])] = vars["Diameters"][i]
        radDict[(dsts[i],srcs[i])] = vars["Diameters"][i]
        edgeIndex[i] = (srcs[i],dsts[i])
    end
    locDict = Dict(i=>(vars["pointX"][i],vars["pointY"][i]) for i in 1:round(Int,vars["n_tot"]))
    lengthDict = Dict()
    for e in edges(G)
        s,d = src(e),dst(e)
        lengthDict[(s,d)] = sqrt(abs2(locDict[s][1]-locDict[d][1])+abs2(locDict[s][2]-locDict[d][2]))
        lengthDict[(d,s)] = lengthDict[(s,d)]
    end
    knownUs = Dict(i=>(vars["P_left"]-vars["P_right"]) for i in round.(Int,vars["left_nodes"]))
    groundNodes = round.(Int,vars["right_nodes"])
    return G,locDict,lengthDict,radDict,knownUs,groundNodes,edgeIndex
end

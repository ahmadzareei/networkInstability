module QDelaunay

export Delaunay3D, delaunay3d, Voronoi3d, voronoi3d

using PyCall
const spatial = PyNULL()

function __init__()
    copy!(spatial, pyimport_conda("scipy.spatial", "scipy"))
end

mutable struct Delaunay3d{T<:Real}
    points::Matrix{T}
    # vertices::Vector{Int}
    simplices::Vector{Vector{Int}}
end
## helper for base-0 / base-1 difference
incone(x) = for i in 1:length(x)
x[i] += 1
end

mutable struct Voronoi3d{T<:Real}
end

function delaunay3d(x::Matrix{T}) where T<:Real
    py = spatial.Delaunay(x,incremental=true, qhull_options="QJ")
    points = convert(Matrix{T}, py."points")
    # vertices = convert(Vector{Int}, py."vertices")
    # incone(vertices)
    simplices = convert(Vector{Vector{Int}}, py."simplices")
    for simplex in simplices
        incone(simplex)
    end
    Delaunay3d(points, simplices,)
end

end

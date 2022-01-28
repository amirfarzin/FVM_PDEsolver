#include("HPDE_options.jl")
##########################################################

mutable struct FVMPDEProblem
    grid::FVMPDEGrid

    """
    Flux function F in x direction

    """
    F::Function

    """
    Flux function G in y direction
    !!! compat
        Will be added in future versions
    """
    G::Function

    """
    Flux function H in z direction
    !!! compat
        Will be added in future versions
    """
    H::Function

    """
    Source function
    !!! compat
        Will be added in future versions
    """
    S::Function

    nvars::Int8

    param::Any

    U0

    U_min
    U_max

    jac
    ode_problem_kwargs
    tspan
    t

    function FVMPDEProblem(grid::FVMPDEGrid, tspan, U0; kwargs...)
        if haskey(kwargs,:jacobian)
            jac = kwargs[:jacobian]
            kwargs = Dict([p for p in pairs(kwargs) if p[1] != :jacobian])
        else
            jac = true
        end
        zeroFunc(U) = zeros(size(U))
        F = kwargs[:F]
        @assert size(U0,1) == grid.nx
        if haskey(kwargs,:G)
            G = kwargs[:G]
            kwargs = Dict([p for p in pairs(kwargs) if p[1] != :G])
            @assert size(U0,2) == grid.ny
            if haskey(kwargs,:H)
                H = kwargs[:H]
                kwargs = Dict([p for p in pairs(kwargs) if p[1] != :H])
                nvars = size(U0,4)
                @assert grid.dimension == _3D
                @assert length(size(U0))-1 == 3
                @assert size(U0,3) == grid.nz
            else
                H = zeroFunc
                nvars = size(U0,3)
                @assert grid.dimension == _2D
                @assert length(size(U0))-1 == 2
            end
        else
            @assert !haskey(kwargs,:H)
            G = zeroFunc
            H = zeroFunc
            nvars = size(U0,2)
            @assert grid.dimension == _1D
            @assert length(size(U0))-1 == 1
        end

        if haskey(kwargs,:source)
            S = kwargs[:source]
            kwargs = Dict([p for p in pairs(kwargs) if p[1] != :source])
        else
            S = zeroFunc
        end
        if jac
            @assert size(F(U0[grid.indices[1],:])) == size(U0[grid.indices[1],:])
            @assert size(G(U0[grid.indices[1],:])) == size(U0[grid.indices[1],:])
            @assert size(H(U0[grid.indices[1],:])) == size(U0[grid.indices[1],:])
        else
            @assert size(F(U0[grid.indices[1],:]),1) == length(U0[grid.indices[1],:])
            @assert size(F(U0[grid.indices[1],:]),2) == length(U0[grid.indices[1],:])
            # TODO: add @assert for G and H
        end
        @assert size(S(U0)) == size(U0)

        #grid.U = reshape(U0,grid.nx,grid.ny,grid.nz,nvars)
        #grid.U = copy(U0)
        if haskey(kwargs,:parameters)
            θ = kwargs[:parameters]
            kwargs = Dict([p for p in pairs(kwargs) if p[1] != :parameters])
        else
            θ = nothing
        end
        if haskey(kwargs,:ul)
            U_min = kwargs[:ul]
            kwargs = Dict([p for p in pairs(kwargs) if p[1] != :ul])
        else
            U_min = -Inf * ones(size(U0))
        end
        if haskey(kwargs,:uh)
            U_max = kwargs[:uh]
            kwargs = Dict([p for p in pairs(kwargs) if p[1] != :uh])
        else
            U_max = Inf * ones(size(U0))
        end

        @assert size(U_min) == size(U0)
        @assert size(U_max) == size(U0)



        #@assert checkDim(D) == size(U0,)

        return new(grid,F, G, H, S, nvars, θ, U0, U_min, U_max,jac,kwargs,tspan,0)
    end


end

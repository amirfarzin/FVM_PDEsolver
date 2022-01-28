
function solve(prob::FVMPDEProblem;kwargs...)
    prob.t = 0
    if haskey(kwargs,:scheme)
        scheme = kwargs[:scheme]
        kwargs = Dict([p for p in pairs(kwargs) if p[1] != :scheme])
    else
        scheme = uw1
    end

    if haskey(kwargs,:showstep)
        dt_show = kwargs[:showstep]
        kwargs = Dict([p for p in pairs(kwargs) if p[1] != :showstep])
    else
        dt_show = Inf
    end

    U0 = copy(prob.U0)

    ODE_Function(u,p,t) = FVMPDE_∂u∂t(prob,u,t,scheme,dt_show)
    ODE_Problem = ODEProblem(ODE_Function,U0,prob.tspan;prob.ode_problem_kwargs...);

    if haskey(kwargs,:algorithm)
        alg = kwargs[:algorithm]
        kwargs = Dict([p for p in pairs(kwargs) if p[1] != :algorithm])
        sol = DifferentialEquations.solve(ODE_Problem,alg;kwargs...)
    else
        sol = DifferentialEquations.solve(ODE_Problem;kwargs...)
    end


    return sol
end

function FVMPDE_∂u∂t(prob::FVMPDEProblem,u,t,scheme,dt_show)
    if t - prob.t > dt_show
        println(t)
        prob.t = t
    end
    θ = prob.param
    indices = prob.grid.indices
    Δx = prob.grid.dX
    Δx_c = prob.grid.dX_C
    ∂F∂x = FVMPDE_∂F∂x(scheme,prob.F,u,prob.grid,1,t,θ,prob.jac)

    if prob.grid.dimension == _2D
        Δy = prob.grid.dY[1]
        ∂G∂y = FVMPDE_∂F∂x(scheme,prob.G,u,prob.grid,2,t,θ,prob.jac)
        ∂H∂z = 0
    elseif prob.grid.dimension == _3D
        Δy = prob.grid.dY[1]
        Δz = prob.grid.dZ[1]
        ∂G∂y = FVMPDE_∂F∂x(scheme,prob.G,u,prob.grid,2,t,θ,prob.jac)
        ∂H∂z = FVMPDE_∂F∂x(scheme,prob.H,u,prob.grid,3,t,θ,prob.jac)
    else
        ∂G∂y = 0
        ∂H∂z = 0
    end

    S = prob.S(u)

    ∂u∂t = S .- (∂F∂x .+ ∂G∂y .+ ∂H∂z)

    return ∂u∂t
end

function FVMPDE_∂F∂x(scheme,func_F,u,grid::FVMPDEGrid,_dim,t,p,jac)
    Uₓᴸ, Uₓᴿ = GeneralFluxLimiter(scheme,u,grid,_dim)


    ∂f∂x = zeros(size(u))

    for i in grid.indices
        A⁺ , A⁻ = LRJacobian(func_F,u[i,:],t,p,jac)
        ∂f∂x[i,:] = (A⁻ * Uₓᴸ[i,:] .+ A⁺ * Uₓᴿ[i,:])
    end
    return ∂f∂x
end


function LRJacobian(func_F,u,t,p,jac)
    if jac
        A = ForwardDiff.jacobian(func_F, u)
    else
        A = func_F(u)
    end

    if size(A,1) == 1
        A = reshape(A,1,1)
        A⁺ = 0.5 .* (A .+ abs.(A))
        A⁻ = 0.5 .* (A .- abs.(A))
    else
        Λ = diagm(eigvals(A))
        @assert typeof(Λ) === Matrix{Float64}
        Λ⁺ = 0.5 .* (Λ .+ abs.(Λ))
        Λ⁻ = 0.5 .* (Λ .- abs.(Λ))
        P = eigvecs(A)
        P⁻¹ = inv(P)
        A⁺ = P * Λ⁺ * P⁻¹
        A⁻ = P * Λ⁻ * P⁻¹


    end

    return A⁺, A⁻
end

function GeneralFluxLimiter(func_Ψ,u,grid::FVMPDEGrid,dim)
    ϵ = 1.0e-8
    indices = grid.indices
    Δx = grid.dX
    Δx_C = grid.dX_C

    func_r(Uₓ₁,Uₓ₂) = (Uₓ₁ .+ ϵ) ./ (Uₓ₂ .+ ϵ)
    func_U(Δxᵢ,r,uₓ) = (Δxᵢ ./ 2) .* func_Ψ(r) .* (uₓ .+ ϵ)

    d = length(size(u)) - 1
    #nx,ny,nz,nvars = size(u)

    minIndex = fill(Tuple(ones(Int8,d)),size(indices))
    maxIndex = fill(Tuple(size(u)[1:end-1]),size(indices))

    adderArray = zeros(Int8,d)
    adderArray[dim] = 1
    adder = CartesianIndex(Tuple(adderArray))

    ViewUᵢ₋₂ = CartesianIndex.(max.(Tuple.(indices .- adder .- adder),minIndex))
    ViewUᵢ₋₁ = CartesianIndex.(max.(Tuple.(indices .- adder         ),minIndex))
    ViewUᵢ₊₁ = CartesianIndex.(min.(Tuple.(indices .+ adder         ),maxIndex))
    ViewUᵢ₊₂ = CartesianIndex.(min.(Tuple.(indices .+ adder .+ adder),maxIndex))

    #println(size( (u[ViewUᵢ₊₂,:] .- u[ViewUᵢ₊₁,:])))
    #println(size(Δx))
    #println(size(Δx_C))

    Δx_extended = [Δx[1];Δx[1];Δx;Δx[end];Δx[end]]

    Uₓi₊¾ = (u[ViewUᵢ₊₂,:] .- u[ViewUᵢ₊₁,:]) ./ Δx_extended[4:end-0]    # ¾ : 3/2
    Uₓi₊½ = (u[ViewUᵢ₊₁,:] .- u) ./ Δx_extended[3:end-1]
    Uₓi₋½ = (u .- u[ViewUᵢ₋₁,:]) ./ Δx_extended[2:end-2]
    Uₓi₋¾ = (u[ViewUᵢ₋₁,:] .- u[ViewUᵢ₋₂,:]) ./ Δx_extended[1:end-3]

    Δx_C_extended = [Δx_C[1];Δx_C;Δx_C[end]]
    Uᴸi₊½ = u[ViewUᵢ₊₁,:] .- func_U(Δx_C_extended[3:end-0],func_r(Uₓi₊½,Uₓi₊¾),Uₓi₊¾)
    Uᴸi₋½ = u             .- func_U(Δx_C_extended[2:end-1],func_r(Uₓi₋½,Uₓi₊½),Uₓi₊½)
    Uᴿi₊½ = u             .+ func_U(Δx_C_extended[2:end-1],func_r(Uₓi₊½,Uₓi₋½),Uₓi₋½)
    Uᴿi₋½ = u[ViewUᵢ₋₁,:] .+ func_U(Δx_C_extended[1:end-2],func_r(Uₓi₋½,Uₓi₋¾),Uₓi₋¾)

    Uₓᴸ = (Uᴸi₊½ .- Uᴸi₋½) ./ Δx_C
    Uₓᴿ = (Uᴿi₊½ .- Uᴿi₋½) ./ Δx_C
    #return Uᴿi₊½, Uᴿi₋½, Uᴸi₊½ , Uᴸi₋½
    #return Uₓi₊¾, Uₓi₊½, Uₓi₋½, Uₓi₋¾
    return Uₓᴸ, Uₓᴿ
end





################## LIMITERS ####################

function κ_Scheme(r,κ)
    α₁ = (1+κ)/2
    α₂ = (1-κ)/2
    return (α₁ .*  r .+ α₂)
end


function uw1(r)
    return zeros(size(r))
end

function uw2(r)
    return ones(size(r))
end

function uw3(r)
    return κ_Scheme(r,1/3)
end

function uw4(r)
    return κ_Scheme(r,1/2)
end

function scd(r)
    return r
end

function fr(r)
    return κ_Scheme(r,0)
end

function kn(r)
    return @. max(0,min(2*r,min(1/3 + 2 /3 * r,2)))
end

function sb(r)
    return @. max(0,min(2 * r,1),min(r,2))
end

function mm(r)
    return @. max(0,min(r,1))
end

function mu(r)
    return @. max(0,min(2 * r,(r+1)/2,2))
end

function ha(r)
    return @. (r + abs(r))/(r+1)
end

function va1(r)
    return @. r*(r + 1)/(r^2+1)
end

function va2(r)
    return @. r*2/(r^2+1)
end

function vl(r)
    return @. (r + abs(r))/(1 + abs(r))
end

function op(r)
    return @. 3*r*(r + 1)/(r^2+r+1)/2
end

function hc(r)
    return @. 1.5*(r + abs(r))/(r+2)
end

function hq(r)
    return @. 2*(r + abs(r))/(r+3)
end

function cm(r)
    return @. max(0,r*(3*r+1)/(r+1)^2)
end

function mc(r)
    return @. max(0,min(2*r,0.5*(r+1),2))
end

function sm(r)
    return @. max(0,min(2*r,(0.75*r+0.25),4))
end

function um(r)
    return @. max(0,min(2*r,(0.25*r+0.75),2))
end

abstract type Dimension end
abstract type _1D <: Dimension end
abstract type _2D <: Dimension end
abstract type _3D <: Dimension end

function checkDim(D::Type{<:Dimension})
    if D === _1D
        return 1
    elseif D === _2D
        return 2
    elseif D === _3D
        return 3
    end
    return 0
end



mutable struct FVMPDEGrid
    X
    Y
    Z
    nx
    ny
    nz
    dX
    dY
    dZ
    dX_C
    dY_C
    dZ_C
    U
    Uboundries
    dimension::Type{<:Dimension}
    indices
    function FVMPDEGrid(X,Y,Z,nx,ny,nz,dX,dY,dZ,dX_C,dY_C,dZ_C,U,Uboundries,dimension,indices)
        return new(X,Y,Z,nx,ny,nz,dX,dY,dZ,dX_C,dY_C,dZ_C,U,Uboundries,dimension,indices)
    end
end

function FVMPDEGrid(_X::AbstractVector, _Y::AbstractVector, _Z::AbstractVector)
    X = _X
    Y = _Y
    Z = _Z
    Xb = [_X[1];(_X[1:end-1].+_X[2:end])./2;_X[end]]
    Yb = [_Y[1];(_Y[1:end-1].+_Y[2:end])./2;_Y[end]]
    Zb = [_Z[1];(_Z[1:end-1].+_Z[2:end])./2;_Z[end]]

    nx = length(X)
    ny = length(Y)
    nz = length(Z)
    dX = _X[2:end] .- _X[1:end-1]
    dY = _Y[2:end] .- _Y[1:end-1]
    dZ = _Z[2:end] .- _Z[1:end-1]
    dX_C = Xb[2:end] .- Xb[1:end-1]
    dY_C = Yb[2:end] .- Yb[1:end-1]
    dZ_C = Zb[2:end] .- Zb[1:end-1]
    #Uboundries = zeros(length(_X),length(_Y))
    dimension = _3D
    indices = CartesianIndices((Base.OneTo(nx),Base.OneTo(ny),Base.OneTo(nz)))
    return FVMPDEGrid(X,Y,Z,nx,ny,nz,dX,dY,dZ,dX_C,dY_C,dZ_C,nothing,nothing,dimension,indices)
end

function FVMPDEGrid(_X::AbstractVector, _Y::AbstractVector)
    X = _X
    Y = _Y
    Xb = [_X[1];(_X[1:end-1].+_X[2:end])./2;_X[end]]
    Yb = [_Y[1];(_Y[1:end-1].+_Y[2:end])./2;_Y[end]]
    nx = length(X)
    ny = length(Y)
    dX = _X[2:end] .- _X[1:end-1]
    dY = _Y[2:end] .- _Y[1:end-1]
    dX_C = Xb[2:end] .- Xb[1:end-1]
    dY_C = Yb[2:end] .- Yb[1:end-1]
    #Uboundries = zeros(length(_X),length(_Y))
    dimension = _2D
    indices = CartesianIndices((Base.OneTo(nx),Base.OneTo(ny)))
    return FVMPDEGrid(X,Y,[],nx,ny,1,dX,dY,nothing,dX_C,dY_C,nothing,nothing,nothing,dimension,indices)
end

function FVMPDEGrid(_X::AbstractVector)
    X = _X
    Xb = [_X[1];(_X[1:end-1].+_X[2:end])./2;_X[end]]
    nx = length(X)
    dX = _X[2:end] .- _X[1:end-1]
    dX_C = Xb[2:end] .- Xb[1:end-1]
    #Uboundries = zeros(length(_X))
    dimension = _1D
    indices = CartesianIndices((Base.OneTo(nx)))
    return FVMPDEGrid(X,[],[],nx,1,1,dX,nothing,nothing,dX_C,nothing,nothing,nothing,nothing,dimension,indices)
end

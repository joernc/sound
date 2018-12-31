# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion
# of Sound Waves" (JPO, 2009)

# next steps:
# - organize channels better?
# - exchange all sides?
# - note: rotational slips assume no-slip BC for now
# - allow stress BC (need to modify rotational split)

@everywhere using Printf
@everywhere using HDF5
@everywhere using Dates
@everywhere using Random

# grid spacing
@everywhere const Δx = 80.
@everywhere const Δy = Δx
@everywhere const Δz = 5.

# inertial frequency
@everywhere const f = -5e-5

# sound speed
@everywhere const c = 1.

# initial stratification
@everywhere const N = 1e-3

# time step
@everywhere const Δt = Δx/c

# aspect ratio
@everywhere const μ = Δz/Δx

# directory for data storage
@everywhere const datadir = "./data/rev"

# exchange of tile edges in x-direction
@everywhere function exchangex!(a::Array{Float64,3},
                                chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  @inbounds begin
    # send edges
    put!(chan_send_w, a[2,:,:])
    put!(chan_send_e, a[end-1,:,:])
    # receive edges
    a[1,:,:] = take!(chan_receive_w)
    a[end,:,:] = take!(chan_receive_e)
  end
end

# exchange of tile edges in y-direction
@everywhere function exchangey!(a::Array{Float64,3},
                                chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  @inbounds begin
    # send edges
    put!(chan_send_s, a[:,2,:])
    put!(chan_send_n, a[:,end-1,:])
    # receive edges
    a[:,1,:] = take!(chan_receive_s)
    a[:,end,:] = take!(chan_receive_n)
  end
end

# exchange of tile edges in z-direction
@everywhere function exchangez!(a::Array{Float64,3},
                                chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  @inbounds begin
    # send edges
    put!(chan_send_b, a[:,:,2])
    put!(chan_send_t, a[:,:,end-1])
    # receive edges
    a[:,:,1] = take!(chan_receive_b)
    a[:,:,end] = take!(chan_receive_t)
  end
end

# sound wave split (x-direction)
@everywhere function Sx!(u::Array{Float64,3}, ϕ::Array{Float64,3}, maskx::Array{UInt8,3},
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(u)
  # fields at initial time
  up = copy(u)
  ϕp = copy(ϕ)
  # loop over grid points
  @inbounds for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      u[i,j,k] = 1/2*(up[i-1,j,k] + up[i+1,j,k]) + 1/2c*(ϕp[i-1,j,k] - ϕp[i+1,j,k])
      ϕ[i,j,k] = 1/2*(ϕp[i-1,j,k] + ϕp[i+1,j,k]) + c/2*(up[i-1,j,k] - up[i+1,j,k])
    elseif maskx[i,j,k] == 2 # western boundary
      u[i,j,k] = 1/2*(-up[i,j,k] + up[i+1,j,k]) + 1/2c*(ϕp[i,j,k] - ϕp[i+1,j,k])
      ϕ[i,j,k] = 1/2*(ϕp[i,j,k] + ϕp[i+1,j,k]) + c/2*(-up[i,j,k] - up[i+1,j,k])
    elseif maskx[i,j,k] == 3 # eastern boundary
      u[i,j,k] = 1/2*(up[i-1,j,k] - up[i,j,k]) + 1/2c*(ϕp[i-1,j,k] - ϕp[i,j,k])
      ϕ[i,j,k] = 1/2*(ϕp[i-1,j,k] + ϕp[i,j,k]) + c/2*(up[i-1,j,k] + up[i,j,k])
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(ϕ, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# sound wave split (y-direction)
@everywhere function Sy!(v::Array{Float64,3}, ϕ::Array{Float64,3}, masky::Array{UInt8,3},
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(v)
  # fields at initial time
  vp = copy(v)
  ϕp = copy(ϕ)
  # loop over grid points
  @inbounds for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      v[i,j,k] = 1/2*(vp[i,j-1,k] + vp[i,j+1,k]) + 1/2c*(ϕp[i,j-1,k] - ϕp[i,j+1,k])
      ϕ[i,j,k] = 1/2*(ϕp[i,j-1,k] + ϕp[i,j+1,k]) + c/2*(vp[i,j-1,k] - vp[i,j+1,k])
    elseif masky[i,j,k] == 2 # southern boundary
      v[i,j,k] = 1/2*(-vp[i,j,k] + vp[i,j+1,k]) + 1/2c*(ϕp[i,j,k] - ϕp[i,j+1,k])
      ϕ[i,j,k] = 1/2*(ϕp[i,j,k] + ϕp[i,j+1,k]) + c/2*(-vp[i,j,k] - vp[i,j+1,k])
    elseif masky[i,j,k] == 3 # northern boundary
      v[i,j,k] = 1/2*(vp[i,j-1,k] - vp[i,j,k]) + 1/2c*(ϕp[i,j-1,k] - ϕp[i,j,k])
      ϕ[i,j,k] = 1/2*(ϕp[i,j-1,k] + ϕp[i,j,k]) + c/2*(vp[i,j-1,k] + vp[i,j,k])
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(ϕ, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# sound wave split (z-direction)
@everywhere function Sz!(w::Array{Float64,3}, ϕ::Array{Float64,3}, maskz::Array{UInt8,3},
                         chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(w)
  # fields at initial time
  wp = copy(w)
  ϕp = copy(ϕ)
  # loop over grid points
  @inbounds for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      w[i,j,k] = 1/2*(wp[i,j,k-1] + wp[i,j,k+1]) + μ/2c*(ϕp[i,j,k-1] - ϕp[i,j,k+1])
      ϕ[i,j,k] = 1/2*(ϕp[i,j,k-1] + ϕp[i,j,k+1]) + c/2μ*(wp[i,j,k-1] - wp[i,j,k+1])
    elseif maskz[i,j,k] == 2 # bottom boundary
      w[i,j,k] = 1/2*(-wp[i,j,k] + wp[i,j,k+1]) + μ/2c*(ϕp[i,j,k] - ϕp[i,j,k+1])
      ϕ[i,j,k] = 1/2*(ϕp[i,j,k] + ϕp[i,j,k+1]) + c/2μ*(-wp[i,j,k] - wp[i,j,k+1])
    elseif maskz[i,j,k] == 3 # top boundary
      w[i,j,k] = 1/2*(wp[i,j,k-1] - wp[i,j,k]) + μ/2c*(ϕp[i,j,k-1] - ϕp[i,j,k])
      ϕ[i,j,k] = 1/2*(ϕp[i,j,k-1] + ϕp[i,j,k]) + c/2μ*(wp[i,j,k-1] + wp[i,j,k])
    end
  end
  # exchange with neighboring tiles
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(ϕ, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

@everywhere function Tx_rhs!(u::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3},
                             kα::Array{Float64,3}, maskx::Array{UInt8,3},
                             chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(u)
  # calculate b
  b = c^2*α./ϕ
  # loop over grid points
  @inbounds for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      kα[i,j,k] = Δt/2Δx*((b[i-1,j,k] + b[i,j,k])*(u[i-1,j,k] + u[i,j,k]) - (b[i,j,k] + b[i+1,j,k])*(u[i,j,k] + u[i+1,j,k]))
    elseif maskx[i,j,k] == 2 # western boundary
      kα[i,j,k] = -Δt/2Δx*(b[i,j,k] + b[i+1,j,k])*(u[i,j,k] + u[i+1,j,k])
    elseif maskx[i,j,k] == 3 # eastern boundary
      kα[i,j,k] = Δt/2Δx*(b[i-1,j,k] + b[i,j,k])*(u[i-1,j,k] + u[i,j,k])
    else
      kα[i,j,k] = 0.
    end
  end
  # exchange with neighboring tiles
  exchangex!(kα, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

@everywhere function Ty_rhs!(v::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3},
                             kα::Array{Float64,3}, masky::Array{UInt8,3},
                             chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(v)
  # calculate b
  b = c^2*α./ϕ
  # loop over grid points
  @inbounds for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      kα[i,j,k] = Δt/2Δy*((b[i,j-1,k] + b[i,j,k])*(v[i,j-1,k] + v[i,j,k]) - (b[i,j,k] + b[i,j+1,k])*(v[i,j,k] + v[i,j+1,k]))
    elseif masky[i,j,k] == 2 # southern boundary
      kα[i,j,k] = -Δt/2Δy*(b[i,j,k] + b[i,j+1,k])*(v[i,j,k] + v[i,j+1,k])
    elseif masky[i,j,k] == 3 # northern boundary
      kα[i,j,k] = Δt/2Δy*(b[i,j-1,k] + b[i,j,k])*(v[i,j-1,k] + v[i,j,k])
    else
      kα[i,j,k] = 0.
    end
  end
  # exchange with neighboring tiles
  exchangey!(kα, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

@everywhere function Tz_rhs!(w::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3},
                             kw::Array{Float64,3}, kα::Array{Float64,3}, maskz::Array{UInt8,3},
                             chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                             chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(w)
  # calculate b
  b = c^2*α./ϕ
  # loop over grid points
  @inbounds for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      kw[i,j,k] = Δt*μ^2/2*(b[i,j,k-1] + 2b[i,j,k] + b[i,j,k+1])
      kα[i,j,k] = Δt/2Δz*((b[i,j,k-1] + b[i,j,k])*(w[i,j,k-1] + w[i,j,k]) - (b[i,j,k] + b[i,j,k+1])*(w[i,j,k] + w[i,j,k+1]))
    elseif maskz[i,j,k] == 2 # bottom boundary
      kw[i,j,k] = Δt*μ^2/2*(b[i,j,k] + b[i,j,k+1])
      kα[i,j,k] = -Δt/2Δz*(b[i,j,k] + b[i,j,k+1])*(w[i,j,k] + w[i,j,k+1])
    elseif maskz[i,j,k] == 3 # top boundary
      kw[i,j,k] = Δt*μ^2/2*(b[i,j,k-1] + b[i,j,k])
      kα[i,j,k] = Δt/2Δz*(b[i,j,k-1] + b[i,j,k])*(w[i,j,k-1] + w[i,j,k])
    else
      kw[i,j,k] = 0.
      kα[i,j,k] = 0.
    end
  end
  # exchange with neighboring tiles
  exchangez!(kw, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(kα, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# buoyancy split (x-direction)
@everywhere function Tx!(u::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3}, maskx::Array{UInt8,3},
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(u)
  # allocate
  kα1 = Array{Float64,3}(undef, nx, ny, nz)
  kα2 = Array{Float64,3}(undef, nx, ny, nz)
  # perform midpoint
  Tx_rhs!(u, ϕ, α, kα1, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  Tx_rhs!(u, ϕ, α + kα1/2, kα2, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  α .+= kα2
end

# buoyancy split (y-direction)
@everywhere function Ty!(v::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3}, masky::Array{UInt8,3},
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(v)
  # allocate
  kα1 = Array{Float64,3}(undef, nx, ny, nz)
  kα2 = Array{Float64,3}(undef, nx, ny, nz)
  # perform midpoint
  Ty_rhs!(v, ϕ, α, kα1, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  Ty_rhs!(v, ϕ, α + kα1/2, kα2, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  α .+= kα2
end

# buoyancy split (z-direction, advection and gravity)
@everywhere function Tz!(w::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3}, maskz::Array{UInt8,3},
                          chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                          chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                          chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                          chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(w)
  # allocate
  kα1 = Array{Float64,3}(undef, nx, ny, nz)
  kα2 = Array{Float64,3}(undef, nx, ny, nz)
  kw1 = Array{Float64,3}(undef, nx, ny, nz)
  kw2 = Array{Float64,3}(undef, nx, ny, nz)
  # perform midpoint
  Tz_rhs!(w, ϕ, α, kw1, kα1, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  Tz_rhs!(w + kw1/2, ϕ, α + kα1/2, kw2, kα2, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  w .+= kw2
  α .+= kα2
end

# rotation split (in y–z plane)
@everywhere function Rx!(v::Array{Float64,3}, w::Array{Float64,3}, ϕ::Array{Float64,3},
                         masky::Array{UInt8,3}, maskz::Array{UInt8,3},
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(v)
  # fields at initial time
  vp = copy(v)
  wp = copy(w)
  # loop over grid points
  @inbounds for k = 2:nz-1, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] > 0 # fluid
      if masky[i,j,k] == 1 # interior
        wy = 1/2Δy*(wp[i,j+1,k] - wp[i,j-1,k])
      elseif masky[i,j,k] == 2 # westerm boundary
        wy = 1/2Δy*(wp[i,j+1,k] + wp[i,j,k])
      elseif masky[i,j,k] == 3 # eastern boundary
        wy = -1/2Δy*(wp[i,j,k] + wp[i,j-1,k])
      else # boundaries on both sides
        wy = 0.
      end
      if maskz[i,j,k] == 1 # interior
        vz = 1/2Δz*(vp[i,j,k+1] - vp[i,j,k-1])
      elseif maskz[i,j,k] == 2 # bottom boundary
        vz = 1/2Δz*(vp[i,j,k+1] + vp[i,j,k])
      elseif maskz[i,j,k] == 3 # top boundary
        vz = -1/2Δz*(vp[i,j,k] + vp[i,j,k-1])
      else # boundaries on both sides
        vz = 0.
      end
      γx = 2Δt*c^2/μ*(wy - μ^2*vz)/ϕ[i,j,k]
      Cx = (1 - 1/4*γx.^2)/(1 + 1/4*γx^2)
      Sx = γx/(1 + 1/4*γx^2)
      v[i,j,k] = vp[i,j,k]*Cx + 1/μ*wp[i,j,k]*Sx
      w[i,j,k] = wp[i,j,k]*Cx - μ*vp[i,j,k]*Sx
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(w, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangez!(v, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# rotation split (in x–z plane)
@everywhere function Ry!(u::Array{Float64,3}, w::Array{Float64,3}, ϕ::Array{Float64,3},
                         maskx::Array{UInt8,3}, maskz::Array{UInt8,3},
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(u)
  # fields at initial time
  up = copy(u)
  wp = copy(w)
  # loop over grid points
  @inbounds for k = 2:nz-1, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] > 0 # fluid
      if maskz[i,j,k] == 1 # interior
        uz = 1/2Δz*(up[i,j,k+1] - up[i,j,k-1])
      elseif maskz[i,j,k] == 2 # bottom boundary
        uz = 1/2Δz*(up[i,j,k+1] + up[i,j,k])
      elseif maskz[i,j,k] == 3 # top boundary
        uz = -1/2Δz*(up[i,j,k] + up[i,j,k-1])
      else # boundaries on both sides
        uz = 0.
      end
      if maskx[i,j,k] == 1 # interior
        wx = 1/2Δx*(wp[i+1,j,k] - wp[i-1,j,k])
      elseif maskx[i,j,k] == 2 # westerm boundary
        wx = 1/2Δx*(wp[i+1,j,k] + wp[i,j,k])
      elseif maskx[i,j,k] == 3 # eastern boundary
        wx = -1/2Δx*(wp[i,j,k] + wp[i-1,j,k])
      else # boundaries on both sides
        wx = 0.
      end
      γy = 2Δt*c^2/μ*(μ^2*uz - wx)/ϕ[i,j,k]
      Cy = (1 - 1/4*γy.^2)/(1 + 1/4*γy^2)
      Sy = γy/(1 + 1/4*γy^2)
      u[i,j,k] = up[i,j,k]*Cy - 1/μ*wp[i,j,k]*Sy
      w[i,j,k] = wp[i,j,k]*Cy + μ*up[i,j,k]*Sy
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(w, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(u, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# rotation split (in x–y plane)
@everywhere function Rz!(u::Array{Float64,3}, v::Array{Float64,3}, ϕ::Array{Float64,3},
                         maskx::Array{UInt8,3}, masky::Array{UInt8,3},
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(u)
  # fields at initial time
  up = copy(u)
  vp = copy(v)
  # loop over grid points
  @inbounds for k = 1:nz, j = 2:ny-1, i = 2:nx-1
    if maskx[i,j,k] > 0 # fluid
      if maskx[i,j,k] == 1 # interior
        vx = 1/2Δx*(vp[i+1,j,k] - vp[i-1,j,k])
      elseif maskx[i,j,k] == 2 # westerm boundary
        vx = 1/2Δx*(vp[i+1,j,k] + vp[i,j,k])
      elseif maskx[i,j,k] == 3 # eastern boundary
        vx = -1/2Δx*(vp[i,j,k] + vp[i-1,j,k])
      else # boundaries on both sides
        vx = 0.
      end
      if masky[i,j,k] == 1 # interior
        uy = 1/2Δy*(up[i,j+1,k] - up[i,j-1,k])
      elseif masky[i,j,k] == 2 # southern boundary
        uy = 1/2Δy*(up[i,j+1,k] + up[i,j,k])
      elseif masky[i,j,k] == 3 # northern boundary
        uy = -1/2Δy*(up[i,j,k] + up[i,j-1,k])
      else # boundaries on both sides
        uy = 0.
      end
      γz = 2Δt*c^2*(f + vx - uy)/ϕ[i,j,k]
      Cz = (1 - 1/4*γz.^2)/(1 + 1/4*γz^2)
      Sz = γz/(1 + 1/4*γz^2)
      u[i,j,k] = up[i,j,k]*Cz + vp[i,j,k]*Sz
      v[i,j,k] = vp[i,j,k]*Cz - up[i,j,k]*Sz
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(v, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(u, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# sound and buoyancy splits (x-direction)
@everywhere function STSx!(u::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3}, maskx::Array{UInt8,3},
                           chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  Sx!(u, ϕ, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  Tx!(u, ϕ, α, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  Sx!(u, ϕ, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# sound and buoyancy splits (y-direction)
@everywhere function STSy!(v::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3}, masky::Array{UInt8,3},
                           chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  Sy!(v, ϕ, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  Ty!(v, ϕ, α, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  Sy!(v, ϕ, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# sound and buoyancy splits (z-direction)
@everywhere function STSz!(w::Array{Float64,3}, ϕ::Array{Float64,3}, α::Array{Float64,3}, maskz::Array{UInt8,3},
                           chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                           chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  Sz!(w, ϕ, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  Tz!(w, ϕ, α, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  Sz!(w, ϕ, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# x-diffusion
@everywhere function Dx!(a::Array{Float64,3}, κ::Array{Float64,3}, maskx::Array{UInt8,3}, diri::Bool,
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  @inbounds for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      a[i,j,k] += Δt/Δz^2*((κ[i-1,j,k] + κ[i,j,k])*(ap[i-1,j,k] - ap[i,j,k]) - (κ[i,j,k] + κ[i+1,j,k])*(ap[i,j,k] - ap[i+1,j,k]))
    elseif maskx[i,j,k] == 2 # western boundary
      if diri # Dirichlet BC
        a[i,j,k] -= Δt/Δz^2*((κ[i-1,j,k] + κ[i,j,k])*2ap[i,j,k] + (κ[i,j,k] + κ[i+1,j,k])*(ap[i,j,k] - ap[i+1,j,k]))
      else # Neumann BC
        a[i,j,k] -= Δt/Δz^2*(κ[i,j,k] + κ[i+1,j,k])*(ap[i,j,k] - ap[i+1,j,k])
      end
    elseif maskx[i,j,k] == 3 # eastern boundary
      if diri
        a[i,j,k] += Δt/Δz^2*((κ[i-1,j,k] + κ[i,j,k])*(ap[i-1,j,k] - ap[i,j,k]) - (κ[i,j,k] + κ[i+1,j,k])*2ap[i,j,k])
      else
        a[i,j,k] += Δt/Δz^2*(κ[i-1,j,k] + κ[i,j,k])*(ap[i-1,j,k] - ap[i,j,k])
      end
    elseif maskx[i,j,k] == 4 # boundaries on both sides
      if diri
        a[i,j,k] -= 2Δt/Δz^2*(κ[i-1,j,k] + 2κ[i,j,k] + κ[i+1,j,k])*ap[i,j,k]
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# y-diffusion
@everywhere function Dy!(a::Array{Float64,3}, κ::Array{Float64,3}, masky::Array{UInt8,3}, diri::Bool,
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  @inbounds for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      a[i,j,k] += Δt/Δz^2*((κ[i,j-1,k] + κ[i,j,k])*(ap[i,j-1,k] - ap[i,j,k]) - (κ[i,j,k] + κ[i,j+1,k])*(ap[i,j,k] - ap[i,j+1,k]))
    elseif masky[i,j,k] == 2 # southern boundary
      if diri
        a[i,j,k] -= Δt/Δz^2*((κ[i,j-1,k] + κ[i,j,k])*2ap[i,j,k] + (κ[i,j,k] + κ[i,j+1,k])*(ap[i,j,k] - ap[i,j+1,k]))
      else
        a[i,j,k] -= Δt/Δz^2*(κ[i,j,k] + κ[i,j+1,k])*(ap[i,j,k] - ap[i,j+1,k])
      end
    elseif masky[i,j,k] == 3 # northern boundary
      if diri
        a[i,j,k] += Δt/Δz^2*((κ[i,j-1,k] + κ[i,j,k])*(ap[i,j-1,k] - ap[i,j,k]) - (κ[i,j,k] + κ[i,j+1,k])*2ap[i,j,k])
      else
        a[i,j,k] += Δt/Δz^2*(κ[i,j-1,k] + κ[i,j,k])*(ap[i,j-1,k] - ap[i,j,k])
      end
    elseif masky[i,j,k] == 4 # boundaries on both sides
      if diri
        a[i,j,k] -= 2Δt/Δz^2*(κ[i,j-1,k] + 2κ[i,j,k] + κ[i,j+1,k])*ap[i,j,k]
      end
    end
  end
  # exchange with neighboring tiles
  exchangey!(a, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# z-diffusion
@everywhere function Dz!(a::Array{Float64,3}, κ::Array{Float64,3}, maskz::Array{UInt8,3}, diri::Bool,
                         chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  @inbounds for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      a[i,j,k] += Δt/Δz^2*((κ[i,j,k-1] + κ[i,j,k])*(ap[i,j,k-1] - ap[i,j,k]) - (κ[i,j,k] + κ[i,j,k+1])*(ap[i,j,k] - ap[i,j,k+1]))
    elseif maskz[i,j,k] == 2 # bottom boundary
      if diri
        a[i,j,k] -= Δt/Δz^2*((κ[i,j,k-1] + κ[i,j,k])*2ap[i,j,k] + (κ[i,j,k] + κ[i,j,k+1])*(ap[i,j,k] - ap[i,j,k+1]))
      else
        a[i,j,k] -= Δt/Δz^2*(κ[i,j,k] + κ[i,j,k+1])*(ap[i,j,k] - ap[i,j,k+1])
      end
    elseif maskz[i,j,k] == 3 # top boundary
      if diri
        a[i,j,k] += Δt/Δz^2*((κ[i,j,k-1] + κ[i,j,k])*(ap[i,j,k-1] - ap[i,j,k]) - (κ[i,j,k] + κ[i,j,k+1])*2ap[i,j,k])
      else
        a[i,j,k] += Δt/Δz^2*(κ[i,j,k-1] + κ[i,j,k])*(ap[i,j,k-1] - ap[i,j,k])
      end
    elseif maskz[i,j,k] == 4 # boundaries on both sides
      if diri
        a[i,j,k] -= 2Δt/Δz^2*(κ[i,j,k-1] + 2κ[i,j,k] + κ[i,j,k+1])*ap[i,j,k]
      end
    end
  end
  # exchange with neighboring tiles
  exchangez!(a, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# exchange all edges of global array
@everywhere function exchange_global!(a)
  a[1,:,:] = a[end-1,:,:]
  a[end,:,:] = a[2,:,:]
  a[:,1,:] = a[:,end-1,:]
  a[:,end,:] = a[:,2,:]
  a[:,:,1] = a[:,:,end-1]
  a[:,:,end] = a[:,:,2]
end

# read topography
@everywhere function read_topo()
  # global domain size
  nx = 256; ny = 1; nz = 256
  # read abyssal hill topography from file
  z = [(k-.5)*Δz for i = 1:nx, j = 1:ny, k = 1:nz]
  h = h5read("abyssal_256x1_80x80.h5", "h")
  # find fluid points
  fluid = BitArray(undef, nx+2, ny+2, nz+2)
  fluid[2:end-1,2:end-1,2:end-1] = z .> h
  fluid[:,:,end-1] .= false
  # fill edges
  exchange_global!(fluid)
  return fluid
end

# identify interior and boundary points
@everywhere function masks(fluid::BitArray{3}, irange, jrange, krange)
  # domain size
  nx, ny, nz = size(fluid)
  # find and categorize boundary points
  maskx = Array{UInt8,3}(undef, nx, ny, nz)
  masky = Array{UInt8,3}(undef, nx, ny, nz)
  maskz = Array{UInt8,3}(undef, nx, ny, nz)
  for k = 2:nz-1, j = 2:ny-1, i = 2:nx-1
    if fluid[i,j,k] # in fluid
      # x-direction
      if fluid[i-1,j,k] & fluid[i+1,j,k]
        maskx[i,j,k] = 1 # interior
      elseif fluid[i+1,j,k]
        maskx[i,j,k] = 2 # western boundary
      elseif fluid[i-1,j,k]
        maskx[i,j,k] = 3 # eastern boundary
      else
        maskx[i,j,k] = 4 # boundaries on both sides
      end
      # y-direction
      if fluid[i,j-1,k] & fluid[i,j+1,k]
        masky[i,j,k] = 1 # interior
      elseif fluid[i,j+1,k]
        masky[i,j,k] = 2 # southern boundary
      elseif fluid[i,j-1,k]
        masky[i,j,k] = 3 # northern boundary
      else
        masky[i,j,k] = 4 # boundaries on both sides
      end
      # z-direction
      if fluid[i,j,k-1] & fluid[i,j,k+1]
        maskz[i,j,k] = 1 # interior
      elseif fluid[i,j,k+1]
        maskz[i,j,k] = 2 # bottom boundary
      elseif fluid[i,j,k-1]
        maskz[i,j,k] = 3 # top boundary
      else
        maskz[i,j,k] = 4 # boundaries on both sides
      end
    else # not in fluid
      maskx[i,j,k] = 0
      masky[i,j,k] = 0
      maskz[i,j,k] = 0
    end
  end
  # exchange edges
  exchange_global!(maskx)
  exchange_global!(masky)
  exchange_global!(maskz)
  # select range of tile
  maskx = maskx[irange[1]:irange[end]+2,jrange[1]:jrange[end]+2,krange[1]:krange[end]+2]
  masky = masky[irange[1]:irange[end]+2,jrange[1]:jrange[end]+2,krange[1]:krange[end]+2]
  maskz = maskz[irange[1]:irange[end]+2,jrange[1]:jrange[end]+2,krange[1]:krange[end]+2]
  return maskx, masky, maskz
end

# set diffusivity/viscosity
@everywhere function read_visc(irange, jrange, krange,
                               chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                               chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # diffusion/viscosity parameters
  ν0 = 1e-2
  ν1 = 0.
#  ν0 = 5.0e-5
#  ν1 = 2.5e-3
  d = 200.
#  H = (512-1)*Δz
#  D = 50.
  # read topography
  h = h5read("abyssal_256x1_80x80.h5", "h")[irange,jrange]
  # tile size
  nx = length(irange); ny = length(jrange); nz = length(krange)
  # initialize viscosity field and its gradients
  ν = Array{Float64}(undef, nx+2, ny+2, nz+2)
  # calculate viscosity
  z = [(k-.5)*Δz for i = irange, j = jrange, k = krange]
  ν[2:end-1,2:end-1,2:end-1] .= ν0 .+ ν1*exp.(-(z.-h)/d) #.+ ν1*exp.((z.-H)/D)
  # exchange with neighboring tiles
  exchangex!(ν, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(ν, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangez!(ν, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  return ν
end

# save tile to file
@everywhere function save_tile(n, i, j, k, u::Array{Float64,3}, v::Array{Float64,3}, w::Array{Float64,3}, ϕ::Array{Float64,3},
                               b::Array{Float64,3}, maskx::Array{UInt8,3}, masky::Array{UInt8,3}, maskz::Array{UInt8,3})
  # tile size
  nx, ny, nz = size(u)
  # discard edges
  us = u[2:nx-1,2:ny-1,2:nz-1]
  vs = v[2:nx-1,2:ny-1,2:nz-1]
  ws = w[2:nx-1,2:ny-1,2:nz-1]
  ϕs = ϕ[2:nx-1,2:ny-1,2:nz-1]
  bs = b[2:nx-1,2:ny-1,2:nz-1]
  # replace missing values with NaNs
  solid = ((maskx .== 0) .& (masky .== 0) .& (maskz .== 0))[2:nx-1,2:ny-1,2:nz-1] # !!!
  us[solid] .= NaN
  vs[solid] .= NaN
  ws[solid] .= NaN
  ϕs[solid] .= NaN
  bs[solid] .= NaN
  # save data
  filename = @sprintf("%s/%010d_%1d_%1d_%1d.h5", datadir, n, i, j, k)
  h5write(filename, "u", us)
  h5write(filename, "v", vs)
  h5write(filename, "w", ws)
  h5write(filename, "ϕ", ϕs)
  h5write(filename, "b", bs)
  # screen print
  us[solid] .= 0.
  vs[solid] .= 0.
  ws[solid] .= 0.
  ϕs[solid] .= 0.
  bs[solid] .= 0.
  println(@sprintf("%s %6i %9.3e %9.3e %9.3e %9.3e", floor(now(), Second), n, 4n*Δt,
                   maximum(abs.(us)), maximum(abs.(vs)), maximum(abs.(ws))))
end

# load tile from file
@everywhere function load_tile(n, i, j, k, maskx::Array{UInt8,3}, masky::Array{UInt8,3}, maskz::Array{UInt8,3},
                              chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(maskx)
  # initialize arrays
  u = Array{Float64,3}(undef, nx, ny, nz)
  v = Array{Float64,3}(undef, nx, ny, nz)
  w = Array{Float64,3}(undef, nx, ny, nz)
  ϕ = Array{Float64,3}(undef, nx, ny, nz)
  b = Array{Float64,3}(undef, nx, ny, nz)
  # load data
  filename = @sprintf("%s/%010d_%1d_%1d_%1d.h5", datadir, n, i, j, k)
  u[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "u")
  v[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "v")
  w[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "w")
  ϕ[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "ϕ")
  b[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "b")
  # replace missing values with zeros
  solid = (maskx .== 0) .& (masky .== 0) .& (maskz .== 0) # !!!
  u[solid] .= 0.
  v[solid] .= 0.
  w[solid] .= 0.
  ϕ[solid] .= 0.
  b[solid] .= 0.
  # exchange edges with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(v, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(w, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(ϕ, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(u, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(w, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(ϕ, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(b, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangez!(u, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(v, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(ϕ, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(b, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # screen print (should really be truncated fields)
  println(@sprintf("%s %6i %9.3e %9.3e %9.3e %9.3e", floor(now(), Second), n, 4n*Δt,
                   maximum(abs.(u)), maximum(abs.(v)), maximum(abs.(w))))
  return u, v, w, ϕ, b
end

# run model on tile
@everywhere function run_tile(i, j, k, irange, jrange, krange, steps,
                              chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                              chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx = length(irange) + 2
  ny = length(jrange) + 2
  nz = length(krange) + 2
  # read topography and generate masks
  fluid = read_topo()
  maskx, masky, maskz = masks(fluid, irange, jrange, krange)
  maski = maskx .> 0
  # read viscosity/diffusivity maps
  ν = read_visc(irange, jrange, krange, chan_send_w, chan_send_e, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
                chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
#  # specify where to impose Dirichlet boundary conditions
#  diriu = .!maski
#  diriv = .!maski
#  diriw = .!maski
#  dirib = falses(nx, ny, nz)
#  diriϕ = falses(nx, ny, nz)
  # coordinates
  x = [(i-.5)*Δy for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1, k = krange[1]-1:krange[end]+1]
  y = [(j-.5)*Δy for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1, k = krange[1]-1:krange[end]+1]
  z = [(k-.5)*Δz for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1, k = krange[1]-1:krange[end]+1]
  if steps[1] == 1 # initialize
    # save masks
    h5write(@sprintf("%s/masks_%1d_%1d_%1d.h5", datadir, i, j, k), "maskx", maskx[2:nx-1,2:ny-1,2:nz-1])
    h5write(@sprintf("%s/masks_%1d_%1d_%1d.h5", datadir, i, j, k), "masky", masky[2:nx-1,2:ny-1,2:nz-1])
    h5write(@sprintf("%s/masks_%1d_%1d_%1d.h5", datadir, i, j, k), "maskz", maskz[2:nx-1,2:ny-1,2:nz-1])
    # initial conditions
    u = zeros(nx, ny, nz)
    v = zeros(nx, ny, nz)
    w = zeros(nx, ny, nz)
    b = N^2*z
    ϕ = c^2*ones(nx, ny, nz) + 1/2*N^2*z.^2
#    exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
#    exchangey!(u, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
#    exchangez!(u, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
#    exchangex!(v, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
#    exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
#    exchangez!(v, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
#    exchangex!(w, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
#    exchangey!(w, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
#    exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    # save initial conditions to file
    save_tile(0, i, j, k, u, v, w, ϕ, b, maskx, masky, maskz)
  else # load
    u, v, w, ϕ, b = load_tile(steps[1] - 1, i, j, k, maskx, masky, maskz, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
                              chan_send_b, chan_send_t, chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n,
                              chan_receive_b, chan_receive_t)
  end
  # time steps
  for n = steps
    
    # Strang splitting
  
    # u-diffusion
    Dx!(u, ν, maskx, true, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dy!(u, ν, masky, true, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dz!(u, ν, maskz, true, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)

    # v-diffusion
    Dx!(v, ν, maskx, true, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dy!(v, ν, masky, true, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dz!(v, ν, maskz, true, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    
    # w-diffusion
    Dx!(w, ν, maskx, true, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dy!(w, ν, masky, true, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dz!(w, ν, maskz, true, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    
    # b-diffusion
    Dx!(b, ν, maskx, false, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dy!(b, ν, masky, false, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dz!(b, ν, maskz, false, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    
    # rotation
    Rx!(v, w, ϕ, masky, maskz,
        chan_send_s, chan_send_n, chan_send_b, chan_send_t, chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
    Ry!(u, w, ϕ, maskx, maskz,
        chan_send_w, chan_send_e, chan_send_b, chan_send_t, chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    Rz!(u, v, ϕ, maskx, masky,
        chan_send_w, chan_send_e, chan_send_s, chan_send_n, chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)

    # sound and buoyancy
    α = ϕ.*b/c^2
    STSx!(u, ϕ, α, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    STSy!(v, ϕ, α, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    STSz!(w, ϕ, α, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    STSz!(w, ϕ, α, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    STSy!(v, ϕ, α, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    STSx!(u, ϕ, α, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    b = c^2*α./ϕ

    # rotation
    Rz!(u, v, ϕ, maskx, masky,
        chan_send_w, chan_send_e, chan_send_s, chan_send_n, chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
    Ry!(u, w, ϕ, maskx, maskz,
        chan_send_w, chan_send_e, chan_send_b, chan_send_t, chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    Rx!(v, w, ϕ, masky, maskz,
        chan_send_s, chan_send_n, chan_send_b, chan_send_t, chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)

    # b-diffusion
    Dz!(b, ν, maskz, false, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dy!(b, ν, masky, false, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dx!(b, ν, maskx, false, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    
    # w-diffusion
    Dz!(w, ν, maskz, true, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dy!(w, ν, masky, true, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dx!(w, ν, maskx, true, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    
    # v-diffusion
    Dz!(v, ν, maskz, true, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dy!(v, ν, masky, true, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dx!(v, ν, maskx, true, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)

    # u-diffusion
    Dz!(u, ν, maskz, true, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dy!(u, ν, masky, true, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dx!(u, ν, maskx, true, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)

    if n % 100 == 0
      # save tile to file
      save_tile(n, i, j, k, u, v, w, ϕ, b, maskx, masky, maskz)
    else
      println(@sprintf("%s %6i %9.3e", floor(now(), Second), n, 4n*Δt))
    end
  end
end

# get index ranges for tiles
tile_range(i, j, k, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
                                   sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j]),
                                   sum(tile_sizes[3][1:k-1])+1:sum(tile_sizes[3][1:k])]

# run the model
function run_model(steps, tile_sizes)
  # number of tiles
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  nk = length(tile_sizes[3])
  # assign tiles to workers
  proc = Array{Int}(undef, ni, nj, nk)
  for k = 1:nk, j = 1:nj, i = 1:ni
    proc_idx = mod1((k-1)*ni*nj + (j-1)*ni + i, nprocs() - 1)
    proc[i,j,k] = workers()[proc_idx]
  end
  # set up remote channels for edge communication
  chan_w = [RemoteChannel(()->Channel{Array{Float64,2}}(1), proc[i,j,k]) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_e = [RemoteChannel(()->Channel{Array{Float64,2}}(1), proc[i,j,k]) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_s = [RemoteChannel(()->Channel{Array{Float64,2}}(1), proc[i,j,k]) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_n = [RemoteChannel(()->Channel{Array{Float64,2}}(1), proc[i,j,k]) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_b = [RemoteChannel(()->Channel{Array{Float64,2}}(1), proc[i,j,k]) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_t = [RemoteChannel(()->Channel{Array{Float64,2}}(1), proc[i,j,k]) for i = 1:ni, j = 1:nj, k = 1:nk]
  # run model on tiles
  @sync begin
    for k = 1:nk, j = 1:nj, i = 1:ni
      # index range in global domain
      irange, jrange, krange = tile_range(i, j, k, tile_sizes)
      println(@sprintf("worker %1d: tile %1d,%1d,%1d", proc[i,j,k], i, j, k))
      println(irange); println(jrange); println(krange)
      # indices of neighboring tiles
      iw = mod1(i-1, ni); ie = mod1(i+1, ni)
      js = mod1(j-1, nj); jn = mod1(j+1, nj)
      kb = mod1(k-1, nk); kt = mod1(k+1, nk)
      # run model on tile
      @async remotecall_fetch(run_tile, proc[i,j,k], i, j, k, irange, jrange, krange, steps, chan_w[i,j,k], chan_e[i,j,k],
                              chan_s[i,j,k], chan_n[i,j,k], chan_b[i,j,k], chan_t[i,j,k], chan_e[iw,j,k], chan_w[ie,j,k],
                              chan_n[i,js,k], chan_s[i,jn,k], chan_t[i,j,kb], chan_b[i,j,kt])
    end
  end
end

## mass
#@everywhere function mass(ϕ, maski, maskx, masky, maskz)
#  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
#  return sum(ϕ[maski])/c^2 + sum(ϕ[maskb])/2c^2
#end
#
## energy
#@everywhere function energy(u, v, w, ϕ, b, z, maski, maskx, masky, maskz)
#  α = ϕ.*b/c^2
#  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
#  return sum(u[maski].^2)/2 + sum(v[maski].^2)/2 + sum(w[maski].^2)/2μ^2 + sum(ϕ[maski].^2)/2c^2 + sum(ϕ[maskb].^2)/4c^2 - sum(α[maski].*z[maski]) - sum(α[maskb].*z[maskb])/2
#end
#
## buoyancy
#@everywhere function buoyancy(ϕ, b, maski, maskx, masky, maskz)
#  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
#  return sum(ϕ[maski].*b[maski]) + sum(ϕ[maskb].*b[maskb])/2
#end
#
## buoyancy variance
#@everywhere function buoyancy_variance(ϕ, b, maski, maskx, masky, maskz)
#  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
#  return sum(ϕ[maski].*b[maski].^2) + sum(ϕ[maskb].*b[maskb].^2)/2
#end

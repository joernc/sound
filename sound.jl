# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion
# of Sound Waves" (JPO, 2009)

# next steps:
# - allow stress BC (need to modify rotational split)
# - organize channels better?

@everywhere using Printf
@everywhere using HDF5
@everywhere using Dates

# grid spacing
@everywhere const Δx = 64e3/256
@everywhere const Δy = Δx
@everywhere const Δz = 1024/256

# vertical viscosity/diffusion
@everywhere const ν = 1e-3

# inertial frequency
@everywhere const f = 5e-5

# sound speed
@everywhere const c = 1.

# time step
@everywhere const Δt = Δx/c

# aspect ratio
@everywhere const μ = Δz/Δx

# constant for diffusion steps
@everywhere const α = ν*2Δt/Δz^2

@everywhere function exchangex!(a::Array{Float64,3},
                                chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  # send edges
  put!(chan_send_w, a[2,:,:])
  put!(chan_send_e, a[end-1,:,:])
  # receive edges
  a[1,:,:] = take!(chan_receive_w)
  a[end,:,:] = take!(chan_receive_e)
end

@everywhere function exchangey!(a::Array{Float64,3},
                                chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # send edges
  put!(chan_send_s, a[:,2,:])
  put!(chan_send_n, a[:,end-1,:])
  # receive edges
  a[:,1,:] = take!(chan_receive_s)
  a[:,end,:] = take!(chan_receive_n)
end

@everywhere function exchangez!(a::Array{Float64,3},
                                chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                                chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # send edges
  put!(chan_send_b, a[:,:,2])
  put!(chan_send_t, a[:,:,end-1])
  # receive edges
  a[:,:,1] = take!(chan_receive_b)
  a[:,:,end] = take!(chan_receive_t)
end

# sound wave split (x-direction)
@everywhere function Sx!(u::Array{Float64,3}, ϕ::Array{Float64,3}, b::Array{Float64,3}, maskx::Array{UInt8,3},
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(u)
  # fields at initial time
  up = copy(u)
  ϕp = copy(ϕ)
  bp = copy(b)
  # loop over grid points
  for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      u[i,j,k] = (up[i-1,j,k] + up[i+1,j,k])/2 + 1/2c*(ϕp[i-1,j,k] - ϕp[i+1,j,k])
      ϕ[i,j,k] = (ϕp[i-1,j,k] + ϕp[i+1,j,k])/2 + c/2*(up[i-1,j,k] - up[i+1,j,k])
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    elseif maskx[i,j,k] == 2 # western boundary
      u[i,j,k] = 0.
      ϕ[i,j,k] = ϕp[i+1,j,k] - c*up[i+1,j,k]
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    elseif maskx[i,j,k] == 3 # eastern boundary
      u[i,j,k] = 0.
      ϕ[i,j,k] = ϕp[i-1,j,k] + c*up[i-1,j,k]
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(ϕ, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# sound wave split (y-direction)
@everywhere function Sy!(v::Array{Float64,3}, ϕ::Array{Float64,3}, b::Array{Float64,3}, masky::Array{UInt8,3},
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(v)
  # fields at initial time
  vp = copy(v)
  ϕp = copy(ϕ)
  bp = copy(b)
  # loop over grid points
  for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      v[i,j,k] = (vp[i,j-1,k] + vp[i,j+1,k])/2 + 1/2c*(ϕp[i,j-1,k] - ϕp[i,j+1,k])
      ϕ[i,j,k] = (ϕp[i,j-1,k] + ϕp[i,j+1,k])/2 + c/2*(vp[i,j-1,k] - vp[i,j+1,k])
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    elseif masky[i,j,k] == 2 # southern boundary
      v[i,j,k] = 0.
      ϕ[i,j,k] = ϕp[i,j+1,k] - c*vp[i,j+1,k]
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    elseif masky[i,j,k] == 3 # northern boundary
      v[i,j,k] = 0.
      ϕ[i,j,k] = ϕp[i,j-1,k] + c*vp[i,j-1,k]
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(ϕ, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(b, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# sound wave split (z-direction)
@everywhere function Sz!(w::Array{Float64,3}, ϕ::Array{Float64,3}, b::Array{Float64,3}, maskz::Array{UInt8,3},
                         chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(w)
  # fields at initial time
  wp = copy(w)
  ϕp = copy(ϕ)
  bp = copy(b)
  # loop over grid points
  for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      w[i,j,k] = (wp[i,j,k-1] + wp[i,j,k+1])/2 + μ/2c*(ϕp[i,j,k-1] - ϕp[i,j,k+1])
      ϕ[i,j,k] = (ϕp[i,j,k-1] + ϕp[i,j,k+1])/2 + c/2μ*(wp[i,j,k-1] - wp[i,j,k+1])
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    elseif maskz[i,j,k] == 2 # bottom boundary
      w[i,j,k] = 0.
      ϕ[i,j,k] = ϕp[i,j,k+1] - c/μ*wp[i,j,k+1]
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    elseif maskz[i,j,k] == 3 # top boundary
      w[i,j,k] = 0.
      ϕ[i,j,k] = ϕp[i,j,k-1] + c/μ*wp[i,j,k-1]
      b[i,j,k] = ϕp[i,j,k]*bp[i,j,k]/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(ϕ, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(b, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# buoyancy split (x-direction, gravity)
@everywhere function Tx!(u::Array{Float64,3}, ϕ::Array{Float64,3}, b::Array{Float64,3}, maskx::Array{UInt8,3},
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(u)
  # fields at midpoint
  ui = Array{Float64,3}(undef, nx, ny, nz)
  bi = Array{Float64,3}(undef, nx, ny, nz)
  # loop over grid points
  for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      ui[i,j,k] = u[i,j,k]
      bi[i,j,k] = b[i,j,k] + Δt/8Δx*((b[i-1,j,k] + b[i,j,k])*(u[i-1,j,k] + u[i,j,k])
                                     - (b[i,j,k] + b[i+1,j,k])*(u[i,j,k] + u[i+1,j,k]))*c^2/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 2) # western boundary
      ui[i,j,k] = 0.
      bi[i,j,k] = b[i,j,k] + Δt/4Δx*(-(b[i,j,k] + b[i+1,j,k])*u[i+1,j,k])*c^2/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 3) # eastern boundary
      ui[i,j,k] = 0.
      bi[i,j,k] = b[i,j,k] + Δt/4Δx*((b[i-1,j,k] + b[i,j,k])*u[i-1,j,k])*c^2/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangex!(ui, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(bi, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # loop over grid points
  for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      b[i,j,k] += Δt/4Δx*((bi[i-1,j,k] + bi[i,j,k])*(ui[i-1,j,k] + ui[i,j,k])
                          - (bi[i,j,k] + bi[i+1,j,k])*(ui[i,j,k] + ui[i+1,j,k]))*c^2/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 2) # western boundary
      u[i,j,k] = 0.
      b[i,j,k] += Δt/2Δx*(-(bi[i,j,k] + bi[i+1,j,k])*ui[i+1,j,k])*c^2/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 3) # eastern boundary
      u[i,j,k] = 0.
      b[i,j,k] += Δt/2Δx*((bi[i-1,j,k] + bi[i,j,k])*ui[i-1,j,k])*c^2/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# buoyancy split (y-direction, no gravity)
@everywhere function Ty!(v::Array{Float64,3}, ϕ::Array{Float64,3}, b::Array{Float64,3}, masky::Array{UInt8,3},
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(v)
  # fields at midpoint
  vi = Array{Float64,3}(undef, nx, ny, nz)
  bi = Array{Float64,3}(undef, nx, ny, nz)
  # loop over grid points
  for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      vi[i,j,k] = v[i,j,k]
      bi[i,j,k] = b[i,j,k] + Δt/8Δy*((b[i,j-1,k] + b[i,j,k])*(v[i,j-1,k] + v[i,j,k])
                                     - (b[i,j,k] + b[i,j+1,k])*(v[i,j,k] + v[i,j+1,k]))*c^2/ϕ[i,j,k]
    elseif (masky[i,j,k] == 2) # southern boundary
      vi[i,j,k] = 0.
      bi[i,j,k] = b[i,j,k] + Δt/4Δy*(-(b[i,j,k] + b[i,j+1,k])*v[i,j+1,k])*c^2/ϕ[i,j,k]
    elseif (masky[i,j,k] == 3) # northern boundary
      vi[i,j,k] = 0.
      bi[i,j,k] = b[i,j,k] + Δt/4Δy*((b[i,j-1,k] + b[i,j,k])*v[i,j-1,k])*c^2/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangey!(vi, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(bi, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # loop over grid points
  for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      b[i,j,k] += Δt/4Δy*((bi[i,j-1,k] + bi[i,j,k])*(vi[i,j-1,k] + vi[i,j,k])
                          - (bi[i,j,k] + bi[i,j+1,k])*(vi[i,j,k] + vi[i,j+1,k]))*c^2/ϕ[i,j,k]
    elseif (masky[i,j,k] == 2) # southern boundary
      v[i,j,k] = 0.
      b[i,j,k] += Δt/2Δy*(-(bi[i,j,k] + bi[i,j+1,k])*vi[i,j+1,k])*c^2/ϕ[i,j,k]
    elseif (masky[i,j,k] == 3) # northern boundary
      v[i,j,k] = 0.
      b[i,j,k] += Δt/2Δy*((bi[i,j-1,k] + bi[i,j,k])*vi[i,j-1,k])*c^2/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(b, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# buoyancy split (z-direction, gravity)
@everywhere function Tz!(w::Array{Float64,3}, ϕ::Array{Float64,3}, b::Array{Float64,3}, maskz::Array{UInt8,3},
                         chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}})
  # tile size
  nx, ny, nz = size(w)
  # fields at midpoint
  wi = Array{Float64,3}(undef, nx, ny, nz)
  bi = Array{Float64,3}(undef, nx, ny, nz)
  # loop over grid points
  for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      wi[i,j,k] = w[i,j,k] + Δt*μ^2/8*(b[i,j,k-1] + 2b[i,j,k] + b[i,j,k+1])
      bi[i,j,k] = b[i,j,k] + Δt/8Δz*((b[i,j,k-1] + b[i,j,k])*(w[i,j,k-1] + w[i,j,k])
                                     - (b[i,j,k] + b[i,j,k+1])*(w[i,j,k] + w[i,j,k+1]))*c^2/ϕ[i,j,k]
    elseif maskz[i,j,k] == 2 # bottom boundary
      wi[i,j,k] = 0.
      bi[i,j,k] = b[i,j,k] + Δt/4Δz*(-(b[i,j,k] + b[i,j,k+1])*w[i,j,k+1])*c^2/ϕ[i,j,k]
    elseif maskz[i,j,k] == 3 # top boundary
      wi[i,j,k] = 0.
      bi[i,j,k] = b[i,j,k] + Δt/4Δz*((b[i,j,k-1] + b[i,j,k])*w[i,j,k-1])*c^2/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangez!(wi, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(bi, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # loop over grid points
  for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      w[i,j,k] += Δt*μ^2/4*(bi[i,j,k-1] + 2bi[i,j,k] + bi[i,j,k+1])
      b[i,j,k] += Δt/4Δz*((bi[i,j,k-1] + bi[i,j,k])*(wi[i,j,k-1] + wi[i,j,k])
                          - (bi[i,j,k] + bi[i,j,k+1])*(wi[i,j,k] + wi[i,j,k+1]))*c^2/ϕ[i,j,k]
    elseif maskz[i,j,k] == 2 # bottom boundary
      w[i,j,k] = 0.
      b[i,j,k] += Δt/2Δz*(-(bi[i,j,k] + bi[i,j,k+1])*wi[i,j,k+1])*c^2/ϕ[i,j,k]
    elseif maskz[i,j,k] == 3 # top boundary
      w[i,j,k] = 0.
      b[i,j,k] += Δt/2Δz*((bi[i,j,k-1] + bi[i,j,k])*wi[i,j,k-1])*c^2/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(b, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# rotation split (in y–z plane)
@everywhere function Rx!(v::Array{Float64,3}, w::Array{Float64,3}, ϕ::Array{Float64,3}, maski::BitArray{3},
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
  for k = 2:nz-1, j = 2:ny-1, i = 1:nx
    if maski[i,j,k] # interior
      ωx = (wp[i,j+1,k] - wp[i,j-1,k])/2Δz - (vp[i,j,k+1] - vp[i,j,k-1])/2Δy # note switched Δ's
      γx = ωx*2Δt*c^2/ϕ[i,j,k]
      Cx = (1 - γx.^2/4)/(1 + γx^2/4)
      Sx = γx/(1 + γx^2/4)
      v[i,j,k] = vp[i,j,k]*Cx + wp[i,j,k]*Sx/μ
      w[i,j,k] = wp[i,j,k]*Cx - vp[i,j,k]*Sx*μ
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(w, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangez!(v, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# rotation split (in x–z plane)
@everywhere function Ry!(u::Array{Float64,3}, w::Array{Float64,3}, ϕ::Array{Float64,3}, maski::BitArray{3},
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
  for k = 2:nz-1, j = 1:ny, i = 2:nx-1
    if maski[i,j,k] # interior
      ωy = (up[i,j,k+1] - up[i,j,k-1])/2Δx - (wp[i+1,j,k] - wp[i-1,j,k])/2Δz # note switched Δ's
      γy = ωy*2Δt*c^2/ϕ[i,j,k]
      Cy = (1 - γy.^2/4)/(1 + γy^2/4)
      Sy = γy/(1 + γy^2/4)
      u[i,j,k] = up[i,j,k]*Cy - wp[i,j,k]*Sy/μ
      w[i,j,k] = wp[i,j,k]*Cy + up[i,j,k]*Sy*μ
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(w, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(u, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# rotation split (in x–y plane)
@everywhere function Rz!(u::Array{Float64,3}, v::Array{Float64,3}, ϕ::Array{Float64,3}, maski::BitArray{3},
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
  for k = 1:nz, j = 2:ny-1, i = 2:nx-1
    if maski[i,j,k] # interior
      ωz = f + (vp[i+1,j,k] - vp[i-1,j,k])/2Δx - (up[i,j+1,k] - up[i,j-1,k])/2Δy
      γz = ωz*2Δt*c^2/ϕ[i,j,k]
      Cz = (1 - γz.^2/4)/(1 + γz^2/4)
      Sz = γz/(1 + γz^2/4)
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

# x-diffusion
@everywhere function Dx!(a::Array{Float64,3}, maskx::Array{UInt8,3}, diri::BitArray{3},
                         chan_send_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_e::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_w::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_e::RemoteChannel{Channel{Array{Float64,2}}};
                         flux::Float64 = 0.)
  # tile size
  nx, ny, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      a[i,j,k] = (1-2α)*ap[i,j,k] + α*(ap[i-1,j,k] + ap[i+1,j,k])
    elseif maskx[i,j,k] == 2 # west boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i+1,j,k] + 4Δt*flux/Δx
      end
    elseif maskx[i,j,k] == 3 # east boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i-1,j,k] - 4Δt*flux/Δx
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# y-diffusion
@everywhere function Dy!(a::Array{Float64,3}, masky::Array{UInt8,3}, diri::BitArray{3},
                         chan_send_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_n::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_s::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_n::RemoteChannel{Channel{Array{Float64,2}}};
                         flux::Float64 = 0.)
  # tile size
  nx, ny, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      a[i,j,k] = (1-2α)*ap[i,j,k] + α*(ap[i,j-1,k] + ap[i,j+1,k])
    elseif masky[i,j,k] == 2 # west boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j+1,k] + 4Δt*flux/Δy
      end
    elseif masky[i,j,k] == 3 # east boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j-1,k] - 4Δt*flux/Δy
      end
    end
  end
  # exchange with neighboring tiles
  exchangey!(a, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# z-diffusion (definition of α takes enhancement by aspect ratio into account)
@everywhere function Dz!(a::Array{Float64,3}, maskz::Array{UInt8,3}, diri::BitArray{3},
                         chan_send_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_send_t::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_b::RemoteChannel{Channel{Array{Float64,2}}},
                         chan_receive_t::RemoteChannel{Channel{Array{Float64,2}}};
                         flux::Float64 = 0.)
  # tile size
  nx, ny, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      a[i,j,k] = (1-2α)*ap[i,j,k] + α*(ap[i,j,k-1] + ap[i,j,k+1])
    elseif maskz[i,j,k] == 2 # bottom boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j,k+1] + 4Δt*flux/Δz
      end
    elseif maskz[i,j,k] == 3 # top boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j,k-1]# - 4Δt*flux/Δz
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
  nx = 256; ny = 256; nz = 256
  # read abyssal hill topography from file
  b = reshape(h5read("abyssal_256x256.h5", "b"), (nx, ny, 1))
  z = [(k-1)*Δz for i = 1:nx, j = 1:ny, k = 1:nz]
  # assign to padded array
  fluid = BitArray(undef, nx+2, ny+2, nz+2)
  fluid[2:end-1,2:end-1,2:end-1] = z .> b .- minimum(b)
  # add solid boundary at the top
  fluid[:,:,end-1] .= false
  # fill edges
  exchange_global!(fluid)
  return fluid
end

# identify interior points (see Fig. 1 in Salmon for a sketch)
@everywhere function interior_mask(fluid::BitArray{3})
  # domain size
  nx, ny, nz = size(fluid)
  # find interior points
  maski = BitArray(undef, nx, ny, nz)
  for k = 2:nz-1, j = 2:ny-1, i = 2:nx-1
    if all(fluid[i-1:i,j-1:j,k-1:k]) # interior
      maski[i,j,k] = true
    else # not interior
      maski[i,j,k] = false
    end
  end
  exchange_global!(maski)
  return maski
end

# boundary masks (see Fig. 1 in Salmon for a sketch)
@everywhere function boundary_masks(maski::BitArray{3})
  # tile size
  nx, ny, nz = size(maski)
  # find and categorize zonal boundary points
  maskx = Array{UInt8,3}(undef, nx, ny, nz)
  masky = Array{UInt8,3}(undef, nx, ny, nz)
  maskz = Array{UInt8,3}(undef, nx, ny, nz)
  for k = 2:nz-1, j = 2:ny-1, i = 2:nx-1
    if maski[i,j,k] # interior
      maskx[i,j,k] = 1
      masky[i,j,k] = 1
      maskz[i,j,k] = 1
    else # not interior
      if maski[i+1,j,k] # western boundary
        maskx[i,j,k] = 2
      elseif maski[i-1,j,k] # eastern boundary
        maskx[i,j,k] = 3
      else # inside corner or solid
        maskx[i,j,k] = 0
      end
      if maski[i,j+1,k] # southern boundary
        masky[i,j,k] = 2
      elseif maski[i,j-1,k] # northern boundary
        masky[i,j,k] = 3
      else # inside corner or solid
        masky[i,j,k] = 0
      end
      if maski[i,j,k+1] # bottom boundary
        maskz[i,j,k] = 2
      elseif maski[i,j,k-1] # top boundary
        maskz[i,j,k] = 3
      else # inside corner or solid
        maskz[i,j,k] = 0
      end
    end
  end
  # exchange edges
  exchange_global!(maskx)
  exchange_global!(masky)
  exchange_global!(maskz)
  return maskx, masky, maskz
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
  solid = ((maskx .== 0) .& (masky .== 0) .& (maskz .== 0))[2:nx-1,2:ny-1,2:nz-1]
  us[solid] .= NaN
  vs[solid] .= NaN
  ws[solid] .= NaN
  ϕs[solid] .= NaN
  bs[solid] .= NaN
  # save data
  filename = @sprintf("data/%010d_%1d_%1d_%1d.h5", n, i, j, k)
  h5write(filename, "u", us)
  h5write(filename, "v", vs)
  h5write(filename, "w", ws)
  h5write(filename, "ϕ", ϕs)
  h5write(filename, "b", bs)
  # screen print
  println(@sprintf("%s %6i %9.3e %9.3e %9.3e %9.3e", floor(now(), Second), n, n*4Δt,
                   maximum(abs.(us[.!solid])), maximum(abs.(vs[.!solid])), maximum(abs.(ws[.!solid]))))
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
  filename = @sprintf("data/%010d_%1d_%1d_%1d.h5", n, i, j, k)
  u[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "u")
  v[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "v")
  w[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "w")
  ϕ[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "ϕ")
  b[2:nx-1,2:ny-1,2:nz-1] = h5read(filename, "b")
  # replace missing values with zeros
  solid = (maskx .== 0) .& (masky .== 0) .& (maskz .== 0)
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
  println(@sprintf("%s %6i %9.3e %9.3e %9.3e %9.3e", floor(now(), Second), n, n*4Δt,
                   maximum(abs.(u[.!solid])), maximum(abs.(v[.!solid])), maximum(abs.(w[.!solid]))))
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
  maski = interior_mask(fluid)
  maskx, masky, maskz = boundary_masks(maski)
  # select range of tile
  maski = maski[irange[1]:irange[end]+2,jrange[1]:jrange[end]+2,krange[1]:krange[end]+2]
  maskx = maskx[irange[1]:irange[end]+2,jrange[1]:jrange[end]+2,krange[1]:krange[end]+2]
  masky = masky[irange[1]:irange[end]+2,jrange[1]:jrange[end]+2,krange[1]:krange[end]+2]
  maskz = maskz[irange[1]:irange[end]+2,jrange[1]:jrange[end]+2,krange[1]:krange[end]+2]
  # specify where Dirichlet boundary conditions are to be imposed
  diriu = .!maski
  diriv = .!maski
  diriw = .!maski
  dirib = falses(nx, ny, nz)
  # coordinate
  z = [(k-1)*Δz for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1, k = krange[1]-1:krange[end]+1]
  if steps[1] == 1 # initialize
    # save masks
    h5write(@sprintf("data/masks_%1d_%1d_%1d.h5", i, j, k), "maskx", maskx[2:nx-1,2:ny-1,2:nz-1])
    h5write(@sprintf("data/masks_%1d_%1d_%1d.h5", i, j, k), "masky", masky[2:nx-1,2:ny-1,2:nz-1])
    h5write(@sprintf("data/masks_%1d_%1d_%1d.h5", i, j, k), "maskz", maskz[2:nx-1,2:ny-1,2:nz-1])
    # initial conditions
    u = zeros(nx, ny, nz)
    v = zeros(nx, ny, nz)
    w = zeros(nx, ny, nz)
    N = 1e-3
    b = N^2*z
    ϕ = c^2 .+ 1/2*N^2*z.^2
#    exchangex!(ϕ, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
#    exchangey!(ϕ, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
#    exchangez!(ϕ, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
#    exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
#    exchangey!(b, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
#    exchangez!(b, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    # save initial conditions to file
    save_tile(0, i, j, k, u, v, w, ϕ, b, maskx, masky, maskz)
  else # load
    u, v, w, ϕ, b = load_tile(steps[1] - 1, i, j, k, maskx, masky, maskz, chan_send_w, chan_send_e, chan_send_s, chan_send_n, chan_send_b, chan_send_t, chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
  end
  # time steps
  for n = steps
    # Strang splitting
    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(v, maskx, diriv, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(w, maskx, diriw, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)

    Dy!(u, masky, diriu, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(v, masky, diriv, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(w, masky, diriw, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(b, masky, dirib, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)

    Dz!(u, maskz, diriu, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(v, maskz, diriv, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(w, maskz, diriw, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(b, maskz, dirib, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)

    Rx!(v, w, ϕ, maski, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
        chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
    Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
        chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
        chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)

    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)

    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)

    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)

    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)

    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)

    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)

    Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
        chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
    Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
        chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    Rx!(v, w, ϕ, maski, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
        chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)

    Dz!(u, maskz, diriu, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(v, maskz, diriv, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(w, maskz, diriw, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(b, maskz, dirib, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)

    Dy!(u, masky, diriu, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(v, masky, diriv, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(w, masky, diriw, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(b, masky, dirib, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)

    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(v, maskx, diriv, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(w, maskx, diriw, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)

#    println(mass(ϕ[2:nx-1,2:ny-1,2:nz-1],
#                 maski[2:nx-1,2:ny-1,2:nz-1],
#                 maskx[2:nx-1,2:ny-1,2:nz-1],
#                 masky[2:nx-1,2:ny-1,2:nz-1],
#                 maskz[2:nx-1,2:ny-1,2:nz-1]))
#    println(energy(u[2:nx-1,2:ny-1,2:nz-1],
#                   v[2:nx-1,2:ny-1,2:nz-1],
#                   w[2:nx-1,2:ny-1,2:nz-1],
#                   ϕ[2:nx-1,2:ny-1,2:nz-1],
#                   b[2:nx-1,2:ny-1,2:nz-1],
#                   z[2:nx-1,2:ny-1,2:nz-1],
#                   maski[2:nx-1,2:ny-1,2:nz-1],
#                   maskx[2:nx-1,2:ny-1,2:nz-1],
#                   masky[2:nx-1,2:ny-1,2:nz-1],
#                   maskz[2:nx-1,2:ny-1,2:nz-1]))
#    println(buoyancy(ϕ[2:nx-1,2:ny-1,2:nz-1],
#                     b[2:nx-1,2:ny-1,2:nz-1],
#                     maski[2:nx-1,2:ny-1,2:nz-1],
#                     maskx[2:nx-1,2:ny-1,2:nz-1],
#                     masky[2:nx-1,2:ny-1,2:nz-1],
#                     maskz[2:nx-1,2:ny-1,2:nz-1]))
    if n % 100 == 0
      # save tile to file
      save_tile(n, i, j, k, u, v, w, ϕ, b, maskx, masky, maskz)
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
      iw = mod1(i - 1, ni); ie = mod1(i + 1, ni)
      js = mod1(j - 1, nj); jn = mod1(j + 1, nj)
      kb = mod1(k - 1, nk); kt = mod1(k + 1, nk)
      # run model on tile
      @async remotecall_fetch(run_tile, proc[i,j,k], i, j, k, irange, jrange, krange, steps, chan_w[i,j,k], chan_e[i,j,k],
                              chan_s[i,j,k], chan_n[i,j,k], chan_b[i,j,k], chan_t[i,j,k], chan_e[iw,j,k], chan_w[ie,j,k],
                              chan_n[i,js,k], chan_s[i,jn,k], chan_t[i,j,kb], chan_b[i,j,kt])
    end
  end
end

# mass
@everywhere function mass(ϕ, maski, maskx, masky, maskz)
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return sum(ϕ[maski])/c^2 + sum(ϕ[maskb])/2c^2
end

# energy
@everywhere function energy(u, v, w, ϕ, b, z, maski, maskx, masky, maskz)
  α = ϕ.*b/c^2
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return sum(u[maski].^2)/2 + sum(v[maski].^2)/2 + sum(w[maski].^2)/2μ^2 + sum(ϕ[maski].^2)/2c^2 + sum(ϕ[maskb].^2)/4c^2 - sum(α[maski].*z[maski]) - sum(α[maskb].*z[maskb])/2 - c^2*sum(maski)/2 - c^2*sum(maskb)/4
end

# buoyancy
@everywhere function buoyancy(ϕ, b, maski, maskx, masky, maskz)
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return sum(ϕ[maski].*b[maski]) + sum(ϕ[maskb].*b[maskb])/2
end

# buoyancy variance
@everywhere function buoyancy_variance(ϕ, b, maski, maskx, masky, maskz)
  maskb = (maskx .> 1) .| (masky .> 1) .| (maskz .> 1)
  return sum(ϕ[maski].*b[maski].^2) + sum(ϕ[maskb].*b[maskb].^2)/2
end

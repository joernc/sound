# Implementation of Rick Salmon's "An Ocean Circulation Model Based on Operator-Splitting, Hamiltonian Brackets, and the Inclusion
# of Sound Waves" (JPO, 2009)

# next steps:
# - extend to third dimension
# - allow specification of nonzero BC
# - allow stress BC (need to modify rotational split)
# - organize channels better?
# - type annotations for better performance?

@everywhere using Printf
@everywhere using HDF5

# grid spacing
@everywhere const Δx = 60e3/1024
@everywhere const Δy = Δx
@everywhere const Δz = 1000/512

# vertical viscosity/diffusion
@everywhere const ν = 1e-3

# forcing amplitude and frequency
@everywhere const u0 = .025
@everywhere const ω = 1.4e-4

# inertial frequency
@everywhere const f = 1e-4

# background stratification
@everywhere const N = 1e-3

# slope angle
@everywhere const θ = 2e-3

# sound speed
@everywhere const c = 10.

# time step
@everywhere const Δt = Δx/c

# aspect ratio
@everywhere const μ = Δz/Δx

# constant for diffusion steps
@everywhere const α = ν*Δt/Δz^2

@everywhere function exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # send edges
  put!(chan_send_w, a[2,:,:])
  put!(chan_send_e, a[end-1,:,:])
  # receive edges
  a[1,:,:] = take!(chan_receive_w)
  a[end,:,:] = take!(chan_receive_e)
end

@everywhere function exchangey!(a, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # send edges
  put!(chan_send_s, a[:,2,:])
  put!(chan_send_n, a[:,end-1,:])
  # receive edges
  a[:,1,:] = take!(chan_receive_s)
  a[:,end,:] = take!(chan_receive_n)
end

@everywhere function exchangez!(a, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # send edges
  put!(chan_send_b, a[:,:,2])
  put!(chan_send_t, a[:,:,end-1])
  # receive edges
  a[:,:,1] = take!(chan_receive_b)
  a[:,:,end] = take!(chan_receive_t)
end

# sound wave split (x-direction)
@everywhere function Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
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
@everywhere function Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
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
@everywhere function Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
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
@everywhere function Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # tile size
  nx, ny, nz = size(u)
  # fields at initial time
  up = copy(u)
  bp = copy(b)
  # fields at midpoint
  ui = Array{Float64, 3}(undef, nx, ny, nz)
  bi = Array{Float64, 3}(undef, nx, ny, nz)
  # loop over grid points
  for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      ui[i,j,k] = up[i,j,k] + Δt/8*(bp[i-1,j,k] + 2bp[i,j,k] + bp[i+1,j,k])*sin(θ)
      bi[i,j,k] = bp[i,j,k] + (c/8*((bp[i-1,j,k] + bp[i,j,k])*(up[i-1,j,k] + up[i,j,k])
                                    - (bp[i,j,k] + bp[i+1,j,k])*(up[i,j,k] + up[i+1,j,k]))
                               - c^2*Δt/2*up[i,j,k]*N^2*sin(θ))/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 2) # western boundary
      ui[i,j,k] = 0.
      bi[i,j,k] = bp[i,j,k] + c/4*(-(bp[i,j,k] + bp[i+1,j,k])*up[i+1,j,k])/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 3) # eastern boundary
      ui[i,j,k] = 0.
      bi[i,j,k] = bp[i,j,k] + c/4*((bp[i-1,j,k] + bp[i,j,k])*up[i-1,j,k])/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangex!(ui, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(bi, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # loop over grid points
  for k = 1:nz, j = 1:ny, i = 2:nx-1
    if maskx[i,j,k] == 1 # interior
      u[i,j,k] = up[i,j,k] + Δt/4*(bi[i-1,j,k] + 2bi[i,j,k] + bi[i+1,j,k])*sin(θ)
      b[i,j,k] = bp[i,j,k] + (c/4*((bi[i-1,j,k] + bi[i,j,k])*(ui[i-1,j,k] + ui[i,j,k])
                                   - (bi[i,j,k] + bi[i+1,j,k])*(ui[i,j,k] + ui[i+1,j,k]))
                              - c^2*Δt*ui[i,j,k]*N^2*sin(θ))/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 2) # western boundary
      u[i,j,k] = 0.
      b[i,j,k] = bp[i,j,k] + c/2*(-(bi[i,j,k] + bi[i+1,j,k])*ui[i+1,j,k])/ϕ[i,j,k]
    elseif (maskx[i,j,k] == 3) # eastern boundary
      u[i,j,k] = 0.
      b[i,j,k] = bp[i,j,k] + c/2*((bi[i-1,j,k] + bi[i,j,k])*ui[i-1,j,k])/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# buoyancy split (y-direction, no gravity)
@everywhere function Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # tile size
  nx, ny, nz = size(v)
  # fields at initial time
  vp = copy(v)
  bp = copy(b)
  # fields at midpoint
  vi = Array{Float64, 3}(undef, nx, ny, nz)
  bi = Array{Float64, 3}(undef, nx, ny, nz)
  # loop over grid points
  for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      vi[i,j,k] = vp[i,j,k]
      bi[i,j,k] = bp[i,j,k] + c/8*((bp[i,j-1,k] + bp[i,j,k])*(vp[i,j-1,k] + vp[i,j,k])
                                   - (bp[i,j,k] + bp[i,j+1,k])*(vp[i,j,k] + vp[i,j+1,k]))/ϕ[i,j,k]
    elseif (masky[i,j,k] == 2) # southern boundary
      vi[i,j,k] = 0.
      bi[i,j,k] = bp[i,j,k] + c/4*(-(bp[i,j,k] + bp[i,j+1,k])*vp[i,j+1,k])/ϕ[i,j,k]
    elseif (masky[i,j,k] == 3) # northern boundary
      vi[i,j,k] = 0.
      bi[i,j,k] = bp[i,j,k] + c/4*((bp[i,j-1,k] + bp[i,j,k])*vp[i,j-1,k])/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangey!(vi, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(bi, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  # loop over grid points
  for k = 1:nz, j = 2:ny-1, i = 1:nx
    if masky[i,j,k] == 1 # interior
      v[i,j,k] = vp[i,j,k]
      b[i,j,k] = bp[i,j,k] + c/4*((bi[i,j-1,k] + bi[i,j,k])*(vi[i,j-1,k] + vi[i,j,k])
                                  - (bi[i,j,k] + bi[i,j+1,k])*(vi[i,j,k] + vi[i,j+1,k]))/ϕ[i,j,k]
    elseif (masky[i,j,k] == 2) # southern boundary
      v[i,j,k] = 0.
      b[i,j,k] = bp[i,j,k] + c/2*(-(bi[i,j,k] + bi[i,j+1,k])*vi[i,j+1,k])/ϕ[i,j,k]
    elseif (masky[i,j,k] == 3) # northern boundary
      v[i,j,k] = 0.
      b[i,j,k] = bp[i,j,k] + c/2*((bi[i,j-1,k] + bi[i,j,k])*vi[i,j-1,k])/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangey!(v, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(b, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# buoyancy split (z-direction, gravity)
@everywhere function Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # tile size
  nx, ny, nz = size(w)
  # fields at initial time
  wp = copy(w)
  bp = copy(b)
  # fields at midpoint
  wi = Array{Float64, 3}(undef, nx, ny, nz)
  bi = Array{Float64, 3}(undef, nx, ny, nz)
  # loop over grid points
  for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      wi[i,j,k] = wp[i,j,k] + μ^2*Δt/8*(bp[i,j,k-1] + 2bp[i,j,k] + bp[i,j,k+1])*cos(θ)
      bi[i,j,k] = bp[i,j,k] + (c/8μ*((bp[i,j,k-1] + bp[i,j,k])*(wp[i,j,k-1] + wp[i,j,k])
                                     - (bp[i,j,k] + bp[i,j,k+1])*(wp[i,j,k] + wp[i,j,k+1]))
                               - c^2*Δt/2*wp[i,j,k]*N^2*cos(θ))/ϕ[i,j,k]
    elseif maskz[i,j,k] == 2 # bottom boundary
      wi[i,j,k] = 0.
      bi[i,j,k] = bp[i,j,k] + c/4μ*(-(bp[i,j,k] + bp[i,j,k+1])*wp[i,j,k+1])/ϕ[i,j,k]
    elseif maskz[i,j,k] == 3 # top boundary
      wi[i,j,k] = 0.
      bi[i,j,k] = bp[i,j,k] + c/4μ*((bp[i,j,k-1] + bp[i,j,k])*wp[i,j,k-1])/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangez!(wi, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(bi, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # loop over grid points
  for k = 2:nz-1, j = 1:ny, i = 1:nx
    if maskz[i,j,k] == 1 # interior
      w[i,j,k] = wp[i,j,k] + μ^2*Δt/4*(bi[i,j,k-1] + 2bi[i,j,k] + bi[i,j,k+1])*cos(θ)
      b[i,j,k] = bp[i,j,k] + (c/4μ*((bi[i,j,k-1] + bi[i,j,k])*(wi[i,j,k-1] + wi[i,j,k])
                                    - (bi[i,j,k] + bi[i,j,k+1])*(wi[i,j,k] + wi[i,j,k+1]))
                              - c^2*Δt*wi[i,j,k]*N^2*cos(θ))/ϕ[i,j,k]
    elseif maskz[i,j,k] == 2 # bottom boundary
      w[i,j,k] = 0.
      b[i,j,k] = bp[i,j,k] + c/2μ*(-(bi[i,j,k] + bi[i,j,k+1])*wi[i,j,k+1])/ϕ[i,j,k]
    elseif maskz[i,j,k] == 3 # top boundary
      w[i,j,k] = 0.
      b[i,j,k] = bp[i,j,k] + c/2μ*((bi[i,j,k-1] + bi[i,j,k])*wi[i,j,k-1])/ϕ[i,j,k]
    end
  end
  # exchange with neighboring tiles
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(b, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# rotation split (in y–z plane)
@everywhere function Rx!(v, w, ϕ, maski, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
                         chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
  # tile size
  nx, ny, nz = size(v)
  # fields at initial time
  vp = copy(v)
  wp = copy(w)
  # loop over grid points
  for k = 2:nz-1, j = 2:ny-1, i = 1:nx
    if maski[i,j,k] # interior
      ωx = (wp[i,j+1,k] - wp[i,j-1,k])/2Δz - (vp[i,j,k+1] - vp[i,j,k-1])/2Δy # note switched Δ's
      γx = ωx*Δt*c^2/ϕ[i,j,k]
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
@everywhere function Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
                         chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
  # tile size
  nx, ny, nz = size(u)
  # fields at initial time
  up = copy(u)
  wp = copy(w)
  # loop over grid points
  for k = 2:nz-1, j = 1:ny, i = 2:nx-1
    if maski[i,j,k] # interior
      ωy = (up[i,j,k+1] - up[i,j,k-1])/2Δx - (wp[i+1,j,k] - wp[i-1,j,k])/2Δz # note switched Δ's
      γy = ωy*Δt*c^2/ϕ[i,j,k]
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
@everywhere function Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
                         chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
  # tile size
  nx, ny, nz = size(u)
  # fields at initial time
  up = copy(u)
  vp = copy(v)
  # loop over grid points
  for k = 1:nz, j = 2:ny-1, i = 2:nx-1
    if maski[i,j,k] # interior
      ωz = f + (vp[i+1,j,k] - vp[i-1,j,k])/2Δx - (up[i,j+1,k] - up[i,j-1,k])/2Δy
      γz = ωz*Δt*c^2/ϕ[i,j,k]
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

# apply forcing
@everywhere function F!(t, u)
  u .+= u0*ω*cos(ω*t)*Δt
end

# x-diffusion
@everywhere function Dx!(a, maskx, diri, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e; flux=0.)
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
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i+1,j,k] + 2Δt*flux/Δx
      end
    elseif maskx[i,j,k] == 3 # east boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i-1,j,k] - 2Δt*flux/Δx
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# y-diffusion
@everywhere function Dy!(a, masky, diri, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n; flux=0.)
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
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j+1,k] + 2Δt*flux/Δy
      end
    elseif masky[i,j,k] == 3 # east boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j-1,k] - 2Δt*flux/Δy
      end
    end
  end
  # exchange with neighboring tiles
  exchangey!(a, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
end

# z-diffusion (definition of α takes enhancement by aspect ratio into account)
@everywhere function Dz!(a, maskz, diri, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t; flux=0.)
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
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j,k+1] + 2Δt*flux/Δz
      end
    elseif maskz[i,j,k] == 3 # top boundary
      if diri[i,j,k] # Dirichlet BC
        a[i,j,k] = 0.
      else # Neumann BC
        a[i,j,k] = (1-2α)*ap[i,j,k] + 2α*ap[i,j,k-1]# - 2Δt*flux/Δz
      end
    end
  end
  # exchange with neighboring tiles
  exchangez!(a, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# boundary masks (see Fig. 1 in Salmon for a sketch)
@everywhere function boundary_masks(fluid, chan_send_w, chan_send_e, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
                                    chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
  nx, ny, nz = size(fluid)
  maski = Array{Bool, 3}(undef, nx, ny, nz)
  maskx = Array{UInt8, 3}(undef, nx, ny, nz)
  masky = Array{UInt8, 3}(undef, nx, ny, nz)
  maskz = Array{UInt8, 3}(undef, nx, ny, nz)
  for k = 2:nz-1, j = 2:ny-1, i = 2:nx-1
    if all(fluid[i-1:i,j-1:j,k-1:k]) # interior
      maski[i,j,k] = true
      maskx[i,j,k] = 1
      masky[i,j,k] = 1
      maskz[i,j,k] = 1
    else
      maski[i,j,k] = false
      if all(fluid[i,j-1:j,k-1:k]) # western boundary
        maskx[i,j,k] = 2
      elseif all(fluid[i-1,j-1:j,k-1:k]) # eastern boundary
        maskx[i,j,k] = 3
      else # inside corner or solid
        maskx[i,j,k] = 0
      end
      if all(fluid[i-1:i,j,k-1:k]) # southern boundary
        masky[i,j,k] = 2
      elseif all(fluid[i-1:i,j-1,k-1:k]) # northern boundary
        masky[i,j,k] = 3
      else # inside corner or solid
        masky[i,j,k] = 0
      end
      if all(fluid[i-1:i,j-1:j,k]) # bottom boundary
        maskz[i,j,k] = 2
      elseif all(fluid[i-1:i,j-1:j,k-1]) # top boundary
        maskz[i,j,k] = 3
      else # inside corner or solid
        maskz[i,j,k] = 0
      end
    end
  end
  # exchange edges
  exchangex!(maski, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(masky, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(maskz, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(maski, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(maskx, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangey!(maskz, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangez!(maski, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(maskx, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(masky, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  return maski, maskx, masky, maskz
end

# read in fluid mask
@everywhere function read_topo(irange, jrange, krange, chan_send_w, chan_send_e, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
                               chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
  # tile size plus two points padding
  nx = length(irange) + 2
  ny = length(jrange) + 2
  nz = length(krange) + 2
  # read abyssal hill topography from file
  b = reshape(h5read("abyssal.h5", "b"), (1024, 1, 1))
  z = [(k-1)*Δz for i = 1:1024, j = 1:5, k = 1:512]
  full_mask = z .> b .- minimum(b)
  # add solid boundary at the top
  full_mask[:,:,end] .= false
  # assign to tile
  fluid = Array{Bool, 3}(undef, nx, ny, nz)
  fluid[2:nx-1,2:ny-1,2:nz-1] = full_mask[irange,jrange,krange]
  # exchange edges
  exchangex!(fluid, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangey!(fluid, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
  exchangez!(fluid, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # interior, x-, and z-masks
  maski, maskx, masky, maskz = boundary_masks(fluid, chan_send_w, chan_send_e, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
                                              chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n, chan_receive_b,
                                              chan_receive_t)
  return maski, maskx, masky, maskz
end

# run model on tile
@everywhere function run_tile(i, j, k, irange, jrange, krange, steps, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
                              chan_send_b, chan_send_t, chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n,
                              chan_receive_b, chan_receive_t)
  # read topography mask
  maski, maskx, masky, maskz = read_topo(irange, jrange, krange, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
                                         chan_send_b, chan_send_t, chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n,
                                         chan_receive_b, chan_receive_t)
  # initialization
  nx = length(irange) + 2
  ny = length(jrange) + 2
  nz = length(krange) + 2
  u = zeros(nx, ny, nz)
  v = zeros(nx, ny, nz)
  w = zeros(nx, ny, nz)
  x = [(i-1)*Δx for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1, k = krange[1]-1:krange[end]+1]
  y = [(j-1)*Δx for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1, k = krange[1]-1:krange[end]+1]
  z = [(k-1)*Δz for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1, k = krange[1]-1:krange[end]+1]
  ϕ = c^2*ones(nx, ny, nz)
  b = zeros(nx, ny, nz)
  # specify where Dirichlet boundary conditions are to be imposed
  diriu = .!maski
  diriv = .!maski
  diriw = .!maski
  dirib = falses(nx, ny, nz)
  # time steps
  for n = 1:steps
    if n % 100 == 1
      println(@sprintf("%6i %9.3e %9.3e %9.3e %9.3e", n-1, (n-1)*2Δt, maximum(abs.(u[2:nx-1,2:ny-1,2:nz-1])),
                       maximum(abs.(v[2:nx-1,2:ny-1,2:nz-1])), maximum(abs.(w[2:nx-1,2:ny-1,2:nz-1]))))
      # discard edges
      us = u[2:nx-1,2:ny-1,2:nz-1]
      vs = v[2:nx-1,2:ny-1,2:nz-1]
      ws = w[2:nx-1,2:ny-1,2:nz-1]
      ϕs = ϕ[2:nx-1,2:ny-1,2:nz-1]
      bs = b[2:nx-1,2:ny-1,2:nz-1]# + N^2*(x[2:nx-1,2:ny-1,2:nz-1]*sin(θ) + z[2:nx-1,2:ny-1,2:nz-1]*cos(θ))
      # replace missing values with NaNs
      us[((maskx .== 0) .& (masky .== 0) .& (maskz .== 0))[2:nx-1,2:ny-1,2:nz-1]] .= NaN
      vs[((maskx .== 0) .& (masky .== 0) .& (maskz .== 0))[2:nx-1,2:ny-1,2:nz-1]] .= NaN
      ws[((maskx .== 0) .& (masky .== 0) .& (maskz .== 0))[2:nx-1,2:ny-1,2:nz-1]] .= NaN
      ϕs[((maskx .== 0) .& (masky .== 0) .& (maskz .== 0))[2:nx-1,2:ny-1,2:nz-1]] .= NaN
      bs[((maskx .== 0) .& (masky .== 0) .& (maskz .== 0))[2:nx-1,2:ny-1,2:nz-1]] .= NaN
      # save data
      filename = @sprintf("data/%010d_%1d_%1d_%1d.h5", n-1, i, j, k)
      h5write(filename, "u", us)
      h5write(filename, "v", vs)
      h5write(filename, "w", ws)
      h5write(filename, "ϕ", ϕs)
      h5write(filename, "b", bs)
    end
    # Strang splitting
    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(v, maskx, diriv, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(w, maskx, diriw, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e; flux=ν*N^2*sin(θ))
    Dy!(u, masky, diriu, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(v, masky, diriv, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(w, masky, diriw, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(b, masky, dirib, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dz!(u, maskz, diriu, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(v, maskz, diriv, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(w, maskz, diriw, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(b, maskz, dirib, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t; flux=ν*N^2*cos(θ))
    F!((2k-1.5)*Δt, u)
    Rx!(v, w, ϕ, maski, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
        chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
    Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
        chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
        chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Sy!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Ty!(v, ϕ, b, masky, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Rz!(u, v, ϕ, maski, chan_send_w, chan_send_e, chan_send_s, chan_send_n,
        chan_receive_w, chan_receive_e, chan_receive_s, chan_receive_n)
    Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
        chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    Rx!(v, w, ϕ, maski, chan_send_s, chan_send_n, chan_send_b, chan_send_t,
        chan_receive_s, chan_receive_n, chan_receive_b, chan_receive_t)
    F!((2k-.5)*Δt, u)
    Dz!(u, maskz, diriu, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(v, maskz, diriv, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(w, maskz, diriw, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(b, maskz, dirib, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t; flux=ν*N^2*cos(θ))
    Dy!(u, masky, diriu, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(v, masky, diriv, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(w, masky, diriw, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dy!(b, masky, dirib, chan_send_s, chan_send_n, chan_receive_s, chan_receive_n)
    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(v, maskx, diriv, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(w, maskx, diriw, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e; flux=ν*N^2*sin(θ))
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
  # set up remote channels for edge communication
  chan_w = [RemoteChannel(()->Channel{Array{Float64, 2}}(1)) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_e = [RemoteChannel(()->Channel{Array{Float64, 2}}(1)) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_s = [RemoteChannel(()->Channel{Array{Float64, 2}}(1)) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_n = [RemoteChannel(()->Channel{Array{Float64, 2}}(1)) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_b = [RemoteChannel(()->Channel{Array{Float64, 2}}(1)) for i = 1:ni, j = 1:nj, k = 1:nk]
  chan_t = [RemoteChannel(()->Channel{Array{Float64, 2}}(1)) for i = 1:ni, j = 1:nj, k = 1:nk]
  # spawn work on tiles
  a = Array{Future}(undef, ni, nj, nk)
  for i in 1:ni
    for j in 1:nj
      for k in 1:nk
        irange, jrange, krange = tile_range(i, j, k, tile_sizes)
        worker_no = mod1((k-1)*ni*nj + (j-1)*ni + i, nprocs() - 1)
        p = workers()[worker_no]
        println(@sprintf("worker %1d: tile %1d,%1d,%1d", p, i, j, k))
        println(irange); println(jrange); println(krange)
        # indices of neighboring tiles
        iw = mod1(i - 1, ni); ie = mod1(i + 1, ni)
        js = mod1(j - 1, nj); jn = mod1(j + 1, nj)
        kb = mod1(k - 1, nk); kt = mod1(k + 1, nk)
        # run on the tiles
        a[i,j,k] = @spawnat p run_tile(i, j, k, irange, jrange, krange, steps, chan_w[i,j,k], chan_e[i,j,k], chan_s[i,j,k],
                                       chan_n[i,j,k], chan_b[i,j,k], chan_t[i,j,k], chan_e[iw,j,k], chan_w[ie,j,k], chan_n[i,js,k],
                                       chan_s[i,jn,k], chan_t[i,j,kb], chan_b[i,j,kt])
      end
    end
  end
  # wait for results
  wait.(a)
end

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
@everywhere const Δz = 1000/512

# vertical viscosity/diffusion
@everywhere const ν = 1e-2

# forcing amplitude and frequency
@everywhere const u0 = .025
@everywhere const ω = 1.4e-4

# background stratification
@everywhere const N = 1e-3

# slope angle
@everywhere const θ = 2e-3

# sound speed
@everywhere const c = 1.

# time step
@everywhere const Δt = Δx/c

# aspect ratio
@everywhere const μ = Δz/Δx

# constant for diffusion steps
@everywhere const α = ν*Δt/Δz^2

@everywhere function exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # send edges
  put!(chan_send_w, a[2,:])
  put!(chan_send_e, a[end-1,:])
  # receive edges
  a[1,:] = take!(chan_receive_w)
  a[end,:] = take!(chan_receive_e)
end

@everywhere function exchangez!(a, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # send edges
  put!(chan_send_b, a[:,2])
  put!(chan_send_t, a[:,end-1])
  # receive edges
  a[:,1] = take!(chan_receive_b)
  a[:,end] = take!(chan_receive_t)
end

# sound wave split (x-direction)
@everywhere function Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # tile size
  nx, nz = size(u)
  # fields at initial time
  up = copy(u)
  ϕp = copy(ϕ)
  bp = copy(b)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:nz
      if maskx[i,j] == 1 # interior
        u[i,j] = (up[i-1,j] + up[i+1,j])/2 + 1/2c*(ϕp[i-1,j] - ϕp[i+1,j])
        ϕ[i,j] = (ϕp[i-1,j] + ϕp[i+1,j])/2 + c/2*(up[i-1,j] - up[i+1,j])
        b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif maskx[i,j] == 2 # western boundary
        u[i,j] = 0.
        ϕ[i,j] = ϕp[i+1,j] - c*up[i+1,j]
        b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif maskx[i,j] == 3 # eastern boundary
        u[i,j] = 0.
        ϕ[i,j] = ϕp[i-1,j] + c*up[i-1,j]
        b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(ϕ, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# sound wave split (z-direction)
@everywhere function Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # tile size
  nx, nz = size(w)
  # fields at initial time
  wp = copy(w)
  ϕp = copy(ϕ)
  bp = copy(b)
  # loop over grid points
  for i = 1:nx
    for j = 2:nz-1
      if maskz[i,j] == 1 # interior
        w[i,j] = (wp[i,j-1] + wp[i,j+1])/2 + μ/2c*(ϕp[i,j-1] - ϕp[i,j+1])
        ϕ[i,j] = (ϕp[i,j-1] + ϕp[i,j+1])/2 + c/2μ*(wp[i,j-1] - wp[i,j+1])
        b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif maskz[i,j] == 2 # bottom boundary
        w[i,j] = 0.
        ϕ[i,j] = ϕp[i,j+1] - c/μ*wp[i,j+1]
        b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      elseif maskz[i,j] == 3 # top boundary
        w[i,j] = 0.
        ϕ[i,j] = ϕp[i,j-1] + c/μ*wp[i,j-1]
        b[i,j] = ϕp[i,j]*bp[i,j]/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(ϕ, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(b, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# buoyancy split (x-direction, no gravity)
@everywhere function Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # tile size
  nx, nz = size(u)
  # fields at initial time
  up = copy(u)
  bp = copy(b)
  # fields at midpoint
  ui = Array{Float64, 2}(undef, nx, nz)
  bi = Array{Float64, 2}(undef, nx, nz)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:nz
      if maskx[i,j] == 1 # interior
        ui[i,j] = up[i,j] + Δt/8*(bp[i-1,j] + 2bp[i,j] + bp[i+1,j])*sin(θ)
        bi[i,j] = bp[i,j] + (c/8*((bp[i-1,j] + bp[i,j])*(up[i-1,j] + up[i,j]) - (bp[i,j] + bp[i+1,j])*(up[i,j] + up[i+1,j]))
                             - Δt/2*up[i,j]*N^2*sin(θ))/ϕ[i,j]
      elseif (maskx[i,j] == 2) # western boundary
        ui[i,j] = 0.
        bi[i,j] = bp[i,j] + c/4*(-(bp[i,j] + bp[i+1,j])*up[i+1,j])/ϕ[i,j]
      elseif (maskx[i,j] == 3) # eastern boundary
        ui[i,j] = 0.
        bi[i,j] = bp[i,j] + c/4*((bp[i-1,j] + bp[i,j])*up[i-1,j])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(ui, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(bi, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:nz
      if maskx[i,j] == 1 # interior
        u[i,j] = up[i,j] + Δt/4*(bi[i-1,j] + 2bi[i,j] + bi[i+1,j])*sin(θ)
        b[i,j] = bp[i,j] + (c/4*((bi[i-1,j] + bi[i,j])*(ui[i-1,j] + ui[i,j]) - (bi[i,j] + bi[i+1,j])*(ui[i,j] + ui[i+1,j]))
                            - Δt*ui[i,j]*N^2*sin(θ))/ϕ[i,j]
      elseif (maskx[i,j] == 2) # west boundary
        u[i,j] = 0.
        b[i,j] = bp[i,j] + c/2*(-(bi[i,j] + bi[i+1,j])*ui[i+1,j])/ϕ[i,j]
      elseif (maskx[i,j] == 3) # east boundary
        u[i,j] = 0.
        b[i,j] = bp[i,j] + c/2*((bi[i-1,j] + bi[i,j])*ui[i-1,j])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangex!(b, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# buoyancy split (z-direction, gravity)
@everywhere function Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # tile size
  nx, nz = size(w)
  # fields at initial time
  wp = copy(w)
  bp = copy(b)
  # fields at midpoint
  wi = Array{Float64, 2}(undef, nx, nz)
  bi = Array{Float64, 2}(undef, nx, nz)
  # loop over grid points
  for i = 1:nx
    for j = 2:nz-1
      if maskz[i,j] == 1 # interior
        wi[i,j] = wp[i,j] + μ^2*Δt/8*(bp[i,j-1] + 2bp[i,j] + bp[i,j+1])*cos(θ)
        bi[i,j] = bp[i,j] + (c/8μ*((bp[i,j-1] + bp[i,j])*(wp[i,j-1] + wp[i,j]) - (bp[i,j] + bp[i,j+1])*(wp[i,j] + wp[i,j+1]))
                             - Δt/2*wp[i,j]*N^2*cos(θ))/ϕ[i,j]
      elseif maskz[i,j] == 2 # bottom boundary
        wi[i,j] = 0.
        bi[i,j] = bp[i,j] + c/4μ*(-(bp[i,j] + bp[i,j+1])*wp[i,j+1])/ϕ[i,j]
      elseif maskz[i,j] == 3 # top boundary
        wi[i,j] = 0.
        bi[i,j] = bp[i,j] + c/4μ*((bp[i,j-1] + bp[i,j])*wp[i,j-1])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangez!(wi, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(bi, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # loop over grid points
  for i = 1:nx
    for j = 2:nz-1
      if maskz[i,j] == 1 # interior
        w[i,j] = wp[i,j] + μ^2*Δt/4*(bi[i,j-1] + 2bi[i,j] + bi[i,j+1])*cos(θ)
        b[i,j] = bp[i,j] + (c/4μ*((bi[i,j-1] + bi[i,j])*(wi[i,j-1] + wi[i,j]) - (bi[i,j] + bi[i,j+1])*(wi[i,j] + wi[i,j+1]))
                            - Δt*wi[i,j]*N^2*cos(θ))/ϕ[i,j]
      elseif maskz[i,j] == 2 # bottom boundary
        w[i,j] = 0.
        b[i,j] = bp[i,j] + c/2μ*(-(bi[i,j] + bi[i,j+1])*wi[i,j+1])/ϕ[i,j]
      elseif maskz[i,j] == 3 # top boundary
        w[i,j] = 0.
        b[i,j] = bp[i,j] + c/2μ*((bi[i,j-1] + bi[i,j])*wi[i,j-1])/ϕ[i,j]
      end
    end
  end
  # exchange with neighboring tiles
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangez!(b, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# rotation split (in x–z plane)
@everywhere function Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
                         chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
  # tile size
  nx, nz = size(u)
  # fields at initial time
  up = copy(u)
  wp = copy(w)
  # loop over grid points
  for i = 2:nx-1
    for j = 2:nz-1
      if maski[i,j] # interior
        ωy = (up[i,j+1] - up[i,j-1])/2Δx - (wp[i+1,j] + wp[i-1,j])/2Δz # note switched Δ's
        γy = ωy*Δt*c^2/ϕ[i,j]
        Cy = (1 - γy.^2/4)/(1 + γy^2/4)
        Sy = γy/(1 + γy^2/4)
        u[i,j] = up[i,j]*Cy - wp[i,j]*Sy/μ
        w[i,j] = wp[i,j]*Cy + up[i,j]*Sy*μ
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(u, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(u, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangex!(w, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(w, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# apply forcing
@everywhere function F!(t, u)
  u .+= u0*ω*cos(ω*t)*Δt
end

# x-diffusion
@everywhere function Dx!(a, maskx, diri, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e; flux=0.)
  # tile size
  nx, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  for i = 2:nx-1
    for j = 1:nz
      if maskx[i,j] == 1 # interior
        a[i,j] = (1-2α)*ap[i,j] + α*(ap[i-1,j] + ap[i+1,j])
      elseif maskx[i,j] == 2 # west boundary
        if diri[i,j] # Dirichlet BC
          a[i,j] = 0.
        else # Neumann BC
          a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i+1,j] + 2Δt*flux/Δx
        end
      elseif maskx[i,j] == 3 # east boundary
        if diri[i,j] # Dirichlet BC
          a[i,j] = 0.
        else # Neumann BC
          a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i-1,j] - 2Δt*flux/Δx
        end
      end
    end
  end
  # exchange with neighboring tiles
  exchangex!(a, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
end

# z-diffusion (definition of α takes enhancement by aspect ratio into account)
@everywhere function Dz!(a, maskz, diri, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t; flux=0.)
  # tile size
  nx, nz = size(a)
  # field at initial time
  ap = copy(a)
  # loop over grid points
  for i = 1:nx
    for j = 2:nz-1
      if maskz[i,j] == 1 # interior
        a[i,j] = (1-2α)*ap[i,j] + α*(ap[i,j-1] + ap[i,j+1])
      elseif maskz[i,j] == 2 # bottom boundary
        if diri[i,j] # Dirichlet BC
          a[i,j] = 0.
        else # Neumann BC
          a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i,j+1] + 2Δt*flux/Δz
        end
      elseif maskz[i,j] == 3 # top boundary
        if diri[i,j] # Dirichlet BC
          a[i,j] = 0.
        else # Neumann BC
          a[i,j] = (1-2α)*ap[i,j] + 2α*ap[i,j-1] - 2Δt*flux/Δz
        end
      end
    end
  end
  # exchange with neighboring tiles
  exchangez!(a, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
end

# boundary masks (see Fig. 1 in Salmon for a sketch)
@everywhere function boundary_masks(fluid, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
                                    chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
  nx, nz = size(fluid)
  maskx = zeros(UInt8, nx, nz)
  maskz = zeros(UInt8, nx, nz)
  for i = 2:nx-1
    for j = 2:nz-1
      # east and west boundaries
      if fluid[i-1,j-1] & fluid[i,j-1] & fluid[i-1,j] & fluid[i,j] # interior
        maskx[i,j] = 1
      elseif (!fluid[i-1,j-1] | !fluid[i-1,j]) & fluid[i,j-1] & fluid[i,j] # western boundary
        maskx[i,j] = 2
      elseif (!fluid[i,j-1] | !fluid[i,j]) & fluid[i-1,j-1] & fluid[i-1,j] # eastern boundary
        maskx[i,j] = 3
      end
      # top and bottom boundaries
      if fluid[i-1,j-1] & fluid[i-1,j] & fluid[i,j-1] & fluid[i,j] # interior
        maskz[i,j] = 1
      elseif (!fluid[i-1,j-1] | !fluid[i,j-1]) & fluid[i-1,j] & fluid[i,j] # bottom boundary
        maskz[i,j] = 2
      elseif (!fluid[i-1,j] | !fluid[i,j]) & fluid[i-1,j-1] & fluid[i,j-1] # top boundary
        maskz[i,j] = 3
      end
    end
  end
  # interior points
  maski = (maskx .== 1) .& (maskz .== 1)
  # exchange edges
  exchangex!(maski, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(maski, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangex!(maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(maskx, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  exchangex!(maskz, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  return maski, maskx, maskz
end

# read in fluid mask
@everywhere function read_topo(irange, jrange, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
                               chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
  # tile size plus two points padding
  nx = length(irange) + 2
  nz = length(jrange) + 2
  # read abyssal hill topography from file
  b = reshape(h5read("abyssal.h5", "b"), (1024, 1))
  z = [(j-1)*Δz for i = 1:1024, j = 1:512]
  full_mask = z .> b .- minimum(b)
  # add solid boundary at the top
  full_mask[:,end] .= false
  # assign to tile
  fluid = Array{Bool}(undef, nx, nz)
  fluid[2:nx-1,2:nz-1] = full_mask[irange,jrange]
  # exchange edges
  exchangex!(fluid, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
  exchangez!(fluid, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
  # interior, x-, and y-masks
  maski, maskx, maskz = boundary_masks(fluid, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
                                       chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
  return maski, maskx, maskz
end

# run model on tile
@everywhere function run_tile(i, j, irange, jrange, steps, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
                              chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
  # read topography mask
  maski, maskx, maskz = read_topo(irange, jrange, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
                                  chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
  # initialization
  nx = length(irange) + 2
  nz = length(jrange) + 2
  u = zeros(nx, nz)
  w = zeros(nx, nz)
  x = [(i-1)*Δx for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1]
  z = [(j-1)*Δz for i = irange[1]-1:irange[end]+1, j = jrange[1]-1:jrange[end]+1]
  ϕ = c^2*ones(nx, nz)
  b = zeros(nx, nz)
  # specify where Dirichlet boundary conditions are to be imposed
  diriu = (maskx .> 1) .| (maskz .> 1)
  diriw = (maskx .> 1) .| (maskz .> 1)
  dirib = falses(nx, nz)
  # time steps
  for k = 1:steps
    if k % 100 == 1
      println(@sprintf("%6i %9.3e %9.3e %9.3e", k-1, (k-1)*2Δt, maximum(abs.(u[2:nx-1,2:nz-1])), maximum(abs.(w[2:nx-1,2:nz-1]))))
      # discard edges
      us = u[2:nx-1,2:nz-1]
      ws = w[2:nx-1,2:nz-1]
      ϕs = ϕ[2:nx-1,2:nz-1]
      bs = b[2:nx-1,2:nz-1] + N^2*(x[2:nx-1,2:nz-1]*sin(θ) + z[2:nx-1,2:nz-1]*cos(θ))
      # replace missing values with NaNs
      us[((maskx .== 0) .& (maskz .== 0))[2:nx-1, 2:nz-1]] .= NaN
      ws[((maskx .== 0) .& (maskz .== 0))[2:nx-1, 2:nz-1]] .= NaN
      ϕs[((maskx .== 0) .& (maskz .== 0))[2:nx-1, 2:nz-1]] .= NaN
      bs[((maskx .== 0) .& (maskz .== 0))[2:nx-1, 2:nz-1]] .= NaN
      # save data
      filename = @sprintf("data/%010d_%1d_%1d.h5", k-1, i, j)
      h5write(filename, "u", us)
      h5write(filename, "w", ws)
      h5write(filename, "ϕ", ϕs)
      h5write(filename, "b", bs)
    end
    # Strang splitting
    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(w, maskx, diriw, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e; flux=ν*N^2*sin(θ))
    Dz!(u, maskz, diriu, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(w, maskz, diriw, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(b, maskz, dirib, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t; flux=ν*N^2*cos(θ))
    #F!((2k-1.5)*Δt, u)
    Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
        chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Tz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Sz!(w, ϕ, b, maskz, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Sx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Tx!(u, ϕ, b, maskx, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Ry!(u, w, ϕ, maski, chan_send_w, chan_send_e, chan_send_b, chan_send_t,
        chan_receive_w, chan_receive_e, chan_receive_b, chan_receive_t)
    #F!((2k-.5)*Δt, u)
    Dz!(u, maskz, diriu, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(w, maskz, diriw, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t)
    Dz!(b, maskz, dirib, chan_send_b, chan_send_t, chan_receive_b, chan_receive_t; flux=ν*N^2*cos(θ))
    Dx!(u, maskx, diriu, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(w, maskx, diriw, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e)
    Dx!(b, maskx, dirib, chan_send_w, chan_send_e, chan_receive_w, chan_receive_e; flux=ν*N^2*sin(θ))
  end
end

# get index ranges for tiles
tile_range(i, j, tile_sizes) = [sum(tile_sizes[1][1:i-1])+1:sum(tile_sizes[1][1:i]),
                                sum(tile_sizes[2][1:j-1])+1:sum(tile_sizes[2][1:j])]

# run the model
function run_model(steps, tile_sizes)
  # number of tiles
  ni = length(tile_sizes[1])
  nj = length(tile_sizes[2])
  # set up remote channels for edge communication
  chan_w = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  chan_e = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  chan_b = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  chan_t = [RemoteChannel(()->Channel{Array{Float64, 1}}(1)) for i = 1:ni, j = 1:nj]
  # spawn work on tiles
  a = Array{Future}(undef, ni, nj)
  for i in 1:ni
    for j in 1:nj
      irange, jrange = tile_range(i, j, tile_sizes)
      worker_no = mod1((i-1)*nj + j, nprocs() - 1)
      p = workers()[worker_no]
      println(p, ": ", i, " ", j, " ", irange, " ", jrange)
      a[i,j] = @spawnat p run_tile(i, j, irange, jrange, steps, chan_w[i,j], chan_e[i,j], chan_b[i,j], chan_t[i,j],
                                   chan_e[mod1(i-1,ni),j], chan_w[mod1(i+1,ni),j], chan_t[i,mod1(j-1,nj)], chan_b[i,mod1(j+1,nj)])
    end
  end
  # wait for results
  wait.(a)
end

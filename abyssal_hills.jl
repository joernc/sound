using Random
using FFTW
using HDF5
using Printf

# domain size
Lx = 32768.
Ly = 32768.
Lz = 1024.

# number of grid points
nx = 256
ny = 256
nz = 256

Δz = Lz/nz

# topography parameters
k0 = 1.0e-3
l0 = 2.2e-4
hm = 110.
ν = .9

# diffusion parameters
ν0 = 5.2e-5
ν1 = 1.8e-3
d = 230.

# 1D spectrum (adapted from Nikurashin and Legg, 2011, eq. 1)
H(k, k0, hm, ν) = 2π*ν*hm^2*gamma(ν+1/2)*k0^2ν/(sqrt(π)*gamma(ν+1))*(k^2 + k0^2)^(-(ν+1/2))

# 2D spectrum (adapted from Nikurashin and Legg, 2011, eq. 1)
H(k, l, k0, l0, hm, ν) = 4π*ν*hm^2/(k0*l0)*((k/k0)^2 + (l/l0)^2 + 1)^(-(ν+1))

# grid spacing
Δx = Lx/nx
Δy = Ly/ny

# wavenumber vector
Δk = 1/Lx
Δl = 1/Ly
k = reshape(2π*(0:Δk:1/2Δx), (nx÷2+1, 1))
l = reshape(fftshift(2π*(-1/2Δy:Δl:1/2Δy-Δl)), (1, ny))

# spectral coefficients (random phase)
Random.seed!(42)
S = sqrt.(H.(k, l, k0, l0, hm, ν)).*exp.(2π*1im*rand(nx÷2+1, ny))

# bottom height
h = irfft(S, nx)*nx*ny/sqrt(Lx*Ly)
h .-= minimum(h)

# slopes
hx = irfft(1im*k.*S, nx)*nx*ny/sqrt(Lx*Ly)
hy = irfft(1im*l.*S, nx)*nx*ny/sqrt(Lx*Ly)

# find fluid points
z = [(k-.5)*Δz for i = 1:nx, j = 1:ny, k = 1:nz]
fluid = z .> h
fluid[:,:,end] .= false

# viscosity/diffusivity map
z = [(k-1)*Δz for i = 1:nx, j = 1:ny, k = 1:nz]
ν = ν0 .+ ν1*exp.(-(z.-h)/d)

# derivatives
νx = ν1*hx/d.*exp.(-(z.-h)/d)
νy = ν1*hy/d.*exp.(-(z.-h)/d)
νz = -ν1/d.*exp.(-(z.-h)/d)

# save 2D data
h5write(@sprintf("abyssal_%dx%d.h5", nx, ny), "h", h)
h5write(@sprintf("abyssal_%dx%d.h5", nx, ny), "hx", hx)
h5write(@sprintf("abyssal_%dx%d.h5", nx, ny), "hy", hy)

# save to file
#h5write(@sprintf("abyssal_%dx%dx%d.h5", nx, ny, nz), "fluid", convert.(UInt8, fluid))
#h5write("abyssal_512x512x512.h5", "ν", ν)
#h5write("abyssal_512x512x512.h5", "νx", νx)
#h5write("abyssal_512x512x512.h5", "νy", νy)
#h5write("abyssal_512x512x512.h5", "νz", νz)

using Random
using FFTW
using HDF5
using Printf
using SpecialFunctions

# number of grid points
nx = 512
ny = 1

# grid spacing
Δx = 40.
Δy = 40.

# domain size
Lx = nx*Δx
Ly = ny*Δy

# topography parameters
k0 = 1.0e-3
l0 = 2.2e-4
hm = 110.
ν = .9

# 1D spectrum (adapted from Nikurashin and Legg, 2011, eq. 1)
H(k, k0, hm, ν) = 2π*ν*hm^2*gamma(ν+1/2)*k0^2ν/(sqrt(π)*gamma(ν+1))*(k^2 + k0^2)^(-(ν+1/2))

# 2D spectrum (adapted from Nikurashin and Legg, 2011, eq. 1)
H(k, l, k0, l0, hm, ν) = 4π*ν*hm^2/(k0*l0)*((k/k0)^2 + (l/l0)^2 + 1)^(-(ν+1))

# wavenumber vector
Δk = 1/Lx
Δl = 1/Ly
k = reshape(2π*(0:Δk:1/2Δx), (nx÷2+1, 1))
l = reshape(fftshift(2π*(-1/2Δy:Δl:1/2Δy-Δl)), (1, ny))

# spectral coefficients (random phase)
Random.seed!(42)
#S = sqrt.(H.(k, l, k0, l0, hm, ν)).*exp.(2π*1im*rand(nx÷2+1, ny))
S = sqrt.(H.(k, k0, hm, ν)).*exp.(2π*1im*rand(nx÷2+1))

# bottom height
#h = irfft(S, nx)*nx*ny/sqrt(Lx*Ly)
h = irfft(S, nx)*nx/sqrt(Lx)
h .-= minimum(h)

# slopes
#hx = irfft(1im*k.*S, nx)*nx*ny/sqrt(Lx*Ly)
#hy = irfft(1im*l.*S, nx)*nx*ny/sqrt(Lx*Ly)
hx = irfft(1im*k.*S, nx)*nx/sqrt(Lx)
hy = irfft(1im*l.*S, nx)*nx/sqrt(Lx)

# save 2D data
h5write(@sprintf("abyssal_%dx%d_%dx%d.h5", nx, ny, Int(round(Δx)), Int(round(Δy))), "h", h)
h5write(@sprintf("abyssal_%dx%d_%dx%d.h5", nx, ny, Int(round(Δx)), Int(round(Δy))), "hx", hx)
h5write(@sprintf("abyssal_%dx%d_%dx%d.h5", nx, ny, Int(round(Δx)), Int(round(Δy))), "hy", hy)

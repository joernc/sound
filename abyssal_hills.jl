using FFTW
using HDF5

# domain size
Lx = 64e3
Ly = 64e3

# number of grid points
nx = 1024
ny = 1024

# parameters
k0 = 1.0e-3
l0 = 2.2e-4
h = 110.
ν = .9

# 1D spectrum (adapted from Nikurashin and Legg, 2011, eq. 1)
H(k, k0, h, ν) = 2π*ν*h^2*gamma(ν+1/2)*k0^2ν/(sqrt(π)*gamma(ν+1))*(k^2 + k0^2)^(-(ν+1/2))

# 2D spectrum (adapted from Nikurashin and Legg, 2011, eq. 1)
H(k, l, k0, l0, h, ν) = 4π*ν*h^2/(k0*l0)*((k/k0)^2 + (l/l0)^2 + 1)^(-(ν+1))

# grid spacing
Δx = Lx/nx
Δy = Ly/ny

# wavenumber vector
Δk = 1/Lx
Δl = 1/Ly
k = reshape(2π*(0:Δk:1/2Δx), (nx÷2+1, 1))
l = reshape(fftshift(2π*(-1/2Δy:Δl:1/2Δy-Δl)), (1, ny))

# spectral coefficients (random phase)
S = sqrt.(H.(k, l, k0, l0, h, ν)).*exp.(2π*1im*rand(nx÷2+1, ny))

# bottom height
b = irfft(S, nx)*nx*ny/sqrt(Lx*Ly)

# save to file
h5write("abyssal_1024x1024.h5", "b", b)

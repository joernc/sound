using FFTW
using SpecialFunctions
using HDF5

# domain width
L = 60e3

# number of grid points
n = 1024

# parameters
k0 = 1.0e-3
l0 = 2.2e-4
h = 110.
ν = .9

# 1D spectrum (adapted from Nikurashin and Legg, 2011, eq. 1)
H(k, k0, h, ν) = 2π*ν*h^2*gamma(ν+1/2)*k0^2ν/(sqrt(π)*gamma(ν+1))*(k^2 + k0^2)^(-(ν+1/2))

# grid spacing
Δx = L/n

# wavenumber vector
Δk = 1/L
k = 2π*(0:Δk:1/2Δx)

# spectral coefficients (random phase)
S = sqrt.(H.(k, k0, h, ν)).*exp.(2π*1im*rand(n÷2+1))

# bottom height
b = irfft(S, n)*n/sqrt(L)

# save to file
h5write("abyssal.h5", "b", b)

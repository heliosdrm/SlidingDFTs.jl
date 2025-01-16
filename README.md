> [!NOTE]
> This package is not maintained. Its functionality has been transferred to the registered package [FourierTools](https://github.com/bionanoimaging/FourierTools.jl)

# SlidingDFTs

A Julia package to compute [Sliding Discrete Fourer Transforms](https://en.wikipedia.org/wiki/Sliding_DFT) recursively, over one-dimensional series of values.

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://heliosdrm.github.io/SlidingDFTs.jl/dev/)
[![Build Status](https://github.com/heliosdrm/SlidingDFTs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/heliosdrm/SlidingDFTs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/heliosdrm/SlidingDFTs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/heliosdrm/SlidingDFTs.jl)

## Installation

This package is not yet registered. Add it with `]add https://github.com/heliosdrm/SlidingDFTs.jl` in the "pkg mode" of the REPL, or in the "standard REPL":

```julia
using Pkg
Pkg.add("https://github.com/heliosdrm/SlidingDFTs.jl")
```

## Basic usage

This package must be used with some implementation of [AbstractFFTs](https://github.com/JuliaMath/AbstractFFTs.jl/) such as [FFTW](https://github.com/JuliaMath/FFTW.jl), [FastTransforms](https://github.com/JuliaApproximation/FastTransforms.jl) or [RustFFT](https://github.com/Taaitaaiger/RustFFT.jl).

The basic Sliding Discrete Fourier Transform (SDFT) of a one-dimensional series of values `x`, using a window of length `n`, is calculated as follows:

**Step 1**: Setup the method for a SDFT of length `n`:

```julia
using SlidingDFTs
using FFTW # or another package implementing AbstractFFTs

method = SDFT(n)
```

See the [API](@ref) for other methods to compute sliding DFTs besides the basic `SDFT`.

**Step 2**: Create an iterator of the SDFT over `x`, with the function `sdft`. This is typically used in a loop:

```julia
for spectrum in sdft(method, x)
    # `spectrum` is a `Vector{Complex(eltype(x))}` of length `n`
end
```

See the [docs](https://heliosdrm.github.io/SlidingDFTs.jl) for further details.

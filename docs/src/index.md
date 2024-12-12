```@meta
CurrentModule = SlidingDFTs
```

# SlidingDFTs

[SlidingDFTs](https://github.com/heliosdrm/SlidingDFTs.jl) is a Julia package to compute
[Sliding Discrete Fourer Transforms](https://en.wikipedia.org/wiki/Sliding_DFT) recursively, over one-dimensional series of values.

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

## Improve performance with "unsafe" iterations

The iterator created by `sdft` produces a new vector of complex values at each iteration. For better performance, it is possible to call this function with the keyword argument `safe` set to `false`:

```julia
iterator = sdft(method, x, safe=false)
```

This iterator will produce only one vector that will be mutated at each iteration. Use this with caution, as any modification of the resulting vector may lead to unexpected results in subsequent iterations.

## Using SlidingDFTs with stateful iterators

By default, this package computes sliding DFTs traversing sequentially the data series `x`, which can be any kind of iterator. In the case of stateful iterators (i.e. those that are modified upon each iteration, like `Base.Channel`s), that computation will "consume" as many items of `x` as the length of the computed DFT in the first iteration, and one additional item in every subsequent iteration.

Apart from that consideration, it is safe to apply sliding DFTs to stateful iterators, since past samples of `x` already used in previous iterations, which are often required for the computations, are temporarily stored in an array — in internal variables that users do not need to deal with.

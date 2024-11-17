```@meta
CurrentModule = SlidingDFTs
```

# SlidingDFTs

[SlidingDFTs](https://github.com/heliosdrm/SlidingDFTs.jl) is a Julia package to compute
[Sliding Discrete Fourer Transforms](https://en.wikipedia.org/wiki/Sliding_DFT) recursively, on one-dimensional series of values.

## Basic usage

This package must be used with some implementation of [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl/) such as [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) or [FastTransforms.jl](https://github.com/JuliaApproximation/FastTransforms.jl).

The basic Sliding Discrete Fourier Transform (SDFT) of a one-dimensional series of values `x`, using a window of length `n`, is calculated as follows:

**Step 1**: Setup the method for a SDFT of length `n`:

```julia
using SlidingDFTs

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

## Sliding DFTs for stateful iterators

The data series `x` over which sliding DFTs are computed can be any kind of iterator, but most methods require reusing past samples of the data, so `sdft` assumes that `x` is a "stateless" iterator (i.e. it is not altered by iterating over it, so there may be simultaneous loops traversing it without interfering with each other).

That is usually the case in Julia, but if it is not, e.g. in the case of `Base.Channel`s, it is possible to use `stateful_sdft` instead of `sdft`. This will traverse the input `x` only once, storing past values internally to reuse them when needed, at the expense of additional memory allocations.

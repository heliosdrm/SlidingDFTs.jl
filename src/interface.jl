using AbstractFFTs

# exports

## Required functions

"""
    SlidingDFTs.windowlength(method)

Return the length of the window used by `method`.
"""
function windowlength end

"""
    updatedft!(dft, x, method, state)

Update the values of a sliding Discrete Fourier Transform (DFT) of a data series,
according to the algorithm of the provided method, for a given state of the sliding DFT.

`dft` is a mutable collection of complex values with length equal to `windowlength(method)`,
containing the value returned by the last iteration of the sliding DFT.

`x` is the data series for which the sliding DFT is computed, at least as long as
`windowlength(method)`.

`method` is the object that defines the method to compute the sliding DFT.

`state` is an object generated automatically at each iteration of an object created by
`sdft(method, x)` or `stateful_sdft(method, x)`, containing the information that is needed
to compute the new values of the sliding DFT, together with `x`.
That information can be extracted from `state` with the following functions:

* `SlidingDFTs.nextdata`(@ref) to get the next value of the data series.
* `SlidingDFTs.previousdata`(@ref) to get a previous value of the data series.
* `SlidingDFTs.previousdft`(@ref) to get the DFTs of previous iterations.
"""
function updatedft! end

## Conditionally required functions

"""
    dataoffsets(method)

Return an integer or a vector of integers with the offsets of data samples
that are needed by the given method to compute a sliding DFT.

If the code of `SlidingDFTs.updatepdf!` that dispatches on the type of `method` uses the function `SlidingDFTs.previousdata`,
this function must return the integers that are used as the third argument (`offset`) of that function.
If that function is not needed (no past samples are used), this one may return `nothing` to reduce memory allocations.
"""
dataoffsets(::Any) = nothing

"""
    dftback(method)

Return an integer or a vector of positive integers with the indices of the previous iterations
that are needed by the given method to compute a sliding DFT.

If the code of `SlidingDFTs.updatepdf!` for the type of `method` uses the function `SlidingDFTs.previousdft`,
this function must return the integers that are used as the third argument (`back`) of that function.

If that function is not needed, this one may return `nothing` to reduce memory allocations.
"""
dftback(::Any) = nothing


## Iterators

struct SDFTIterator{M, T, S}
    method::M
    data::T
    safe::Bool
    statefulness::S
end

issafe(iterator::SDFTIterator) = iterator.safe
statefultrait(iterator::SDFTIterator) = iterator.statefulness

struct Stateless end
struct Stateful end

"""
    sdft(method, x[, safe=true])

Return an iterator to produce a sliding DFT of `x` using the given `method`.
If `safe == true` (the default behavior), this iterator produces a new vector at each iteration.

Set the optional argument `safe=false` to improve performance by reducing allocations,
at the expense of unexpected behavior if the resulting vector is mutated between iterations.

It is assummed that `x` is a stateless iterator. To work with a stateful iterator, use `stateful_sdft`(@ref) instead.
"""
sdft(method, x, safe=true) = SDFTIterator(method, x, safe, Stateless())

"""
    stateful_sdft(method, x[, safe=true])

Return an iterator to produce a sliding DFT of `x`, using the given `method`,
considering that `x` is a stateful iterator.

See `sdft`(@ref) for more details.
"""
stateful_sdft(method, x, safe=true) = SDFTIterator(method, x, safe, Stateful())


## States

"""
    previousdata(state, x[, offset=0])

Return the first value of the fragment of the data series `x` that was used
in the most recent iteration of the sliding DFT represented by `state`,
or at `offset` positions after the beginning of that fragment.

If the most recent iteration computed the DFT of the fragment of `x`
corresponding to the range `i : i+n`, then this function returns
the `i+offset`-th value of `x`.
"""
function previousdata end

"""
    previousdft(state, x[, back=0])

Return the DFT computed in the most recent iteration
of the sliding DFT represented by `state`, or in a number of previous
iterations equal to `back`, using the input data `x`.

If the most recent iteration computed the DFT of the fragment of `x`
corresponding to the range `i : i+n`, then this function returns the DFT of `x`
in the range `i-back : i+n-back`.
"""
function previousdft end

"""
    nextdata(state, x)

Return the next value after the fragment of the data series `x` that was used
in the most recent iteration of the sliding DFT represented by `state`.

If the most recent iteration computed the DFT of the fragment of `x`
corresponding to the range `i : i+n`, then this function returns
the `i+n+1`-th value of `x`.
"""
function nextdata end
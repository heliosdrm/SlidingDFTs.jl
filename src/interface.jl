using AbstractFFTs
import Base: iterate

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


## States

# State for a stateless data series
struct SDFTState{T, DH, LS, PS}
    dft::Vector{Complex{T}}        # dft returned
    dfthistory::DH                 # record of back dfts if needed, otherwise nothing
    lastdatastate::LS              # state for the last position of the data fragment
    previousdatastates::PS         # states for previous positions of the data fragment if needed, otherwise nothing
    iteration::Int                 # number of the iteration
end

# State for a stateful data series
struct SDFTStatefulState{T, DH, LS, D}
    dft::Vector{Complex{T}}        # dft returned
    dfthistory::DH                 # record of back dfts if needed, otherwise nothing
    lastdatastate::LS              # state for the last position of the data fragment
    datafragment::Vector{D}        # full record of the data fragment
    iteration::Int                 # iteration
end

getdft(state) = state.dft
getdfthistory(state) = state.dfthistory

# State data used by `updatedft!`
struct StateData{H, P, N}
    dfthistory::H
    previousdata::P
    nextdatapoint::N
    windowlength::Int
    iteration::Int
end

hasdfthistory(::StateData{Nothing}) = false
hasdfthistory(::StateData) = true
haspreviousdata(::StateData{<:Any, Nothing}) = false
haspreviousdata(::StateData) = true

# Return the updated state of the iterator, or `nothing` if the data series is consumed.
function updatestate(state::SDFTState, method, x)
    nextiter = iterate(x, state.lastdatastate)
    if nextiter === nothing
        return nothing
    end
    dft = getdft(state)
    dfthistory = getdfthistory(state)
    previousstates, updatedstates = getpreviousdata(x, state.previousdatastates)
    nextdatapoint, nextdatastate = nextiter
    n = windowlength(method)
    statedata = StateData(dfthistory, previousstates, nextdatapoint, n, state.iteration)
    updatedft!(dft, x, method, statedata)
    updatedfthistory!(dfthistory, dft, n, state.iteration)
    return SDFTState(dft, dfthistory, nextdatastate, updatedstates, state.iteration + 1)
end

function getpreviousdata(x, states_dict)
    previousdata = Dict{Int, eltype(x)}()
    updatedstates = copy(states_dict)
    for (i, state) in states_dict
        value, state = iterate(x, state...)
        previousdata[i] = value
        updatedstates[i] = state
    end
    return previousdata, updatedstates
end

getpreviousdata(::Any, ::Nothing) = nothing, nothing

function updatestate(state::SDFTStatefulState, method, x)
    nextiter = iterate(x, state.lastdatastate)
    if nextiter === nothing
        return nothing
    end
    dft = getdft(state)
    dfthistory = getdfthistory(state)
    fragment = state.datafragment
    nextdatapoint, nextdatastate = nextiter
    n = windowlength(method)
    statedata = StateData(dfthistory, fragment, nextdatapoint, n, state.iteration)
    updatedft!(dft, x, method, statedata)
    updatedfthistory!(dfthistory, dft, n, state.iteration)
    updatefragment!(fragment, nextdatapoint, n, state.iteration)
    return SDFTStatefulState(dft, dfthistory, nextdatastate, fragment, state.iteration + 1)
end

"""
    previousdata(state[, offset=0])

Return the first value of the fragment of the data series that was used
in the most recent iteration of the sliding DFT represented by `state`,
or at `offset` positions after the beginning of that fragment.

If the most recent iteration computed the DFT of the fragment of a data series
corresponding to the range `i : i+n`, then this function returns its `i+offset`-th value.
"""
function previousdata(state::StateData{<:Any, <:AbstractVector}, offset=0)
    fragment = state.previousdata
    adjustedoffset = rem(offset + state.iteration - 1, length(fragment))
    return fragment[firstindex(fragment) + adjustedoffset]
end

function previousdata(state::StateData{<:Any, <:AbstractDict}, offset=0)
    fragment = state.previousdata
    return fragment[offset]
end

previousdata(::StateData{<:Any, Nothing}, args...) = throw(ErrorException("define `dataoffsets` as something different to `nothing`"))

"""
    previousdft(state[, back=0])

Return the DFT computed in the most recent iteration
of the sliding DFT represented by `state`, or in a number of previous
iterations equal to `back`.

If the most recent iteration computed the DFT of the fragment of a data series
corresponding to the range `i : i+n`, then this function returns its DFT
for the range `i-back : i+n-back`.
"""
function previousdft(state::StateData, back=0)
    !hasdfthistory(state) && throw(ErrorException("define `dftback` as something different to `nothing`"))
    dfthistory = state.dfthistory
    n = state.windowlength
    nh = length(dfthistory) ÷ n
    offset = rem(state.iteration - back - 1, nh)
    rg = (offset * n) .+ (1:n)
    return view(dfthistory, rg)
end


"""
    nextdata(state)

Return the next value after the fragment of the data series that was used
in the most recent iteration of the sliding DFT represented by `state`.

If the most recent iteration computed the DFT of the fragment of a data series
corresponding to the range `i : i+n`, then this function returns its `i+n+1`-th value.

There is no defined behavior if such value does not exist
(i.e. if the end of a finite data series was reached).
"""
nextdata(state::StateData) = state.nextdatapoint

"""
    iterationcount(state)

Return the number of iterations done for the sliding DFT represented by `state`.

If the most recent iteration computed the DFT of the fragment of `x`
corresponding to the range `i : i+n`, then this function returns  the number `i`.
"""
iterationcount(state) = state.iteration

function updatedfthistory!(dfthistory, dft, n, iteration)
    offset = rem((iteration - 1) * n, length(dfthistory))
    dfthistory[(1:n) .+ offset] .= dft
end

function updatedfthistory!(::Nothing, args...) end

function updatefragment!(fragment, nextdatapoint, n, iteration)
    offset = rem(iteration - 1, n)
    fragment[1 + offset] = nextdatapoint
end

function updatefragment!(::Nothing, args...) end

## Iterators

struct SDFTIterator{M, T, S}
    method::M
    data::T
    safe::Bool
    statefulness::S
end

getmethod(iterator::SDFTIterator) = iterator.method
getdata(iterator::SDFTIterator) = iterator.data
issafe(iterator::SDFTIterator) = iterator.safe
statefultrait(iterator::SDFTIterator) = iterator.statefulness

function Base.length(iterator::SDFTIterator)
    method = getmethod(iterator)
    data = getdata(iterator)
    return length(data) - windowlength(method) + 1
end

Base.eltype(::SDFTIterator{M,T,S}) where {M,T,S} = Vector{Complex{eltype(T)}}

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

function iterate(itr::SDFTIterator)
    windowed_data, previousstates, laststate = initialize(itr)
    dft = fft(windowed_data)
    method = getmethod(itr)
    state = firststate(method, dft, windowed_data, previousstates, laststate, statefultrait(itr))
    returned_dft = itr.safe ? copy(dft) : dft
    return returned_dft, state
end

function iterate(itr::SDFTIterator, state)
    method = getmethod(itr)
    x = getdata(itr)
    newstate = updatestate(state, method, x)
    if newstate === nothing
        return nothing
    end
    dft = getdft(newstate)
    returned_dft = itr.safe ? copy(dft) : dft
    return returned_dft, newstate
end

# Helper functions

# Get the window of the first chunk of data and the state of the data iterator at the end
function initialize(itr)
    x = getdata(itr)
    firstiteration = iterate(x)
    if firstiteration === nothing
        throw(ErrorException("insufficient data"))
    end
    datapoint, datastate = firstiteration
    n = windowlength(getmethod(itr))
    windowed_data = fill(datapoint, n)
    offsetlist = dataoffsets(getmethod(itr))
    previousstates = init_statedict(offsetlist, datastate)
    i = 1
    while i < n
        iteration = iterate(x, datastate)
        if iteration === nothing
            throw(ErrorException("insufficent data"))
        end
        datapoint, datastate = iteration
        if i in offsetlist
            update_statedict!(previousstates, i, datastate)
        end
        i += 1
        windowed_data[i] = datapoint
    end
    return windowed_data, previousstates, datastate
end

init_statedict(::Nothing, ::Any) = nothing
init_statedict(::Any, ::S) where S = Dict{Int, Union{S, Tuple{}, Nothing}}(0 => ())

function update_statedict!(::Nothing, ::Any, ::Any) end
function update_statedict!(dict::AbstractDict, i, state)
    dict[i] = state
end

function firststate(method, dft, ::Any, previousstates, laststate, ::Stateless)
    backindices = dftback(method)
    dfthistory = create_dfthistory(dft, backindices)
    return SDFTState(dft, dfthistory, laststate, previousstates, 1)
end

function firststate(method, dft, windowed_data, ::Any, laststate, ::Stateful)
    backindices = dftback(method)
    dfthistory = create_dfthistory(dft, backindices)
    return SDFTStatefulState(dft, dfthistory, laststate, windowed_data, 1)
end

create_dfthistory(::Any, ::Nothing) = nothing
create_dfthistory(dft, n::Integer) = repeat(dft, n+1)
create_dfthistory(dft, indices) = create_dfthistory(dft, maximum(indices))


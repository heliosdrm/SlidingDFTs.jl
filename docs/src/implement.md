```@meta
CurrentModule = SlidingDFTs
```

# Impementing SlidingDFTs

 SlidingDFTs provides an interface to define new methods to compute Sliding Discrete Fourier Transforms (SDFT), which can be implemented in pull requests to this package, but also in independent packages.

## Theoretical basis

A SDFT is generally implemented as a recursive equation, such that if $X_{i}[k]$ is the DFT of $x[j]$ for $j = i, \ldots, i+n$ at frequency $k$, the next iteration is:

$X_{i+1}[k] = f(X_{i}[k], X_{i-1}[k], \ldots, X_{1}[k], x[i], \ldots x[i+n], x[i+n+1])$

Such an equation depends on:

* The values of the DFT in one or more previous iterations, $X_{p}[k]$ for $p = 1, \ldots, i$.
* The values of the data series used in the most recent iteration, $x[j]$ for $j = i, \ldots, i+n$.
* The next value of the data series after the fragment used in the most recent iteration, $x[i+n+1]$.

For instance, the [basic definition of the SDFT](https://www.researchgate.net/publication/3321463_The_sliding_DFT) uses the following formula: 

$X_{i+1}[k] = W[k] \cdot (X_{i}[k] + x[i+n] - x[i]),$

which depends only on the most recent DFT ($X_{i}[k]$), the first data point of the fragment used in that DFT ($x[i]$), and the next data point after that fragment ($x[i+n+1]$), plus a "twiddle" factor $W[k]$ that only depends on the frequency $k$, equal to $exp(j2{\pi}k/n)$.
Other variations of the SDFT may use formulas that depend on previous iterations or other values of the data series in that fragment.

## Implementation in object types

A method to compute an SDFT is defined by three object types:

* One for the method, which contains the fixed parameters that are needed by the algorithm to compute the SDFT.
* Another for the iterator created by the functions `sdft` or `stateful_sdft`, which binds a method with the target data series.
* And another for the state of the iterator, which complements the method with the variable information that changes at each iteration.

The internals of SlidingDFTs take care of the design of the iterator and state types and the creation of their instances. The only thing that has to be defined to create a new kind of SDFT is the `struct` of the method with the fixed parameters, and a few function methods that dispatching on that type. One of those functions is the one that implements the recursive formula, using also the information stored in the state.

## Extract information from the state of SDFT iterators

Before explaining how to define a new type for SDFTs, it is convenient to know the functions that can be used to extract the information that is stored in the state of SDFT iterators. There is a function for each of the three kinds of variables used in the general recursive equation presented above.

```@docs; canonical=false
previousdft
previousdata
nextdata
```

For instance, the values used in the formula of the basic SDFT may be obtained from a `state` object and its corresponding data series `x` as:
* `SlidingDFTs.previousdft(state, x, 0)` for $X_{i}$.
* `SlidingDFTs.previousdata(state, x, 0)` for $x[i]$.
* `SlidingDFTs.nextdata(state, x)` for $x[i+n+1]$.

Notice that the third argument of `previousdft` and `previousdata` might have been ommited in this case, since it is zero by default.

# Definition of new SDFT types

The design of the `struct` representing a new SDFT type is free, but it is required to implement the following functions dispatching on that type:

```@docs; canonical=false
windowlength
updatedft!
```

The formula of the basic SDFT formula could be implemented for a type `MyBasicSDFT` as follows:

```julia
import SlidingDFTs: updatedft!, windowlength, nextdata, previousdata

function udpatedft!(dft, x, method::MyBasicSDFT, state)
    n = windowlength(method)
    for k in eachindex(dft)
        X_i = dft[k]
        x_iplusn = nextdata(state, x)
        x_i = previousdata(state, x)
        Wk = exp(2π*im*k/n)
        dft[k] = Wk * (X_i + x_iplusn - x_i)
    end
end
```

(The type [`SDFT`]@ref implemented in `SlidingDFTs` actually has as a similar, but not identic definition.)

Notice that it was not necessary to use `previousdft` to retrieve the most recent iteration of the DFT, since that is assummed to be already contained in the first argument `dft`.

Depending on the functions that are used in the particular implementation of `updatepdf!` for a given type, the following methods should be defined too:

```@docs; canonical=false
dataoffsets
dftback
```

The fallback methods of those functions return `nothing`, as if neither `previousdata` or `previousdft` had to be used.

The implementation of `updatepdf!` given in the previous example does use `previousdata` - with the default offset value, equal to zero - so the following is required in this case:

```julia
SlidingDFTs.dataoffsets(::MyBasicSDFT) = 0
```

On the other hand there is no need to define `SlidingDFTs.dftback` in this case, since as it has been noted above, the interface of of `updatedft!` assumes that the most recent DFT is already contained in its first argument, soit is not necessary to use the function `previousdft` to get it.

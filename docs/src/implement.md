# Impementing SlidingDFTs

SlidingDFTs provides an interface to define new methods to compute sliding DFTs, which can be implemented in pull requests to this package, but also in independent packages.

A method to compute sliding DFTs is defined by an object type (`struct`), which contains the parameters that are needed to initialize the DFT with the first fragment of data or update it in subsequent iterations, and has an associated algorithm to perform those calculations. That information is complemented by another type of objects, defined internally in `SlidingDFTs`, that represents the state of the sliding DFT at each iteration.

## Requirements for sliding DFT types

In order to create a new type of sliding DFT, only the `struct` for the method has to be defined. The design of that `struct` is free, but it is required to implement the following functions dispatching on that type:

```@docs
windowlength
updatedft!
```

For instance, the [basic definition of the Sliding DFT](https://www.researchgate.net/publication/3321463_The_sliding_DFT) uses the following recursive formula: 

$X_{i+1}[k] = W[k] \cdot (X_{i}[k] + x[i+n] - x[i])$

with $X_{i}[k]$ being the DFT of the previous fragment evaluated at the frequency $k$, $n$ the length of the window, and $W[k]$ the "twiddle factor" for that frequency (equal to $exp(j2{\pi}k/n)$).

That formula could be implemented for a type `MyBasicSDFT` as follows:

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

Depending on the functions that are used in the particular implementation of `updatepdf!` for a given type, the following methods should be defined too:

```@docs
dataoffsets
dftback
```

The fallback methods of those functions return `nothing`, as if neither `previousdata` or `previousdft` had to be used.

The implementation of `updatepdf!` given in the previous example does use `previousdata` - with the default offset value, equal to zero - so the following is required in this case:

```julia
SlidingDFTs.dataoffsets(::MyBasicSDFT) = 0
```

On the other hand, `previousdft` is not used, so there is no need to define `SlidingDFTs.dftback` in this case. It might have been used to retrieve the most recent value of the DFT (also with the default `back` argument equal to zero), which is used in the recursive formula, but the interface of `updatepdf!` assumes that the most recent DFT is already contained in its first argument (`dft`), so it is not necessary to use the function `previousdft` to get it.

## Basic SDFT

"""
    SDFT(n)

Basic method to compute a Sliding Discrete Fourier Transform of window length `n` [1],
through the recursive formula:

\$X_{i+1}[k] = \\exp(j2{\\pi}k/n) \\cdot (X_{i}[k] + x[i+n] - x[i])\$

The transfer function for the `k`-th bin of this method is:

\$ H(z) = \\frac{1 - z^{-n}}{1 - \\exp(j2{\\pi}k/n) z^{-1}} \$

## References

[1] Jacobsen, E. and Lyons, R. "The sliding DFT," *IEEE Signal Processing Magazine*, 20(2), 74-80 (2003),
DOI: [10.1109/MSP.2003.1184347](https://doi.org/10.1109/MSP.2003.1184347)
"""
struct SDFT{T,C}
    n::T
    factor::C
end

function SDFT(n)
    factor = exp(2π*im/n)
    SDFT(n, factor)
end

# Required functions

windowlength(method::SDFT) = method.n

function updatedft!(dft, x, method::SDFT{T,C}, state) where {T,C}
    twiddle = one(C)
    for k in eachindex(dft)
        dft[k] = twiddle * (dft[k] + nextdata(state) - previousdata(state))
        twiddle *= method.factor
    end
end

dataoffsets(::SDFT) = 0

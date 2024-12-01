## Basic SDFT

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

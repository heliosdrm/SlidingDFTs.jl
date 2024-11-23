using SlidingDFTs
using Test
using FFTW

# Piecewise sinusoidal signal

function signal(x)
    if x < 1
        5*cos(4π*x)
    elseif x < 2
        (-2x+7)*cos(2π*(x^2+1))
    else
        3*cos(10π*x)
    end
end

@testset "SlidingDFTs.jl" begin
    y = signal.(range(0, 3, length=61))
    n = 15
    method = SlidingDFTs.SDFT(n)
    dfty = collect(sdft(method, y))
end

using SlidingDFTs
using Test
using RustFFT

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

y = signal.(range(0, 3, length=61))
n = 20
sample_offsets = (0, 20, 40)
dfty_sample = [fft(view(y, (1:n) .+ offset)) for offset in sample_offsets]

# Compare SDFT
@testset "SDFT" begin
    method = SDFT(n)
    dfty = collect(sdft(method, y))
    @testset "stateless" for i in eachindex(sample_offsets)
        @test dfty[1 + sample_offsets[i]] ≈ dfty_sample[i]
    end
    dfty_stateful_1 = collect(stateful_sdft(method, y))
    dfty_stateful_2 = collect(stateful_sdft(method, Iterators.Stateful(y)))
    dfty_stateful_wrong = collect(sdft(method, Iterators.Stateful(y)))
    @testset "stateful" for i in eachindex(sample_offsets)
        @test dfty_stateful_1[1 + sample_offsets[i]] ≈ dfty_sample[i]
        @test dfty_stateful_2[1 + sample_offsets[i]] ≈ dfty_sample[i]
        if sample_offsets[i] != 0
            @test dfty_stateful_wrong[1 + sample_offsets[i]] ≉ dfty_sample[i]
        end
    end
end

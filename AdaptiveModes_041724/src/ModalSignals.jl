#= Custom signal types designed for easily plugging into adaptive modes. =#

module HarmonicSignals
    import CtrlVQE: Signals, ParametricSignals, ParametricSignal, parameters

    mutable struct Harmonic{F} <: ParametricSignal{F,Complex{F}}
        A::F    # AMPLITUDE - REAL PART
        B::F    # AMPLITUDE - IMAG PART
        n::Int  # HARMONIC NUMBER
        T::F    # FUNDAMENTAL PERIOD (aka T such that ν=n/2T)
    end

    # """ Hard-code the fact that T is not meant to be a variational parameter.

    # (Alternatively, we could use Constrained(my_harmonic, :T),
    #     but this seems simpler.)
    # """
    # function ParametricSignals.parameters(::Harmonic)
    #     return [:A, :B]
    # end

    # TODO: Hard-coding didn't work. Why not?!

    function Signals.valueat(signal::Harmonic{F}, t::Real) where {F}
        A = signal.A; B = signal.B; n = signal.n; T = signal.T
        ν = F(2π * n / 2T)
        return Complex(A, B) * sin(ν*t)
    end

    function Signals.partial(i::Int, signal::Harmonic{F}, t::Real) where {F}
        field = parameters(signal)[i]
        A = signal.A; B = signal.B; n = signal.n; T = signal.T
        ν = F(2π * n / 2T)
        sine = sin(ν*t)
        return (field == :A ?   Complex(sine, zero(F))
            :   field == :B ?   Complex(zero(F), sine)
            :                   error("Not Implemented"))
    end

    function Base.string(::Harmonic, names::AbstractVector{String})
        A, B, n, T = names
        return "($A+i$B) sin(2t⋅$n/$T)"
    end
end
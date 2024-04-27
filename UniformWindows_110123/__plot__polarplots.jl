#= The objective of this file is to play around with radial plots,
    to visualize pulses in different ways.

1. Plot Ω=A⋅exp(iϕ) on complex polar plot, where t/T gives color gradient.
2. Plot ϕ as angle on a polar plot, where r=t. Use A/ΩMAX for color gradient.

Note that this requires 4 distinct color gradients.
Easiest solution there is actually to scale α from 0 to 1, I think?
Then we can just use the four colors we're accustomed to using...

We need to make these for a standard windowed plot,
    and also for a nice continuous plot.

Once we like these, transfer them to the plotting library in `src`.

=#

jobdir = "jobs/lih30_resonant_perns20/T17.0_W340"
jobdir = "jobs/lih30_resonant_ΔsMAX3.0/T21.0_W7"

import Plots

# LOAD DATA
load!(jobdir; run=false)
try include("$jobdir/script.jl") catch end
n = _!.work.system.n
T = _!.setup.T
ΩMAX = _!.setup.ΩMAX
t = CtrlVQE.lattice(_!.work.grid)
Ω, α, β = JOB.make_pulses(_!)

# PREP COLOR GRADIENTS
import ColorSchemes: RGB
cgrads = [
    Plots.cgrad([RGB(0.8,0.8,0.8), RGB(0.0,0.0,0.0)]),
    Plots.cgrad([RGB(0.0,0.0,1.0), RGB(0.0,0.0,0.3)]),
    Plots.cgrad([RGB(0.0,1.0,0.0), RGB(0.0,0.3,0.0)]),
    Plots.cgrad([RGB(1.0,0.0,0.0), RGB(0.3,0.0,0.0)]),
    Plots.cgrad([RGB(0.8,0.8,0.0), RGB(0.2,0.2,0.0)]),
    Plots.cgrad([RGB(0.8,0.0,0.8), RGB(0.2,0.0,0.2)]),
    Plots.cgrad([RGB(0.0,0.8,0.8), RGB(0.0,0.2,0.2)]),
][1:n]


# PLOT t(ϕ)
tϕ = Plots.scatter(; proj=:polar, yticks=0.0:3.0:T, ylims=[0,T], legend=false)
for q in 1:n
    Plots.scatter!(tϕ, angle.(Ω[:,q]), t;
        marker_z = abs.(Ω[:,q])./ΩMAX,
        ms=2,
        msw=0,
        color=Plots.cgrad(cgrads[q], rev=true),
        # colorbar=true,
    )
end

# PLOT A(ϕ)
Aϕ = Plots.plot(; proj=:polar, legend=false)
Plots.plot!(Aϕ, [ϕ for ϕ in range(0,2π,1000)], [ΩMAX for _ in 1:1000]; color=:black)
for q in 1:n
    Plots.plot!(Aϕ, angle.(Ω[:,q]), abs.(Ω[:,q]);
        marker_z = t./T,
        ms=5,
        msw=0,
        color=Plots.cgrad(cgrads[q], rev=true),
        # colorbar=true,
    )
end

Plots.gui()
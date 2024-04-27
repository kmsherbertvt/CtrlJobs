#= My optimal pool is supposed to literally be optimal. It's not. What gives? =#

load!("jobs/lih30_brute/T24.0")
unarchive!(JOB.adaptid(1))
_!.run = false
try include("brute.jl") catch end

ϕ, _, _ = JOB.make_gradientsignals()

pO = JOB.makepool_optimal(_!, ϕ)
pN = JOB.makepool_nodes(_!, ϕ)
pB = JOB.makepool_brute(_!, ϕ)

sO = JOB.score_candidates(_!, pO, ϕ)
sN = JOB.score_candidates(_!, pN, ϕ)
sB = JOB.score_candidates(_!, pB, ϕ)

cO = JOB.select_oneperpulse(pO, sO)
cN = JOB.select_oneperpulse(pN, sN)
cB = JOB.select_oneperpulse(pB, sB)


cd(@__DIR__) # changes the directory to the current directory, the default I guess is the HOME
using Pkg; Pkg.activate("."); Pkg.instantiate()

#= Pkg.add("NeuralPDE")
Pkg.add("Flux")
Pkg.add("ModelingToolkit")
Pkg.add("GalacticOptim")
Pkg.add("Optim")
Pkg.add("DiffEqFlux")
Pkg.add("Plots") =#

using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
import ModelingToolkit: Interval, infimum, supremum

@parameters t, x
@variables u(..)
Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dt = Differential(t)

# PDE
C=1
eq = Dtt(u(t,x)) ~ C^2*Dxx(u(t,x))

# Initial & boundary conditions
bcs = [
    u(t,0) ~ 0., # for all t > 0
    u(t,1) ~ 0., # for all t > 0
    u(0,x) ~ x*(1. - x), # for all 0 < x < 1
    Dt(u(0,x)) ~ 0. # for all 0 < x < 1
]

# Space-time domains
domains = [
    t ∈ Interval(0.0, 1.0),
    x ∈ Interval(0.0, 1.0)
]

# Discretization
dx = 0.1

println("antes de chain")

# Neural Network
chain = FastChain(FastDense(2,16,Flux.σ),FastDense(16,16,Flux.σ),FastDense(16,1))
initθ = Float64.(DiffEqFlux.initial_params(chain))
discretization = PhysicsInformedNN(chain, GridTraining(dx); init_params = initθ)

@named pde_system = PDESystem(eq,bcs,domains,[t,x],[u(t,x)])
prob = discretize(pde_system,discretization)

cb = function (p,l)
    println("Current loss is: $l")
    return false
end

println("antes de optimizer")

# optimizer
opt = Optim.BFGS()
res = GalacticOptim.solve(prob,opt; cb = cb, maxiters=1200)
phi = discretization.phi


using Plots

ts,xs = [infimum(d.domain):dx:supremum(d.domain) for d in domains]
analytic_sol_func(t,x) =  sum([(8/(k^3*pi^3)) * sin(k*pi*x)*cos(C*k*pi*t) for k in 1:2:50000])

u_predict = reshape([first(phi([t,x],res.minimizer)) for t in ts for x in xs],(length(ts),length(xs)))
u_real = reshape([analytic_sol_func(t,x) for t in ts for x in xs], (length(ts),length(xs)))

diff_u = abs.(u_predict .- u_real)
p1 = plot(ts, xs, u_real, linetype=:contourf,title = "analytic");
p2 =plot(ts, xs, u_predict, linetype=:contourf,title = "predict");
p3 = plot(ts, xs, diff_u,linetype=:contourf,title = "error");
plot(p1,p2,p3)
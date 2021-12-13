cd(@__DIR__) # changes the directory to the current directory, the default I guess is the HOME
using Pkg; Pkg.activate("."); Pkg.instantiate()

#Pkg.add("OrdinaryDiffEq")
#Pkg.add("ModelingToolkit")
#Pkg.add("DiffEqOperators")
#Pkg.add("DomainSets")
#Pkg.add("Plots")
# neuralpde
#Pkg.add("NeuralPDE")
#Pkg.add("Flux")
#Pkg.add("GalacticOptim")
#Pkg.add("DiffEqFlux")
#Pkg.add("Optim")

using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, DomainSets
# Method of Manufactured Solutions: exact solution
# dot notation, vectorized, see https://julialang.org/blog/2017/01/moredots/
u_exact = (x,t) -> exp.(-t) * cos.(x)

# Parameters, variables, and derivatives
@parameters t x
@variables u(..)
Dt = Differential(t)
Dxx = Differential(x)^2


# 1D PDE and boundary conditions
eq  = Dt(u(t,x)) ~ Dxx(u(t,x))
bcs = [u(0,x) ~ cos(x),
        u(t,0) ~ exp(-t),
        u(t,1) ~ exp(-t) * cos(1)]

# Space and time domains
domains = [t ∈ Interval(0.0,1.0),
           x ∈ Interval(0.0,1.0)]

# PDE system - ModelingToolkit
@named pdesys = PDESystem(eq,bcs,domains,[t,x],[u(t,x)])

using NeuralPDE, Flux, GalacticOptim, Optim, DiffEqFlux
# Neural Network
chain = FastChain(FastDense(2,16, Flux.σ),FastDense(16,16,Flux.σ),FastDense(16,1))
initθ = Float64.(DiffEqFlux.initial_params(chain))
discretization_nn = PhysicsInformedNN(chain, GridTraining(dx); init_params = initθ)

prob_nn = discretize(pdesys, discretization_nn)

cb = function (p,l)
    println("Current loss is: $l")
    return false
end

# optimizer
opt = Optim.BFGS()
res = GalacticOptim.solve(prob_nn,opt; cb = cb, maxiters=1200)
phi = discretization_nn.phi

# Method of lines discretization
dx = 0.1

discretization = MOLFiniteDifference([x=>dx],t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob,Tsit5(),saveat=0.2)

# Plot results and compare with exact solution
x = (0:dx:1)[2:end-1]
t = sol.t

using Plots
plt = plot()

for i in 1:length(t)
    plot!(x,sol.u[i],label="Numerical, t=$(t[i])")
    plot!(x,res.u[i],label="NN, t=$(t[i])")
    scatter!(x, u_exact(x, t[i]),label="Exact, t=$(t[i])")
    
end
display(plt)
savefig("plot.png")
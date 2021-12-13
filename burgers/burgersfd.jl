cd(@__DIR__) # changes the directory to the current directory, the default I guess is the HOME
using Pkg; Pkg.activate("."); Pkg.instantiate()

using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, DomainSets

@parameters t, x
@variables u(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# PDE
eq  = Dt(u(t,x)) + u(t,x)*Dx(u(t,x)) - (0.01/pi)*Dxx(u(t,x)) ~ 0

# Initial and boundary conditions
bcs = [u(0,x) ~ -sin(pi*x),
       u(t,-1) ~ 0.,
       u(t,1) ~ 0.,
       u(t,-1) ~ u(t,1)]

# Space and time domains
domains = [t ∈ Interval(0.0,1.0),
           x ∈ Interval(-1.0,1.0)]
# Discretization
dx = 0.1

# PDE system - ModelingToolkit
@named pdesys = PDESystem(eq,bcs,domains,[t,x],[u(t,x)])

order = 2
discretization = MOLFiniteDifference([x=>dx],t;centered_order=order)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob,Tsit5(),saveat=0.1)

#= # Plot results and compare with exact solution
x = (0:dx:1)[2:end-1]
t = sol.t

using Plots
plt = plot()

for i in 1:length(t)
    plot!(x,sol.u[i],label="Numerical, t=$(t[i])")
end
display(plt)
savefig("plot.png") =#
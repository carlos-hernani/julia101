# Automated Finite Difference Solution to the Heat Equation
## by [SciML-DiffEqOperators](https://github.com/SciML/DiffEqOperators.jl)

## __Heat Equation__

$$
\Huge
\partial_t u(t, x) = \alpha\Delta u
\\
\alpha=1
$$

___Boundary Conditions___
$$
\large
u(0,x) = cos(x) \\ u(t,0)= exp(-t) \\ u(t, \pi)=-exp(-t)
$$

___Space-time Domain___
$$
\large
t \in [0,1] \\
x \in [0,\pi]
$$
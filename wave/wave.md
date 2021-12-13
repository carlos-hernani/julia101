# 1D Wave Equation with Dirichlet boundary conditions
## by [SciML-NeuralPDE](https://neuralpde.sciml.ai/dev/pinn/wave/#D-Wave-Equation-with-Dirichlet-boundary-conditions)

## __Wave Equation__

$$
\Huge
\partial_{t}^2 u(t, x) = c^2\Delta u(t,x)

$$

___Boundary Conditions___
$$
\large
u(0,x) = 0 \\ \partial_t u(0,x) = x(1-x) \\ u(t,0)= u(t, 1)= 0
$$

___Space-time Domain___
$$
\large
t \in [0,+\inf) \\
x \in [0,1]
$$
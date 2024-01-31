# Nonlinear difference equations
This Rstudio function allows you to study the stability of the equilibrium points of non-linear, first-order and autonomous difference equations: $x_{t+1} = f(x_t)$, where $f(x_t)$ is a non-linear function. At the moment, this Rstudio function only allows $f(x_t)$ to be a polynomial, as in these cases:: $x_{t+1} = \frac{x^3 + x^2 + 2x - 2}{x^2 + 1}$, $x_{t+1} = x^3 + x^2 + 2x - 2$.

Particularly, the function, which can be found in [code.R](/code.R), performs the following tasks:
- calculate the equilibrium points and identify whether they are stable (locally asymptotically stable) or unstable
- graph the phase diagram and represent the trajectories, based on initial conditions

An example of the expected output of this Rstudio function can be found in [documento.pdf](/documento.pdf) (in Spanish).

# Integro-Differential Equation Solver using Euler-Cromer Method

This project numerically solves a second-order integro-differential equation using the Euler-Cromer method. The system models dynamic behavior with memory effects, where the acceleration depends on both the current state and a convolution integral involving incomplete gamma functions.

## 🧠 Mathematical Formulation

We consider the following integro-differential equation:

$$
\frac{d^2U}{dt^2} + \lambda^2 U(t) - \epsilon \lambda^2 \int_0^t k(t - S) U(S) \, dS = f(t)
$$

Where:
- \( U(t) \) is the unknown function (displacement)
- \( \lambda \) is a stiffness-like parameter
- \( \epsilon \) controls the strength of the memory term
- \( k(t - S) \) is a kernel function involving incomplete gamma functions
- \( f(t) = \sigma_0 \) is a constant forcing term

The kernel \( k(t - S) \) is defined as:

$$
k(t - S) = \frac{1}{\beta^\alpha} \left[ \Gamma(\alpha) \left( \text{gammainc}(\alpha, \beta(t - S)) \right)' \right]
$$

Additionally, a velocity-dependent term is included via:

$$
\int_0^t h(t - S) \frac{dU}{dS} \, dS
$$

Where \( h(t - S) \) also involves incomplete gamma functions with shifted parameters.

## 🔢 Numerical Method: Euler-Cromer

The Euler-Cromer method is a semi-implicit scheme ideal for second-order systems:

- Update velocity:  
  $$ V_{i} = V_{i-1} + \left( \text{RHS of equation} \right) \cdot dt $$
- Update position:  
  $ U_{i} = U_{i-1} + V_{i} \cdot dt $

The integral terms are approximated using discrete summation over previous time steps.

## ⚙️ Parameters

```python
epsilon = 0.1       # Memory strength
alpha = 0.074       # Shape parameter for gamma function
beta = 0.05         # Scale parameter
lambda_ = 100       # Stiffness-like coefficient
sigma0 = 1.0        # Constant forcing
phi0 = 0.0          # Initial displacement
phi1 = 1.0          # Initial velocity

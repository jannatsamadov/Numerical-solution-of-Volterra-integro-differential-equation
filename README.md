# Numerical-solution-of-Volterra-integro-differential-equation
# Integro-Differential Equation Solver using Euler-Cromer Method

This project numerically solves a second-order integro-differential equation using the Euler-Cromer method. The equation models a dynamic system with memory effects, where the evolution of the function \( U(t) \) depends on both its current state and a convolution-like integral involving incomplete gamma functions.

## 🧠 Mathematical Model

The equation is of the form:

$$
\frac{d^2U}{dt^2} + \lambda^2 U(t) - \epsilon \lambda^2 \int_0^t k(t - S) U(S) \, dS = f(t)
$$

\[
\frac{d^2U}{dt^2} = \sigma_0 - \lambda^2 U(t) + \lambda^2 \epsilon \int_0^t \left[ A(t - \tau) U(\tau) + B(t - \tau) \frac{dU}{d\tau} \right] d\tau
\]



Where:
- \( A(t - \tau) \) and \( B(t - \tau) \) involve incomplete gamma functions
- \( \epsilon, \alpha, \beta, \lambda, \sigma_0 \) are system parameters

The integral terms capture hereditary effects, making this a **Volterra-type integro-differential equation**.

## 🔢 Numerical Method

The solution is computed using the **Euler-Cromer method**, a variant of Euler's method suitable for second-order systems:

- Velocity \( V(t) \) is updated using the derivative
- Position \( U(t) \) is updated using the new velocity

Time is discretized with a fixed step size \( dt \), and the integral is approximated using a Riemann sum.

## 📐 Parameters

```python
epsilon = 0.1
alpha = 0.074
beta = 0.05
lambda_ = 100
sigma0 = 1.0
phi0 = 0.0
phi1 = 1.0

# Integro-Differential Equation Solver using Euler-Cromer Method

This project numerically solves a second-order integro-differential equation using the Euler-Cromer method. The equation models a dynamic system with memory effects, where the evolution of the function \( U(t) \) depends on both its current state and a convolution-like integral involving incomplete gamma functions.

## 🧠 Mathematical Model

## Problem Overview

### Governing Equation
The system is described by:

$$
\frac{d^2U}{dt^2} + \lambda^2 U(t) - \epsilon \lambda^2 \int_0^t k(t - S) U(S) \, dS = f(t)
$$

Where:
- `U(t)` is the unknown function
- `k(t) = ε t^(α-1) e^(-β t)` is the memory kernel
- `f(t) = σ₀·H(t)` is the forcing function (step function)
- Initial conditions: `U(0) = 0`, `U'(0) = 1`

### Physical Significance
This equation models systems with:
- Harmonic oscillator term (`λ² U(t)`)
- Memory-dependent effects (integral term)
- Step-function external forcing
- Applications in viscoelastic materials, nuclear reactor dynamics, and systems with hereditary effects

## Numerical Approach

### Solution Method
1. **System Transformation**:
   - Convert to system of first-order ODEs:
     ```
     dU/dt = V
     dV/dt = f(t) - λ² U(t) + ε λ² ∫₀ᵗ k(t-S) U(S) dS
     ```

2. **Time Discretization**:
   - Uniform grid: `tᵢ = i·Δt` (i=0,1,...,N)
   - Time step: `Δt = 0.001`

3. **Integral Approximation**:
   - Break integral into subintervals `[t_j, t_{j+1}]`
   - Use linear approximation for `U(S)` on each subinterval
   - Analytical evaluation using incomplete gamma functions

4. **Time Integration**:
   - Euler-Cromer method (semi-implicit Euler)
   - Update equations:
     ```
     Vᵢ = Vᵢ₋₁ + [σ₀ - λ² Uᵢ₋₁ + Iᵢ] Δt
     Uᵢ = Uᵢ₋₁ + Vᵢ Δt
     ```

### Key Features
- Handles weakly singular kernel `t^(α-1)`
- Efficient computation of memory integral
- Stable time integration for oscillatory systems
- Accurate treatment of initial conditions

## Code Implementation

### Requirements
- Python 3.6+
- NumPy
- SciPy
- Matplotlib

### Parameters
| Parameter | Value | Description |
|----------|-------|-------------|
| `ε`      | 0.1   | Nonlinearity strength |
| `α`      | 0.074 | Kernel exponent |
| `β`      | 0.05  | Exponential decay rate |
| `λ`      | 100   | Natural frequency |
| `σ₀`     | 1.0   | Forcing amplitude |
| `t_max`  | 1.0   | Simulation time |
| `Δt`     | 0.001 | Time step |

### Usage
1. Install dependencies:
   ```bash
   pip install numpy scipy matplotlib
   ```

2. Run the script:
   ```bash
   python integro_diff_solver.py
   ```

3. Output:
   - Plot of `U(t)` vs `t`
   - Solution displayed in matplotlib window

## Results
The solution exhibits:
- Initial transient response
- Oscillatory behavior at frequency `λ/(2π)`
- Gradual decay due to memory effects
- Asymptotic approach to steady state

![Solution Plot](solution_plot.png)

## References
1. Podlubny, I. (1998). Fractional Differential Equations
2. Mainardi, F. (2010). Fractional Calculus and Waves in Linear Viscoelasticity
3. Linz, P. (1985). Analytical and Numerical Methods for Volterra Equations

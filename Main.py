
#This code provides U(t) vs (t) plot

import numpy as np

from scipy.special import gamma, gammainc

import matplotlib.pyplot as plt

# Parameters

epsilon = 0.1

alpha = 0.074

beta = 0.05

lambda_ = 100

sigma0 = 1.0

phi0 = 0.0

phi1 = 1.0

Gamma_alpha = gamma(alpha)

Gamma_alpha_plus1 = gamma(alpha + 1)

  

# Time parameters

t_max = 1.0  # Simulation time

dt = 0.001  # Time step

num_steps = int(t_max / dt) + 1

t = np.linspace(0, t_max, num_steps)

  

# Initialize arrays

U = np.zeros(num_steps)

V = np.zeros(num_steps)

U[0] = phi0

V[0] = phi1

  

# Precompute beta powers

beta_alpha = beta**alpha

beta_alpha_p1 = beta**(alpha + 1)

  

# Compute the solution using Euler-Cromer method

for i in range(1, num_steps):

    integral = 0.0

    for j in range(i):

        x_lower = max((i - j - 1) * dt, 0)  # Avoid negative x_lower

        x_upper = (i - j) * dt

        a_lower = beta * x_lower

        a_upper = beta * x_upper

        # Compute incomplete gamma functions

        Gamma_lower = Gamma_alpha * (1 - gammainc(alpha, a_lower))

        Gamma_upper = Gamma_alpha * (1 - gammainc(alpha, a_upper))

        A_j = epsilon * U[j] * (Gamma_lower - Gamma_upper) / beta_alpha

        c = (i - j) * dt

        term1 = c * (Gamma_lower - Gamma_upper) / beta_alpha

        Gamma_plus_lower = Gamma_alpha_plus1 * (1 - gammainc(alpha + 1, a_lower))

        Gamma_plus_upper = Gamma_alpha_plus1 * (1 - gammainc(alpha + 1, a_upper))

        term2 = (Gamma_plus_lower - Gamma_plus_upper) / beta_alpha_p1

        B_j = epsilon * (U[j+1] - U[j]) / dt * (term1 - term2)

        integral += A_j + B_j

    integral *= lambda_**2  # Multiply by λ²

    V_prime = sigma0 - (lambda_**2) * U[i-1] + integral

    # Update velocity and position using Euler-Cromer

    V[i] = V[i-1] + V_prime * dt

    U[i] = U[i-1] + V[i] * dt

  

# Plotting

plt.figure(figsize=(10, 6))

plt.plot(t, U, label='U(t)')

plt.xlabel('t')

plt.ylabel('U(t)')

plt.title('Solution of the Integro-Differential Equation')

plt.grid(True)

plt.legend()

plt.show()

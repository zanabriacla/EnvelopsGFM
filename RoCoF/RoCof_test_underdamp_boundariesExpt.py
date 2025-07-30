import numpy as np
import matplotlib.pyplot as plt

# Parameters

A = -0.145
B = 1.5458233232247594e-07
C = 0.14498454176676773
D = 0.519347758204538
xi = 0.3265440677878888  # Damping ratio (epsilon)
omega_d = 5.168750951105869  # Damped frequency (wd)
omega_n = 5.468524655221178  # Natural frequency (wn)
T_pll = 0.01  # Assumed (adjust if needed)

# Time vector
t = np.linspace(0, 5, 1000)  # Time from 0 to 5 seconds

# Compute -ΔP(t)
exp_decay_pll = (B / T_pll) * np.exp(-t / T_pll)
oscillatory_part = C * np.exp(-xi * omega_n * t) * np.cos(omega_d * t)
damping_correction = (C * (D / C - xi * omega_n) / omega_d) * np.exp(-xi * omega_n * t) * np.sin(omega_d * t)
deltaP = A + exp_decay_pll + oscillatory_part + damping_correction
deltaP = deltaP*-1

# Compute the envelopes
exp_decay_pll = (B / T_pll) * np.exp(-t / T_pll)
R = np.sqrt(C**2 + ((D - C * xi * omega_n) / omega_d)**2)  # Amplitude factor
upper_env = A + exp_decay_pll + R * np.exp(-xi * omega_n * t)
lower_env = A + exp_decay_pll - R * np.exp(-xi * omega_n * t)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(t, deltaP, label=r'$-\Delta P(t)$', color='blue')
plt.plot(t, -upper_env, '--', label='Upper Envelope', color='red')
plt.plot(t, -lower_env, '--', label='Lower Envelope', color='green')
plt.xlabel('Time (s)')
plt.ylabel(r'$-\Delta P(t)$')
plt.title('Function $-\Delta P(t)$ with Upper and Lower Envelopes')
plt.legend()
plt.grid(True)
plt.show()
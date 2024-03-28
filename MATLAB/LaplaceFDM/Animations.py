import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.animation import PillowWriter


metadata = dict(title='Movie', artist='Cauchy')
writer = PillowWriter(fps=30, metadata=metadata)

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

# Define constants
bSteps = 79  # Amount of steps on the boundary
tSteps = 1000  # Amount of time steps
dr = (2 * np.pi) / (bSteps + 1)  # Delta r
dtheta = (2 * np.pi) / (bSteps + 1)  # Delta theta
k = 1e1  # Diffusivity Constant
dt = (k / 2) * ((dr ** 2) + (dtheta ** 2)) * 0.001  # Delta time

# Initialize tensor u(x,y,t)
u = np.zeros((bSteps + 1, bSteps + 1, tSteps))

# Initial Condition function
def IC(r, theta):
    return np.cos(r)


# Boundary Condition function
def BC(t, theta):
    return 1


# Initialize initial conditions
for i in range(bSteps + 1):
    for j in range(bSteps + 1):
        u[i, j, 0] = IC(i * dr, j * dtheta)


# Main loop
for t in range(1, tSteps):
    for i in range(bSteps + 1)[::-1]:
        for j in range(bSteps + 1):
            ri = dr * (i+1)

            psi = ((-2 * k * dt) * (((ri ** 2) * (dtheta ** 2)) + (dr ** 2))) / ((ri ** 2) * (dr ** 2) * (dtheta ** 2))
            alpha = ((2 * ri + dr) * k * dt) / (2 * ri * (dr ** 2))
            beta = ((2 * ri - dr) * k * dt) / (2 * ri * (dr ** 2))
            phi = (k * dt) / ((ri ** 2) * (dtheta ** 2))

            if i == 0:
                u[i, j, t] = (((-4 * dt * k) / (dr ** 2)) + 1) * u[0, 0, t - 1] + (dt * k / (dr ** 2)) * (
                        u[1, 0, t - 1] + u[1, (bSteps // 4), t - 1] + u[1, (bSteps // 2), t - 1] +
                        u[1, bSteps, t - 1])
            elif i == bSteps and j == 0:
                u[i, j, t] = ((psi + 1) * u[i, j, t - 1]) + (alpha * BC(dt * t, dtheta * j)) + (
                            beta * u[i - 1, j, t - 1]) + (
                                     phi * u[i, j + 1, t - 1]) + (phi * u[i, bSteps, t - 1])
            elif i == bSteps and j == bSteps:
                u[i, j, t] = ((psi + 1) * u[i, j, t - 1]) + (alpha * BC(dt * t, dtheta * j)) + (
                            beta * u[i - 1, j, t - 1]) + (
                                     phi * u[i, 0, t - 1]) + (phi * u[i, j - 1, t - 1])
            elif i == bSteps:
                u[i, j, t] = ((psi + 1) * u[i, j, t - 1]) + (alpha * BC(dt * t, dtheta * j)) + (
                            beta * u[i - 1, j, t - 1]) + (
                                     phi * u[i, j + 1, t - 1]) + (phi * u[i, j - 1, t - 1])
            elif j == 0:
                u[i, j, t] = ((psi + 1) * u[i, j, t - 1]) + (alpha * u[i + 1, j, t - 1]) + (
                            beta * u[i - 1, j, t - 1]) + (
                                     phi * u[i, j + 1, t - 1]) + (phi * u[i, bSteps, t - 1])
            elif j == bSteps:
                u[i, j, t] = ((psi + 1) * u[i, j, t - 1]) + (alpha * u[i + 1, j, t - 1]) + (
                            beta * u[i - 1, j, t - 1]) + (
                                     phi * u[i, 0, t - 1]) + (phi * u[i, j - 1, t - 1])
            else:
                u[i, j, t] = ((psi + 1) * u[i, j, t - 1]) + (alpha * u[i + 1, j, t - 1]) + (
                            beta * u[i - 1, j, t - 1]) + (
                                     phi * u[i, j + 1, t - 1]) + (phi * u[i, j - 1, t - 1])

# Convert polar coordinates to Cartesian coordinates for all slices
theta = np.linspace(0, 2 * np.pi, bSteps + 1)
r = np.linspace(0, 1, bSteps + 1)
Theta, R = np.meshgrid(theta, r)
X = R * np.cos(Theta)
Y = R * np.sin(Theta)

with writer.saving(fig, "aniTest.gif", 200):
    for t in range(0, tSteps, 5):
        z = u[:, :, t]
        ax.set_zlim(-1, 1)
        ax.plot_surface(X, Y, z, cmap=cm.viridis)
        writer.grab_frame()
        plt.cla()

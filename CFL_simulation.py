import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
dt_range = np.linspace(0, 0.01, 10)
h_range = np.linspace(0, 0.003, 10)


print(dt_range, h_range)

# create mesh using dt and h
dt, h = np.meshgrid(dt_range, h_range)

# plot scatter of dt and h
ax.scatter(dt, h, s=10, c='lightgreen', marker='o')

# # plot the boundary line
D = 0.0002
f = lambda dt: (D*dt/0.25)**0.5
ax.plot(np.linspace(0, 0.01, 100), f(np.linspace(0, 0.01, 100)), c='black', label='D*dt/h^2 = 0.25')

ax.set_xlabel('dt')
ax.set_ylabel('h')
ax.set_xlim(0, 0.011)
ax.set_ylim(0, 0.0035)
ax.legend()
ax.set_title('dt vs h')
plt.show()
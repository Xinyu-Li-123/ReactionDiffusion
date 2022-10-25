import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage import convolve

"""
Bimolecular reaction with diffusion on a torus
A + B <-> AB    on a torus T \in R^2

dA  / dt = -k1*A*B + k2*AB + D_A*Laplacian(A) + f_A(x, t)
dB  / dt = -k1*A*B + k2*AB + D_B*Laplacian(B) + f_B(x, t)
dAB / dt = k1*A*B - k2*AB + D_AB*Laplacian(AB) + f_AB(x, t)
where
    A is the concentration of A [m^-3]
    B is the concentration of B [m^-3]
    AB is the concentration of AB [m^-3]
    k1 is an coefficient [m^3*s^-1]
    k2 is an coefficient [s^-1]
    D_A is the diagonal matrix of diffusion coefficient of A
"""
""" Initialization """

## Initialize
species = np.zeros(3, dtype=float)   # [A, B, AB]
index2name = ["A", "B", "AB"]
# reaction coefficients
k1 = 1.
k2 = 0.
# diffusion coefficients of A, B, AB at x and y direction, a value in [0, 1]
D_A = 0.23
D_B = 0.23
D_AB = 0.23
# feed / kill species A / B / AB at certain time
#TODO: implement this feature

# Define a 2d rectangle T. We will wrap it into a torus in the algorithm
# T is a 3-dim array. Fix time t, T[0, x1, x2] is the concentration of species[0] 
#   at location (x1, x2)
T_size = (100, 100)
h = 0.2       # gap of the grid, the smaller h, the more refined the grid is.
print(f"A grid of size {T_size[0]*h:.2f}m*{T_size[1]*h:.2f}m cutted into {T_size[0]}*{T_size[1]} blocks, gap h = {h:.2f}m")
T = np.zeros((len(species), *T_size), dtype=float)
# T[0, :70, :70] = 5*np.random.random((70,70))
# T[1, 30:, 30:] = 5*np.random.random((70,70))
T[0, 20:50, 20:50] = 5*np.ones((30,30))
T[1, 40:70, 40:70] = 5*np.ones((30,30))
# T[0,100:150,100:150] = 50*np.ones((50,50))
# T[1,200:250,200:250] = 70*np.ones((50,50))
# T[2,3:13,3:13] = 0*np.ones((10,10))


## Animation setting
save_animation = True
colormaps = ["Reds", "Blues", "Purples"]      # colormap for each species

# cmap_min = -np.max(T)     # use this setting to visually check for negative concentration (numerical instability)
# cmap_min = np.min(T)
# cmap_max = np.max(T)
cmap_min = 0
cmap_max = 5


## Numerical solution setting
dt = 0.01   # second
interval = 10   # interval between frames, milisecond
playback_speed = dt * 1000 / interval
total_frame = 40000
print(f"dt = {dt}s, Playback speed: x{playback_speed:.1f}, Total time: {dt*total_frame:.2f}s")
print(f"CFD number is A: {D_A*dt/h**2:.2f}, B={D_B*dt/h**2:.2f}, AB={D_AB*dt/h**2:.2f}")


Astars = np.zeros(total_frame)
Bstars = np.zeros(total_frame)
Astars[0] = np.sum(T[[0,2],:,:])
Bstars[0] = np.sum(T[[1,2],:,:])

# laplacian_matrix = np.array([     
#     [1, 1, 1],
#     [1, -8, 1],
#     [1, 1, 1],
# ])
laplacian_matrix = 1 / h**2 * np.array([     
    [0, 1, 0],
    [1, -4, 1],
    [0, 1, 0],
])


def update(frame_index):
    global T
    # Dummy update
    # T = np.random.random((len(species), *T_size))
    

    # change caused by A + B -> AB
    delta_reaction = -k1*(T[0,:,:]*T[1,:,:]) + k2*T[2,:,:]
    delta_diffusion = np.zeros((len(species), *T_size))
    for index in range(len(species)):
        delta_diffusion[index,:,:] = convolve(T[index,:,:], laplacian_matrix, mode="wrap")

    T[0,:,:] += ( delta_reaction + D_A  * delta_diffusion[0, :, :]) * dt
    T[1,:,:] += ( delta_reaction + D_B  * delta_diffusion[1, :, :]) * dt
    T[2,:,:] += (-delta_reaction + D_AB * delta_diffusion[2, :, :]) * dt
    

    # for index in range(len(species)):
    #     # axes[index].set_title(f"Concentration of Species {index2name[index]} at {dt*frame_index:.2f}s, avg={np.average(T[index, :, :]):.2f}")
    #     axes[index].set_title(f"[{index2name[index]}] at {dt*frame_index:.2f}s, avg={np.average(T[index, :, :]):.2e}")
    #     mats[index].set_data(T[index, :, :])

    # axes[-1].set_title("[A*]={Astar_cur:.2e}/{Astar_init:.2e}, [B*]={Bstar_cur:.2e}/{Bstar_init:.2e}".format(
    #     Astar_cur = Astars[frame_index], Astar_init = Astars[0],
    #     Bstar_cur = Bstars[frame_index], Bstar_init = Bstars[0]
    # ))
    # Astars[frame_index] = np.sum(T[[0,2],:,:])
    # Bstars[frame_index] = np.sum(T[[1,2],:,:])

    # axes[-1].set_title("[A*]={Astar_cur:.2e}/{Astar_init:.2e}, [B*]={Bstar_cur:.2e}/{Bstar_init:.2e}".format(
    #     Astar_cur = Astars[frame_index], Astar_init = Astars[0],
    #     Bstar_cur = Bstars[frame_index], Bstar_init = Bstars[0]
    # ))
    # A_line.set_xdata(np.linspace(0, dt*frame_index, frame_index))
    # A_line.set_ydata(Astars[:frame_index])
    # B_line.set_xdata(np.linspace(0, dt*frame_index, frame_index))
    # B_line.set_ydata(Bstars[:frame_index])


# fig, axes = plt.subplots(1, 4, figsize=(10, 10))
# mats = []

# ## initialize the plot
# # initialize animation plot
# for index in range(len(species)):
#     axes[index].set_xlim(0, T_size[0])
#     axes[index].set_ylim(0, T_size[1])
#     axes[index].set_title(f"Concentration of Species {index2name[index]} at {0:.2f}s")
#     axes[index].set_aspect("equal")
#     axes[index].axis("off")
#     mat = axes[index].matshow(
#         T[index, :, :],
#         cmap = colormaps[index],
#         vmin = cmap_min,
#         vmax = cmap_max,
#     )
#     mats.append(mat)

# # initialize analysis plot
# axes[-1].set_xlim(0, total_frame * dt * 1.1)
# axes[-1].set_ylim(
#     0.5*np.min([Astars[0], Bstars[0]]), 
#     1.2*np.max([Astars[0], Bstars[0]]))
# axes[-1].set_aspect("auto")
# axes[-1].set_title("[A*]={Astar_cur:.2e}/{Astar_init:.2e}, [B*]={Bstar_cur:.2e}/{Bstar_init:.2e}".format(
#     Astar_cur = Astars[0], Astar_init = Astars[0],
#     Bstar_cur = Bstars[0], Bstar_init = Bstars[0]
# ))
# A_line, = axes[-1].plot(range(1), Astars[0], label="A*")
# B_line, = axes[-1].plot(range(1), Bstars[0], label="B*")
# axes[-1].legend()

# ## run animation
# ani = animation.FuncAnimation(
#     fig, 
#     update, 
#     frames=total_frame, 
#     interval=interval, 
#     repeat = False,
#     save_count=50)

# if save_animation:
#     # save the animation
#     ani.save('BimolecularReactionDiffusion.mp4', writer='ffmpeg', fps=30)
# else:
#     plt.show()


def simulate_no_plot():
    print("Simulating without animation...")    
    fig, axes = plt.subplots(1, 4, figsize=(10, 10))
    axes[0].set_title("CDF of A")
    # axes[0].set_xlim(0, np.max(T[0,:,:])*1.1)
    # axes[0].set_ylim(0, 1)
    # axes[0].set_aspect("auto")

    axes[1].set_title("CDF of B")
    # axes[1].set_xlim(0, np.max(T[1,:,:])*1.1)
    # axes[1].set_ylim(0, 1)
    # axes[1].set_aspect("auto")

    for frame_index in range(1, total_frame+1):
        if frame_index % 2000 == 0:
            print(f"frame {frame_index} / {total_frame}")

        update(frame_index)

        # plot CDF of A and B on axes[0] and axes[1] with label
        if frame_index > 19999 and frame_index % 5000 == 0:
            for species_index in range(2):
                # The calculation of CDF is based on the following link:
                #    https://www.geeksforgeeks.org/how-to-calculate-and-plot-a-cumulative-distribution-function-with-matplotlib-in-python/
                count, bins_count = np.histogram(T[0,:,:].flatten(), bins=100)
                
                # finding the PDF of the histogram using count values
                pdf = count / sum(count)
                
                # using numpy np.cumsum to calculate the CDF
                # We can also find using the PDF values by looping and adding
                cdf = np.cumsum(pdf)
                
                # plotting PDF and CDF
                axes[species_index].plot(bins_count[1:], cdf, label=f"CDF at {frame_index*dt:.2f}s", linewidth=2)
    


    axes[0].legend()
    axes[1].legend()

    # axes[0].hist(T[0,:,:].flatten(), bins=100)
    # axes[0].set_title("Histogram of [A]")
    # axes[1].hist(T[1,:,:].flatten(), bins=100)
    # axes[1].set_title("Histogram of [B]")

    axes[2].matshow(T[0,:,:], cmap="Reds")
    axes[2].set_title("Concentration of [A]")
    axes[3].matshow(T[1,:,:], cmap="Blues")
    axes[3].set_title("Concentration of [B]")
    plt.show()

simulate_no_plot()
print(f"In the end, average of [A] = {np.mean(T[0,:,:])}")
print(f"In the end, average of [B] = {np.mean(T[1,:,:])}")
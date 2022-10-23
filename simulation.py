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
k1 = 0.5
k2 = 0.5
# diffusion coefficients of A, B, AB at x and y direction, a value in [0, 1]
D_A = 1.
D_B = 1.
D_AB = 1.
# feed / kill species A / B / AB at certain time
#TODO: implement this feature

# Define a 2d rectangle T. We will wrap it into a torus in the algorithm
# T is a 3-dim array. Fix time t, T[0, x1, x2] is the concentration of species[0] 
#   at location (x1, x2)
T_size = (40, 40)
T = np.zeros((len(species), *T_size), dtype=float)
# T[0, :20, :20] = 100*np.random.random((20,20))
# T[1, :20, :20] = 100*np.random.random((20,20))
T[0,15:25,15:25] = 50*np.ones((10,10))
T[1,5:15,5:15] = 70*np.ones((10,10))
T[2,3:13,3:13] = 0*np.ones((10,10))


## Animation setting
save_animation = False
colormaps = ["Reds", "Blues", "Purples"]      # colormap for each species
# cmap_min = -np.max(T)
cmap_min = np.min(T)
cmap_max = np.max(T)


## Numerical solution setting
dt = 0.05   # second
interval = 10   # interval between frames, milisecond
playback_speed = dt * 1000 / interval
total_frame = 5000
print(f"dt = {dt}s, Playback speed: x{playback_speed:.1f}, Total time: {dt*total_frame:.2f}s")

Astars = np.zeros(total_frame)
Bstars = np.zeros(total_frame)
Astars[0] = np.sum(T[[0,2],:,:])
Bstars[0] = np.sum(T[[1,2],:,:])

# laplacian_matrix = np.array([     # numerically unstable?
#     [1, 1, 1],
#     [1, -8, 1],
#     [1, 1, 1],
# ])

# laplacian_matrix = np.array([
#     [0.05, 0.2, 0.05],
#     [0.2,  -1 , 0.2 ],
#     [0.05, 0.2, 0.05],
# ])

laplacian_matrix = 0.05 * np.array([     # larger laplacian
    [1, 1, 1  , 1, 1],
    [1, 1, 1  , 1, 1],
    [1, 1, -24, 1, 1],
    [1, 1, 1  , 1, 1],
    [1, 1, 1  , 1, 1],
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
    

    for index in range(len(species)):
        # axes[index].set_title(f"Concentration of Species {index2name[index]} at {dt*frame_index:.2f}s, avg={np.average(T[index, :, :]):.2f}")
        axes[index].set_title(f"[{index2name[index]}] at {dt*frame_index:.2f}s, avg={np.average(T[index, :, :]):.2e}")
        mats[index].set_data(T[index, :, :])

    # axes[-1].set_title("[A*]={Astar_cur:.2e}/{Astar_init:.2e}, [B*]={Bstar_cur:.2e}/{Bstar_init:.2e}".format(
    #     Astar_cur = Astars[frame_index], Astar_init = Astars[0],
    #     Bstar_cur = Bstars[frame_index], Bstar_init = Bstars[0]
    # ))
    Astars[frame_index] = np.sum(T[[0,2],:,:])
    Bstars[frame_index] = np.sum(T[[1,2],:,:])

    axes[-1].set_title("[A*]={Astar_cur:.2e}/{Astar_init:.2e}, [B*]={Bstar_cur:.2e}/{Bstar_init:.2e}".format(
        Astar_cur = Astars[frame_index], Astar_init = Astars[0],
        Bstar_cur = Bstars[frame_index], Bstar_init = Bstars[0]
    ))
    A_line.set_xdata(np.linspace(0, dt*frame_index, frame_index))
    A_line.set_ydata(Astars[:frame_index])
    B_line.set_xdata(np.linspace(0, dt*frame_index, frame_index))
    B_line.set_ydata(Bstars[:frame_index])


fig, axes = plt.subplots(1, 4, figsize=(10, 10))
mats = []

## initialize the plot
# initialize animation plot
for index in range(len(species)):
    axes[index].set_xlim(0, T_size[0])
    axes[index].set_ylim(0, T_size[1])
    axes[index].set_title(f"Concentration of Species {index2name[index]} at {0:.2f}s")
    axes[index].set_aspect("equal")
    axes[index].axis("off")
    mat = axes[index].matshow(
        T[index, :, :],
        cmap = colormaps[index],
        vmin = cmap_min,
        vmax = cmap_max,
    )
    mats.append(mat)

# initialize analysis plot
axes[-1].set_xlim(0, total_frame * dt * 1.1)
axes[-1].set_ylim(
    0.5*np.min([Astars[0], Bstars[0]]), 
    1.2*np.max([Astars[0], Bstars[0]]))
axes[-1].set_aspect("auto")
axes[-1].set_title("[A*]={Astar_cur:.2e}/{Astar_init:.2e}, [B*]={Bstar_cur:.2e}/{Bstar_init:.2e}".format(
    Astar_cur = Astars[0], Astar_init = Astars[0],
    Bstar_cur = Bstars[0], Bstar_init = Bstars[0]
))
A_line, = axes[-1].plot(range(1), Astars[0], label="A*")
B_line, = axes[-1].plot(range(1), Bstars[0], label="B*")
axes[-1].legend()

## run animation
ani = animation.FuncAnimation(
    fig, 
    update, 
    frames=total_frame, 
    interval=interval, 
    repeat = False,
    save_count=50)

if save_animation:
    # save the animation
    ani.save('BimolecularReactionDiffusion.mp4', writer='ffmpeg', fps=30)
else:
    plt.show()




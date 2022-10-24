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
k2 = 1.
# diffusion coefficients of A, B, AB at x and y direction, a value in [0, 1]
D_A = 1.
D_B = 1.
D_AB = 1.
# external force
#TODO: implement this feature
# velocity field


# Define a 2d rectangle T. We will wrap it into a torus in the algorithm
# T is a 3-dim array. Fix time t, T[0, x1, x2] is the concentration of species[0] 
#   at location (x1, x2)
T_size = (5, 5)
dx = 0.1   # meter, the gap of the grid 
T = np.zeros((len(species), *T_size), dtype=float)
# T[0, 1:3, 1:3] = np.ones((2, 2))
T[0, 1, 1] = 1
T[1, 2, 2] = 1
# T[0, 6:11, 6:11] = 10*np.ones((5,5))
# T[1, 30:70, 30:70] = 5*np.ones((40,40))
# T[2, 7:12, 7:12] = 50*np.random.random((5,5))
# T[0,:50,:50] = 5*np.ones((50,50))
# T[1,30:80,30:80] = 7*np.ones((50,50))
# T[2,3:13,3:13] = 0*np.ones((10,10))

# b = np.ones((2, *T_size))   # velocity field in R2 on the plane T
b = 0.1*np.stack([np.zeros(T_size), np.ones(T_size)])
# b = np.zeros((2,*T_size))


## Animation setting
save_animation = False
colormaps = ["Reds", "Blues", "Purples"]      # colormap for each species

# cmap_min = -np.max(T)     # use this setting to visually check for negative concentration (numerical instability)
# cmap_min = np.min(T)
# cmap_max = np.max(T)
cmap_min = -1
cmap_max = 2

## Numerical solution setting
dt = 0.01   # second
interval = 10   # interval between frames, milisecond
playback_speed = dt * 1000 / interval
total_frame = 500
print(f"dt = {dt}s, Playback speed: x{playback_speed:.1f}, Total time: {dt*total_frame:.2f}s")

Astars = np.zeros(total_frame)
Bstars = np.zeros(total_frame)
Astars[0] = np.sum(T[[0,2],:,:])
Bstars[0] = np.sum(T[[1,2],:,:])

laplacian_matrix = np.array([     
    [1, 1, 1],
    [1, -8, 1],
    [1, 1, 1],
])


def update(frame_index):
    global T
    # Dummy update
    # T = np.random.random((len(species), *T_size))
    

    # change caused by A + B -> AB
    delta_reaction = -k1*(T[0,:,:]*T[1,:,:]) + k2*T[2,:,:]
    
    delta_diffusion = np.zeros((len(species), *T_size))
    for index in range(len(species)):
        # delta_diffusion[index,:,:] = convolve(T[index,:,:], laplacian_matrix)   # square
        delta_diffusion[index,:,:] = convolve(T[index,:,:], laplacian_matrix, mode="wrap")    # torus

    delta_convection = np.zeros((len(species), *T_size))
    ## Non-compressible, gradient-based
    # for index in range(len(species)):
    #     # gradient of current species
    #     grad = np.stack(
    #             np.gradient(np.pad(T[index, :, :], 1, mode="wrap"), edge_order=1),   # velocity_field dot gradient
    #             axis = 0)[:,1:-1,1:-1] / dx
    #     # velocity field dot product gradient
    #     delta_convection[index, :, :] = np.sum(b * grad, axis=0)
    
    ## Direct discretization of divergence
    # bu = np.zeros((len(species), *T_size))
    # for index in range(len(species)):
    #     bu = np.pad(b * np.stack([T[index,:,:], T[index,:,:]]), 1, mode="wrap")
    #     delta_convection[index, :, :]  =  1/8*(bu[0, 2:, 1:-1] + bu[0, :-2, 1:-1] + bu[1, 1:-1, 2:] + bu[1, 1:-1, :-2])

    T[0,:,:] += ( delta_reaction + D_A  * delta_diffusion[0, :, :] - delta_convection[0,:,:]) * dt
    T[1,:,:] += ( delta_reaction + D_B  * delta_diffusion[1, :, :] - delta_convection[1,:,:]) * dt
    T[2,:,:] += (-delta_reaction + D_AB * delta_diffusion[2, :, :] - delta_convection[2,:,:]) * dt
    

    for index in range(len(species)):
        # axes[index].set_title(f"Concentration of Species {index2name[index]} at {dt*frame_index:.2f}s, avg={np.average(T[index, :, :]):.2f}")
        axes[index].set_title(f"[{index2name[index]}] at {dt*frame_index:.2f}s, min={np.min(T[index, :, :]):.2e}")
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
fig.patch.set_facecolor('xkcd:gray')
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




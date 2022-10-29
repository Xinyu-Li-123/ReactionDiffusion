from calendar import c
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage import convolve
import warnings


"""
Bimolecular reaction with diffusion on a 1m*1m torus
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
warnings.filterwarnings("error")
np.random.seed(42)

# binary cmap where zero is black and nonzero is white
# binary_cmap = matplotlib.colors.ListedColormap(['black', 'white'])
binary_cmap = matplotlib.cm.get_cmap("gray")

# cmap_A = matplotlib.cm.get_cmap("Reds")
# cmap_B = matplotlib.cm.get_cmap("Blues")
cmap_A = binary_cmap
cmap_B = binary_cmap

def initialize():
    ## Initialize
    species = np.zeros(3, dtype=float)   # [A, B, AB]
    index2name = ["A", "B", "AB"]
    # reaction coefficients
    k1 = .23
    k2 = .6
    # diffusion coefficients of A, B, AB at x and y direction, a value in [0, 1]
    D_A = 0.0002
    D_B = 0.0002
    D_AB = 0.0002
    # feed / kill species A / B / AB at certain time
    #TODO: implement this feature

    # Define a 2d rectangle T. We will wrap it into a torus in the algorithm
    # T is a 3-dim array. Fix time t, T[0, x1, x2] is the concentration of species[0] 
    #   at location (x1, x2)
    T_size = (200, 200)
    h = 0.0029       # gap of the grid, the smaller h, the more refined the grid is.
    T = np.zeros((len(species), *T_size), dtype=float)
    # T[0, :70, :70] = 5*np.random.random((70,70))
    # T[1, 30:, 30:] = 5*np.random.random((70,70))
    # T[0, 20:50, 20:50] = 10*np.ones((30,30))
    # T[1, 30:60, 30:60] = 10*np.ones((30,30))
    
    # # random
    # T[0, :, :] = 1*np.random.random(T_size)
    # T[1, :, :] = 1*np.random.random(T_size)    


    # # Tetris
    # T[0,40:80,40:120] = 1
    # T[0,80:120, 40:80] = 1
    # T[0,120:160, 40:120] = 1

    # T[1,40:80,120:160] = 1
    # T[1,80:120, 80:160] = 1
    # T[1,120:160, 120:160] = 1

    # # Long strip
    # T[0, 40:80, :] = 1
    # T[1, 70:110, :] = 1

    # Overlapping square
    T[0, 50:130, 40:120] = 1
    T[1, 70:150, 70:150] = 1


    # T[2,3:13,3:13] = 0*np.ones((10,10))

    # T[0,:,:] = 1*np.random.random(T_size)
    # T[1,:,:] = 1*np.random.random(T_size)

    # T[0,:,:] = 0
    # T[1,:,:] = 0

    # generate checkboard pattern with blocks of size 50*50
    # n = 20
    # for i in range(0, T_size[0], n):
    #     for j in range(0, T_size[1], n):
    #         T[0, i:i+n//2, j:j+n//2] = 1
    #         T[1, i+n//2:i+n, j+n//2:j+n] = 1
            
    # generate two concentric rings on T of different sizes
    # x, y = np.ogrid[0:T_size[0], 0:T_size[1]]
    # mask = (x - T_size[0]/2)**2 + (y - T_size[1]/2)**2 < (T_size[0]/2)**2
    # mask = mask & ((x - T_size[0]/2)**2 + (y - T_size[1]/2)**2 > (T_size[0]/4)**2)
    # T[0, mask] = 1
    # mask = (x - T_size[0]/2)**2 + (y - T_size[1]/2)**2 < (T_size[0]/3)**2
    # T[1, mask] = 1


    ## Animation setting
    save_animation = False
    animation_name = "ReactionDiffusion_A_B_concentric"
    colormaps = ["Reds", "Blues", "Purples"]      # colormap for each species
    binary_cmap = matplotlib.colors.ListedColormap(['black', 'white'])

    # cmap_min = -np.max(T)     # use this setting to visually check for negative concentration (numerical instability)
    # cmap_min = np.min(T)
    # cmap_max = np.max(T)
    # cmap_min = 0
    # cmap_max = 5


    ## Numerical solution setting
    dt = 0.01   # second
    interval = 10   # interval between frames, milisecond
    playback_speed = dt * 1000 / interval
    total_frame = 2000

    Astars = np.zeros(total_frame)
    Bstars = np.zeros(total_frame)
    Astars[0] = np.sum(T[[0,2],:,:])
    Bstars[0] = np.sum(T[[1,2],:,:])

    laplacian_matrix = 1 / h**2 * np.array([     
        [0, 1, 0],
        [1, -4, 1],
        [0, 1, 0],
    ])
    
    return species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix


def simulate_without_animation(species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix):
    """
    return if has overflow
    """
    # global species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix
    
    def update_without_animation(frame_index):

        # change caused by A + B -> AB
        try:
            delta_reaction = -k1*(T[0,:,:]*T[1,:,:]) + k2*T[2,:,:]
            delta_diffusion = np.zeros((len(species), *T_size))
            for index in range(len(species)):
                delta_diffusion[index,:,:] = convolve(T[index,:,:], laplacian_matrix, mode="wrap")

            T[0,:,:] += ( delta_reaction + D_A  * delta_diffusion[0, :, :]) * dt
            T[1,:,:] += ( delta_reaction + D_B  * delta_diffusion[1, :, :]) * dt
            T[2,:,:] += (-delta_reaction + D_AB * delta_diffusion[2, :, :]) * dt

        except RuntimeWarning as e:
            print(e)
            return True
        

    print("Simulating without animation...")    

    has_overflow = False

    fig, axes = plt.subplots(2, 2, figsize=(7, 7))
    # set background color to gray
    fig.patch.set_facecolor('gray')
    axes[0, 0].matshow(
        T[0, :, :],
        cmap = cmap_A
        )
    axes[0, 0].set_title("[A] at {time:.2f}s".format(time=0))
    axes[0, 1].matshow(
        T[1, :, :],
        cmap = cmap_B
        )
    axes[0, 1].set_title("[B] at {time:.2f}s".format(time=0))


    for frame_index in range(0, total_frame):
        if has_overflow:
            print("Overflow detected, aborting simulation...")
            return True

        if frame_index % 1000 == 0:
            print(f"frame {frame_index} / {total_frame}")
            print("{:.4f}, {:.4f}".format(T[0,:,:].min(), T[0,:,:].max()))

        has_overflow = update_without_animation(frame_index)

    axes[1, 0].matshow(
        T[0,:,:], 
        cmap = cmap_A
        # cmap="Reds"
        )
    axes[1, 0].set_title("[A] at {time:.2f}s".format(time=frame_index*dt))
    axes[1, 1].matshow(
        T[1,:,:], 
        cmap = cmap_B
        # cmap="Blues"
        )
    axes[1, 1].set_title("[B] at {0:.2f}s".format(frame_index*dt))
    plt.show()

    # save as png
    fig.savefig(f"bd_tetris_{str(total_frame).zfill(5)}.png")

    return False


def simulate_with_animation(species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix):
    # same as simulate_without_animation, but with animation
    # global species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix

    def update_with_animation(frame_index):
        # change caused by A + B -> AB
        if frame_index+1 % 1000 == 0:
            print(f"frame {frame_index} / {total_frame}")
        if frame_index % 1000 == 0:
            print(T[0,:,:].min(), T[0,:,:].max())

        try:
            delta_reaction = -k1*(T[0,:,:]*T[1,:,:]) + k2*T[2,:,:]
            delta_diffusion = np.zeros((len(species), *T_size))
            for index in range(len(species)):
                delta_diffusion[index,:,:] = convolve(T[index,:,:], laplacian_matrix, mode="wrap")

            T[0,:,:] += ( delta_reaction + D_A  * delta_diffusion[0, :, :]) * dt
            T[1,:,:] += ( delta_reaction + D_B  * delta_diffusion[1, :, :]) * dt
            T[2,:,:] += (-delta_reaction + D_AB * delta_diffusion[2, :, :]) * dt

        except RuntimeWarning as e:
            print(e)
            return True

        # update species plots
        for index in range(len(species)-1):
            mats[2+index].set_data(T[index, :, :])
            axes[1, index].set_title("[{name}] at {time:.2f}s".format(name=index2name[index], time=frame_index*dt))
            # print(T[index,:,:].min(), T[index,:,:].max())

        return False
    
    fig, axes = plt.subplots(2, 2, figsize=(7, 7))
    # set background color to gray
    fig.patch.set_facecolor('gray')
    mats = []
    for index in range(len(species)-1):
        mat = axes[0, index].matshow(
            T[index, :, :],
            cmap = cmap_A if index == 0 else cmap_B
            )
        mats.append(mat)
    axes[0, 0].set_title("[A] at {time:.2f}s".format(time=0))
    axes[0, 1].set_title("[B] at {time:.2f}s".format(time=0))    
    for index in range(len(species)-1):
        mat = axes[1, index].matshow(
            T[index, :, :],
            cmap = cmap_A if index == 0 else cmap_B
            )
        mats.append(mat)
    
    # animation
    print("Simulating with animation...")
    ani = animation.FuncAnimation(
        fig, 
        update_with_animation, 
        frames=total_frame, 
        interval=interval, 
        repeat = False,
        save_count=50)
    plt.show()


species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix = initialize()
dt = dt
total_frame = 1000

h = h
laplacian_matrix = 1 / h**2 * np.array([     
    [0, 1, 0],
    [1, -4, 1],
    [0, 1, 0],
])

D = D_A
print(f"dt = {dt}, h = {h}, CFD = {D*dt/h**2}")
# simulate_with_animation(species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix)
simulate_without_animation(species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix)


# fig, axes = plt.subplots(1, 2, figsize=(10, 10))
# axes[0].plot(np.linspace(0, dt*total_frame, total_frame), Astars, label="[A*]")
# axes[0].set_title("Time evolution of [A*], std = {:.2e}".format(np.std(Astars)))
# axes[0].set_xlabel("Time (s)")
# axes[0].set_ylabel("[A*]")
# axes[0].set_aspect("auto")
# axes[0].legend()

# axes[1].plot(np.linspace(0, dt*total_frame, total_frame), Bstars, label="[B*]")
# axes[1].set_title("Time evolution of [B*], std = {:.2e}".format(np.std(Bstars)))
# axes[1].set_xlabel("Time (s)")
# axes[1].set_ylabel("[B*]")
# axes[1].set_aspect("auto")
# axes[1].legend()

# with open("std.csv", "w") as file:
#     # standard deviation of [A*] and [B*] at different time period of size 2000 
#     file.write("t0,t1,std([A*]),std([B*])\n")
#     for i in range(0, total_frame, 8000):
#         file.write(f"{i*dt:.2f},{(i+8000)*dt:.2f},{np.std(Astars[i:i+8000]):.2e},{np.std(Bstars[i:i+8000]):.2e}\n")

# plot the standard deviation of A, B and AB in one figure
# fig, axes = plt.subplots(1, 1, figsize=(10, 10))
# axes.plot(np.linspace(0, dt*total_frame, total_frame), std_species[0,:], label="[A]")
# axes.plot(np.linspace(0, dt*total_frame, total_frame), std_species[1,:], label="[B]")
# axes.plot(np.linspace(0, dt*total_frame, total_frame), std_species[2,:], label="[AB]")
# axes.set_title("Time evolution of standard deviation of [A], [B] and [AB]")
# axes.set_xlabel("Time (s)")
# axes.set_ylabel("Standard deviation")
# axes.set_aspect("auto")
# axes.legend()



# standard deviation of [A*] and [B*] with respect to dt
# fig, axes = plt.subplots(1, 2, figsize=(10, 10))
# axes[0].plot(dt_range, std_Astars, label="std([A*])")
# axes[0].set_title("std([A*])")
# axes[0].set_xlabel("dt")
# axes[0].set_ylabel("std([A*])")
# axes[0].legend()

# axes[1].plot(dt_range, std_Bstars, label="std([B*])")
# axes[1].set_title("std([B*])")
# axes[1].set_xlabel("dt")
# axes[1].set_ylabel("std([B*])")
# axes[1].legend()
# plt.show()
# simulate_without_animation()

# simulate_no_animation()
# print(f"In the end, average of [A] = {np.mean(T[0,:,:])}")
# print(f"In the end, average of [B] = {np.mean(T[1,:,:])}")


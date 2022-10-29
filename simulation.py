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

def initialize():
    ## Initialize
    species = np.zeros(3, dtype=float)   # [A, B, AB]
    index2name = ["A", "B", "AB"]
    # reaction coefficients
    k1 = .23
    k2 = .23
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
    h = 0.0027       # gap of the grid, the smaller h, the more refined the grid is.
    T = np.zeros((len(species), *T_size), dtype=float)
    # T[0, :70, :70] = 5*np.random.random((70,70))
    # T[1, 30:, 30:] = 5*np.random.random((70,70))
    # T[0, 20:50, 20:50] = 10*np.ones((30,30))
    # T[1, 30:60, 30:60] = 10*np.ones((30,30))
    # T[0,100:150,100:150] = 50*np.ones((50,50))
    # T[1,200:250,200:250] = 70*np.ones((50,50))
    # T[2,3:13,3:13] = 0*np.ones((10,10))

    T[0,:,:] = 2*np.random.random(T_size)
    T[1,:,:] = 2*np.random.random(T_size)
    # T[2,:,:] = 0.05*np.random.random(T_size)

    # generate checkboard pattern with blocks of size 50*50
    # n = 20
    # for i in range(0, T_size[0], n):
    #     for j in range(0, T_size[1], n):
    #         T[0, i:i+n//2, j:j+n//2] = 1
    #         T[1, i+n//2:i+n, j+n//2:j+n] = 1
            
    # generate two concentric rings on T of different sizes
    x, y = np.ogrid[0:T_size[0], 0:T_size[1]]
    mask = (x - T_size[0]/2)**2 + (y - T_size[1]/2)**2 < (T_size[0]/2)**2
    mask = mask & ((x - T_size[0]/2)**2 + (y - T_size[1]/2)**2 > (T_size[0]/4)**2)
    T[0, mask] = 1
    mask = (x - T_size[0]/2)**2 + (y - T_size[1]/2)**2 < (T_size[0]/3)**2
    T[1, mask] = 1


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


def simulate_with_animation(species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix):
    # global species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix
    def update_with_animation(frame_index):
        # Dummy update
        # T = np.random.random((len(species), *T_size))
        if frame_index % 100 == 0:
            print(f"Frame {frame_index}/{total_frame}")

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
            if index != 2:
                mats[index].set_data(T[index, :, :])
            else:
                mats[index].set_data(T[index, :, :]*1000 > 0)


        axes[-1].set_title("[A*]={Astar_cur:.2e}/{Astar_init:.2e}, [B*]={Bstar_cur:.2e}/{Bstar_init:.2e}".format(
            Astar_cur = Astars[frame_index], Astar_init = Astars[0],
            Bstar_cur = Bstars[frame_index], Bstar_init = Bstars[0]
        ))
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
    # set background color to gray
    fig.patch.set_facecolor('gray')
    mats = []

    ## initialize the plot
    # initialize animation plot
    for index in range(len(species)):
        axes[index].set_xlim(0, T_size[0])
        axes[index].set_ylim(0, T_size[1])
        axes[index].set_title(f"[{index2name[index]}] at {0:.2f}s")
        axes[index].set_aspect("equal")
        axes[index].axis("off")
        mat = axes[index].matshow(
            T[index, :, :],
            # cmap = colormaps[index],
            # vmin = cmap_min,
            # vmax = cmap_max,
        )
        mat.set_cmap(colormaps[index])
        # mat.set_cmap(binary_cmap)
        mats.append(mat)
    # mats[2].set_data((T[0, :, :] < T[1, :, :])*10)
    # mats[2].set_cmap(binary_cmap)
    # print(T[0,:,:] < T[1,:,:])
    # sys.exit()


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
        update_with_animation, 
        frames=total_frame, 
        interval=interval, 
        repeat = False,
        save_count=50)

    if save_animation:
        # save the animation
        ani.save(f'{animation_name}.mp4', writer='ffmpeg', fps=30)
    else:
        plt.show()


def simulate_without_animation(species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix):
    """
    return if has overflow
    """
    # global species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix
    
    def update_without_animation(frame_index):
        # Dummy update
        # T = np.random.random((len(species), *T_size))
        

        # change caused by A + B -> AB
        try:
            delta_reaction = -k1*(T[0,:,:]*T[1,:,:]) + k2*T[2,:,:]
            delta_diffusion = np.zeros((len(species), *T_size))
            for index in range(len(species)):
                delta_diffusion[index,:,:] = convolve(T[index,:,:], laplacian_matrix, mode="wrap")

            T[0,:,:] += ( delta_reaction + D_A  * delta_diffusion[0, :, :]) * dt
            T[1,:,:] += ( delta_reaction + D_B  * delta_diffusion[1, :, :]) * dt
            T[2,:,:] += (-delta_reaction + D_AB * delta_diffusion[2, :, :]) * dt

            # std_species[0, frame_index] = T[0, :, :].std()
            # std_species[1, frame_index] = T[1, :, :].std()
            # std_species[2, frame_index] = T[2, :, :].std()

            # Astars[frame_index] = np.sum(T[[0,2],:,:])
            # Bstars[frame_index] = np.sum(T[[1,2],:,:])
        except RuntimeWarning as e:
            print(e)
            return True
        

    # print("Simulating without animation...")    

    # fig, axes = plt.subplots(1, 2, figsize=(10, 10))
    # axes[0].set_title("CDF of A")
    # # axes[0].set_xlim(0, np.max(T[0,:,:])*1.1)
    # # axes[0].set_ylim(0, 1)
    # # axes[0].set_aspect("auto")

    # axes[1].set_title("CDF of B")
    # # axes[1].set_xlim(0, np.max(T[1,:,:])*1.1)
    # # axes[1].set_ylim(0, 1)
    # # axes[1].set_aspect("auto")

    has_overflow = False
    for frame_index in range(0, total_frame):
        if has_overflow:
            print("Overflow detected, aborting simulation...")
            return True

        if frame_index+1 % 1000 == 0:
            print(f"frame {frame_index} / {total_frame}")
        if frame_index % 1000 == 0:
            print(T[0,:,:].min(), T[0,:,:].max())

        has_overflow = update_without_animation(frame_index)

        # # plot CDF of A and B on axes[0] and axes[1] with label
        # if frame_index > 19999 and frame_index % 5000 == 0:
        #     for species_index in range(2):
        #         # The calculation of CDF is based on the following link:
        #         #    https://www.geeksforgeeks.org/how-to-calculate-and-plot-a-cumulative-distribution-function-with-matplotlib-in-python/
        #         count, bins_count = np.histogram(T[0,:,:].flatten(), bins=100)
                
        #         # finding the PDF of the histogram using count values
        #         pdf = count / sum(count)
                
        #         # using numpy np.cumsum to calculate the CDF
        #         # We can also find using the PDF values by looping and adding
        #         cdf = np.cumsum(pdf)
                
        #         # plotting PDF and CDF
        #         # axes[species_index].plot(bins_count[1:], cdf, label=f"CDF at {frame_index*dt:.2f}s", linewidth=2)
    # print("mu_A={:.2e}, mu_B={:.2e}, std([A*]): {:.2e}, std([B*]): {:.2e})".format(
    #     D_A*dt/(h**2), D_B*dt/(h**2),
    #     np.std(T[0,:,:]), np.std(T[1,:,:])
    # ))


    # axes[0].legend()
    # axes[1].legend()

    # axes[0].hist(T[0,:,:].flatten(), bins=100)
    # axes[0].set_title("Histogram of [A]")
    # axes[1].hist(T[1,:,:].flatten(), bins=100)
    # axes[1].set_title("Histogram of [B]")

    # axes[0].matshow(
    #     T[0,:,:], 
    #     # cmap="Reds"
    #     )
    # axes[0].set_title("Concentration of [A]")
    # axes[1].matshow(
    #     T[1,:,:], 
    #     # cmap="Blues"
    #     )
    # axes[1].set_title("Concentration of [B]")
    # plt.show()
    return False


# std_Astars = []
# std_Bstars = []
dt_range = np.linspace(0.0001, 0.01, 10)
h_range = np.linspace(0.0001, 0.003, 10)
dt_mesh, h_mesh = np.meshgrid(dt_range, h_range)

# a dictionary that maps (dt, h)
dt_h_has_overflow = {}
for dt in dt_range:
    for h in h_range:
        dt_h_has_overflow[(dt, h)] = False


# for _h in h_range:
#     for _dt in dt_range:
for dt_index in range(len(dt_range)):
    for h_index in range(len(h_range)):
        species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix = initialize()
        dt = dt_range[dt_index]
        total_frame = 200
        # Astars = np.zeros(total_frame)
        # Bstars = np.zeros(total_frame)
        
        h = h_range[h_index]
        laplacian_matrix = 1 / h**2 * np.array([     
            [0, 1, 0],
            [1, -4, 1],
            [0, 1, 0],
        ])

        # std_species = np.zeros((len(species), total_frame))
        # print(f"A grid of size {T_size[0]*h:.2f}m*{T_size[1]*h:.5f}m cutted into {T_size[0]}*{T_size[1]} blocks, gap h = {h:.2f}m")
        # print(f"dt = {dt:.5f}s, Playback speed: x{playback_speed:.1f}, Total time: {dt*total_frame:.2f}s")
        # print(f"h={h:.4f}, dt={dt:.4f}, CFD number is A={D_A*dt/h**2:.5f}, B={D_B*dt/h**2:.5f}, AB={D_AB*dt/h**2:.5f}")
        
        dt_h_has_overflow[(dt, h)] = simulate_without_animation(species, index2name, k1, k2, D_A, D_B, D_AB, T_size, h, T, save_animation, animation_name, colormaps, binary_cmap, dt, interval, playback_speed, total_frame, Astars, Bstars, laplacian_matrix)

fig, ax = plt.subplots(1, 1, figsize=(4, 3))

# plot scatter of dt and h
for dt_index in range(len(dt_range)):
    for h_index in range(len(h_range)):
        dt = dt_range[dt_index]
        h = h_range[h_index]
        ax.scatter(dt, h, c="red" if dt_h_has_overflow[(dt, h)] else "green")

for key in dt_h_has_overflow:
    mu = D_A * key[0] / key[1]**2
    dt_h_has_overflow[key] = (dt_h_has_overflow[key], mu)
print(dt_h_has_overflow)

# # plot the boundary line
D = D_A
f = lambda dt: (D*dt/0.25)**0.5
ax.plot(np.linspace(0, 0.01, 100), f(np.linspace(0, 0.01, 100)), c='black', label='D*dt/h^2 = 0.25')

ax.set_xlabel('dt')
ax.set_ylabel('h')
ax.set_xlim(0, 0.011)
ax.set_ylim(0, 0.0035)
ax.legend()
ax.set_title('dt vs h')
plt.show()


        # print(f"std([A*]) = {np.std(Astars):.2e}, std([B*]) = {np.std(Bstars):.2e}")
#         std_Astars.append(np.std(Astars))
#         std_Bstars.append(np.std(Bstars))

# print(std_Astars)
# print(std_Bstars)



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


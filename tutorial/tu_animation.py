import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# plot a cellular automation

fig, ax = plt.subplots(1, 1)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)

# initialize the cellular automation
# 0: dead, 1: alive
CA = np.zeros((100, 100), dtype=int)
CA[45:56, 45:56] = 1

# define the update function
def update(data):
    global CA
    # for i in range(1, CA.shape[0]-1):
    #     for j in range(1, CA.shape[1]-1):
    #         if CA[i, j] == 1:
    #             if np.sum(CA[i-1:i+2, j-1:j+2]) < 3:
    #                 CA[i, j] = 0
    #         else:
    #             if np.sum(CA[i-1:i+2, j-1:j+2]) == 3:
    #                 CA[i, j] = 1
    I = CA.shape[0]-1
    J = CA.shape[1]-1
    # mod indices by I and J to simulate a torus
    for i in range(0, I):
        for j in range(0, J):
            if CA[i%I, j%J] == 1:
                if np.sum(CA[
                    ((i-1)%I , i%I    , (i+1)%I, 
                     (i-1)%I , i%I    , (i+1)%I, 
                     (i-1)%I , i%I    , (i+1)%I ),
                    ((j-1)%J , (j-1)%J, (j-1)%J,  
                      j%J    , j%J    , j%J    , 
                      (j+1)%J, (j+1)%J, (j+1)%J )
                    ]) < 3:
                    CA[i%I, j%J] = 0
            else:
                if np.sum(CA[
                    ((i-1)%I , i%I    , (i+1)%I, 
                     (i-1)%I , i%I    , (i+1)%I, 
                     (i-1)%I , i%I    , (i+1)%I ),
                    ((j-1)%J , (j-1)%J, (j-1)%J,  
                      j%J    , j%J    , j%J    , 
                      (j+1)%J, (j+1)%J, (j+1)%J )
                    ]) == 3:
                    CA[i%I, j%J] = 1
    mat.set_data(CA)
    return [mat]

# define the animation
mat = ax.matshow(CA)
ani = animation.FuncAnimation(fig, update, frames=1000, interval=50, save_count=50)
plt.show()

# save the animation
# ani.save('CA.mp4', writer='ffmpeg', fps=30)


# Development Log
[2022/10/23 00:03]
# [Fixed] Bug: Negative concentration

For some reason, there is negative concentration when using Discrete Laplacian operator to simulate diffusion of a single substance. This cause the entire simulation to crash.

* Solution
<del>By negating the laplacian term (practically, negative diffusion coefficient), I fix this problem. But I have no idea why that would work. Theoretically, it is NOT diffufion euqation if you add this negative sign.</del>

[Update]
The real reason is that I write Laplacian matrix wrong. I write the negation of the Laplacian matrix...


# Bug: Large initial concentration / diffusion coefficient / Laplacian matrix cause numerical instability

Large value of inital concentration / diffusion or Laplacian of large size will cause numerical instability.

The numerical instability manifests in the form of negative concentration (which is physically impossible in our setting). 


# [Fixed] Bug: Forget about external forces

The gas seems to be too still. It only diffuses within a limited radius.

This is because we ignore all forms of forces (gravity, levitation, etc)

Possible solution: [Convection-diffusion equation](https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation)

Add divergence to the differential equation, so that the convection between gases is also considered

# Bug: Convection not working, negative concentration again.

I add a constant velocity field of (1, 1) to the torus. However, the negative convection problem occurs during the simulation again. Although it's a bit different this time. Th
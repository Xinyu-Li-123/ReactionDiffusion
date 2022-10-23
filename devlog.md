# Development Log
[2022/10/23 00:03]
# Bug: Negative concentration

For some reason, there is negative concentration when using Discrete Laplacian operator to simulate diffusion of a single substance. This cause the entire simulation to crash.

* Solution
By negating the laplacian term (practically, negative diffusion coefficient), I fix this problem. But I have no idea why that would work. Theoretically, it is NOT diffufion euqation if you add this negative sign.

[Update]
The real reason is that I write Laplacian matrix wrong. I write the negation of the Laplacian matrix...


# Bug: Large initial concentration cause numerical instability
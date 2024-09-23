# CoSpeciation


A two-patch, sexually reproducing individual-based model with nonoverlapping generations for investigating ecological and non-ecological factors of the speciation process.  

## Fitness function

$$
W_i = \frac{\lambda \psi(z_i, z_j)}{1 + \frac{C(z_i, z_j)}{K(z_i)}} (1-I_i)
$$

- The fitness of individual i $`W_i`$ is the product of its fertility $\lambda$, reduction in fitness due to incompatibility of its parents $`I_i`$, and its ability to find a mate $`\psi(z_i, z_j)`$, divided by the amount of competition the individual experiences $`(1 + \frac{C(z_i, z_j}{K(z_i})`$
- For simplicity, we assume $`\lambda`$ is equal for all individuals, hence the removal of the subscript.

## Incompatibility

$$
I_i = 1 - (1 - \phi)^{n_i}
$$

- The fitness effect due to the incompatibility of individual i's parents $`\I_i`$ increases geometrically with the number of incompatible BDMI loci $`\n_i`$ between the parents by the factor $`1 - \phi`$.

##

# CoSpeciation


A two-patch, sexually reproducing individual-based model with nonoverlapping generations for investigating ecological and non-ecological factors of the speciation process.  

## Fitness function

$$
W_i = \frac{\lambda \psi(z_i, z_j)}{1 + \frac{C(z_i, z_j)}{K(z_i)}} (1-I_i)
$$

- The fitness of individual i ($`W_i`$) is the product of its fertility $\lambda$, reduction in fitness due to incompatibility of its parents $`I_i`$, and its ability to find a mate $`\psi(z_i, z_j)`$, divided by the amount of competition the individual experiences $`1 + \frac{C(z_i, z_j)}{K(z_i)}`$
- For simplicity, we assume $`\lambda`$ is equal for all individuals, hence the removal of the subscript.

## Incompatibility

$$
I_i = 1 - (1 - \phi)^{n_i}
$$

- The fitness effect due to the incompatibility of individual i's parents $`I_i`$ increases geometrically with the number of incompatible BDMI loci $`n_i`$ between the parents by the factor $`1 - \phi`$.

- Individuals are haploids with $B = 5$ unlinked BDMI loci, with 2 alleles (1 or 0) each. Number of incompatible BDMI loci $n_i$ is computed as the number of mismatched BDMIs between mated pairs.


## Mating probability

-   The probability of female $`i`$ mating with male $`j`$ (e.g., $`\psi (z_i, z_j)`$) is determined by the distance between their ecological trait $`z`$ and the strength of assortative mating ($`\sigma_m`$) normalized over potential male mates of female $`i`$ (i.e., all females mate with probability = 1 given that males are present) according to the equation

$$
\psi (z_i, z_j)  = \frac{e^{\frac{(z_i - z_j)^2}{2\sigma_m^2}}} {\sum_j{e^{\frac{(z_i - z_j)^2}{2\sigma_m^2}}}}
$$

-   This model thus incorporates mate competition such that at each time step, some males may not mate, while some may mate with multiple females
-   **NOT YET INCORPORATED:** $\sigma_m$ currently fixed across individuals

## Resource competition

-   Individuals in a population compete for resources according to their differences in ecological trait $`z_i`$ (e.g., $`C_i(z_i, z_j)`$) scaled by resources available to the individual $`K(z_i)`$ (i.e., $`C(z_i, z_j) / K(z_i)`$)

-   Specifically, $`C(z_i, z_j)`$ is the competition experienced by individual $`i`$ in the presence of individual $`j`$s defined as,

$$
C_{ij} = exp\frac{-\sum_{i \neq j}(z_i - z_j)^2}{2\sigma_z^2}
$$

where $`\sigma_z^2`$ is the competition kernel (or niche width) which we assume is fixed and equal for all individuals

-   In the absence of competitors, an individual with phenotype $`z_i`$ has access to $`K(z_i)`$ amount of resources which is defined as,

$$
K(z_i) =   K_oexp(\frac{-(z_i - \theta_p)^2}{2\sigma_k^2})
$$

where $`K_o`$ is the maximum amount of resources, or the resources available to an individual with an optimum phenotype for that patch, $`\theta_p`$, and $`\sigma_k^2`$ determines the rate in which resource availability decays with the distances between $`z_i`$ and the optimum (i.e., strength of stabilizing selection).

## Inheritance and mutations

-   **BDMIs**: Offsprings inherit each BDMI locus independently via random sampling of parent alleles. Each inherited allele has a probability of flipping $`\mu_i`$

-   **Assortative mating** $`\sigma_m`$: *currently not incorporated yet*

-   **Ecological trait** $`z_i`$: offspring $`z`$ is drawn from a Gaussian distribution with mean equal the average of parent $`z_i`$ and variance, $`\mu_z`$ (i.e., $`\mathcal{N}(z_{parents}, \mu_z)`$)

-   **Location** $`p_i`$: offspring are born into the patch of its parents with probability $`1 - m`$ and migrate to the other patch with probability $`m`$


## R functions

1. `initial_condition(N = 50, bdmi_B = 5)`

- this function initializes a datamframe given starting population size of `N` individuals and `bdmi_B` BMDI loci

2. `ibm(population, num_gens)`

- this is the main simulator function
- the argument `population` specifies the initial condition and takes in the output of `initial_condition()` as its input
- `num_gens` specifies the number of generations to run the simulation

   

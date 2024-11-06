---
---
---

# CoSpeciation

A two-patch, sexually reproducing individual-based model with nonoverlapping generations for investigating ecological and non-ecological factors of the speciation process.

## Fitness function

$$
W_i = \frac{\lambda \psi(z_i^p, z_j^p)}{1 + \frac{C(z_i^p, z_j^p)}{K(z_i^p, \theta_p)}} (1- I_i^p)
$$

-   The fitness of individual *i* in patch *p* ($W_i^p$) is the product of its fertility $\lambda$, reduction in fitness due to incompatibility of its parents $I_i^p$, and its ability to find a mate $\psi(z_i^p, z_j^p)$, divided by the amount of competition the individual experiences $1 + \frac{C(z_i^p, z_j^p)}{K(z_i^p, \theta_p)}$
-   For simplicity, we assume $\lambda$ is equal for all individuals, hence the removal of the subscript.

## Incompatibility

$$
I_i = 1 - (1 - \phi)^{n_i}
$$

-   The fitness effect due to the incompatibility of individual *i*'s parents $`I_i`$ increases geometrically with the number of incompatible BDMI loci $`n_i`$ between the parents by the factor $`1 - \phi`$. Higher values of $\phi$ equates to greater fitness reduction per incompatible locus

-   Individuals are haploids with $B = 5$ unlinked BDMI loci, with 2 alleles (1 or 0) each. Number of incompatible BDMI loci $n_i$ is computed as the number of mismatched BDMIs between mated pairs

## Mating probability

-   The probability of female *i* mating with male *j* in patch *p* (e.g., $\psi (z_i^p, z_j^p)$) is determined by the distance between their ecological trait *z* and the tendency of female *i* to mate assortatively ($\sigma_m$) normalized over the number of males in patch *p* (i.e., all females mate with probability = 1 given that males are present in the patch) according to the equation

$$
\psi (z_i^p, z_j^p)  = \frac{e^{\frac{(z_i^p - z_j^p)^2}{2\sigma_m^2}}} {\sum_j{e^{\frac{(z_i - z_j)^2}{2\sigma_m^2}}}}
$$

-   This model thus incorporates mate competition such that at each time step, some males may not mate, while some may mate with multiple females
-   **Assoratative mating** $\sigma_m$ **varies by individuals and can evolve:** what would be a better notation to capture the fact that individuals vary in $\sigma_m$ ? Use different symbols for different variance parameters

## Resource competition

-   The amount of competition that individual $i$ in patch $p$experiences is determined by the difference in ecological trait $z_i^p$ values between it and individual $j$'s in patch $p$ $z_j^p$ (e.g., $C_i^p(z_i^p, z_j^p)$) scaled by the amount of resources available to the individual *i* $K(z_i^p, \theta_p)$ (i.e., $`C(z_i^p, z_j^p, \sigma_z) / K(z_i^p, \theta_p)`$)

-   Specifically, $`C(z_i^p, z_j^p, \sigma_z)`$ is the competition experienced by individual $i$ in the presence of individual $j$s defined as,

$$
C_{ij} = exp\frac{-\sum_{i \neq j}(z_i^p - z_j^p)^2}{2\sigma_z^2}
$$

where $\sigma_z^2$ is the competition kernel (or niche width) which we assume to be fixed and equal for all individuals

-   In the absence of competitors, an individual with phenotype $`z_i`$ has access to $`K(z_i, \theta_p)`$ amount of resources which is defined as,

$$
K(z_i, \theta_p) =   K_oexp(\frac{-(z_i - \theta_p)^2}{2\sigma_K^2})
$$

where $K_o$ is the maximum amount of resources, or the resources available to an individual with an optimum phenotype for patch $p$, $`\theta_p`$, and $`\sigma_K^2`$ determines the rate at which resource availability decays with the distances between $`z_i`$ and the optimum (i.e., strength of stabilizing selection).

## Inheritance and mutations

-   **BDMIs**: Offsprings inherit each BDMI locus independently via random sampling of parent alleles. Each inherited allele has a probability of mutating = $`\mu_b`$

-   **Assortative mating** $\sigma_m$: Offspring $\sigma_m$ is drawn from a Gaussian distribution with mean equal the average of parent $\overline\sigma_m$ and variance, $\mu_m$ (i.e., $\mathcal{N}(\overline\sigma_m, \mu_m)$)

-   **Ecological trait** $z_i$: offspring $`z`$ is drawn from a Gaussian distribution with mean equal the average of parent $\overline z_i$ and variance, $`\mu_z`$ (i.e., $`\mathcal{N}(\overline z_i, \mu_z)`$)

-   **Location** $p$: offspring are born into the patch of its parents with probability $1 - m$ and migrate to the other patch with probability $m$

-   **NOTE**: parameter names are *extremely confusing* when used in the conjunction with probability distributions

## R functions

1.  `initial_condition(N = 50, assort_sigma = 0.01, bdmi_B = 5)`

-   this function initializes a datamframe given starting population size of `N` individuals, each with same tendency to mate assortatively `assort_sigma` and `bdmi_B` BMDI loci

2.  `ibm(population, num_gens)`

-   this is the main simulator function
-   the argument `population` specifies the initial condition and takes in the output of `initial_condition()` as its input
-   `num_gens` specifies the number of generations to run the simulation

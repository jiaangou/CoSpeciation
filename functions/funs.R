#BDMI fitness reduction -----------
bdmi_I <- function(genotype, phi = 0.2, incomp_type = c('10','01'), group_size = 2){

  #Check number of loci is divisible by group_size
  n_loci <- length(genotype)
  if (n_loci %% group_size != 0) {
    stop("Number of loci is not divisible by group_size")
  }

  #Pair up BDMI loci
  bdmi_groups <- sapply(seq(from = 1, to = n_loci, by = group_size), function(i)
    paste0(genotype[i:(i + (group_size-1))], collapse = ''))


  #Count number of incompatible types
  incomp_n  <- sum(bdmi_groups %in% incomp_type)


  #Calculate fitness reduction due to incompatbility
  I <- 1 - (1 - phi)^incomp_n

  return(I)
}

#Initialize data frame -----
initial_condition <- function(N = 50,assort_sigma = 0.01, bdmi_B = 6, phi = 0.2){
  require(dplyr)
  n0 <- data.frame(zi = runif(N), #starting traits are uniformly distributed
                   sex = sample(c("F", "M"), size = N, replace = TRUE),
                   assort_sigma = assort_sigma,
                   incomp = 0,
                   patch = sample(1:2, N, replace = TRUE))%>%
    bind_cols(matrix(rbinom(n = bdmi_B*N, size = 1, 0.9), ncol = bdmi_B)%>%
                as.data.frame()%>%
                setNames(paste0('bdmi',1:bdmi_B)))

  #Calculate incompatbilities
  n0$incomp <- n0[,paste0('bdmi',1:bdmi_B)]%>%
    apply(1, function(x)x%>%
            bdmi_I(phi = phi))

  return(n0)

}


#Mating function -----
mating <- function(population, competition,
                   birth = 5, phi = 0.2){
  require(dplyr)
  #Check that population comes from only a single patch
  patch_i <- unique(population$patch)
  stopifnot(length(patch_i) == 1)

  #Check if both males and females are present
  sex_present <- all(c('F','M') %in% population$sex) #checks that both sex present

  #Name of BDMIs
  bdmi_names <- grep("^bdmi", names(population), value = TRUE) #used for sampling later

  #Split up male and females
  females <- population[population$sex == 'F', ]
  males <- population[population$sex == 'M',]

  #Calculate the fecundity of each female
  female_fecund <- (birth * (1 - females$incomp)) / competition[population$sex == 'F']
  #set.seed(10)
  num_offsprings <- rpois(n = length(female_fecund), lambda = female_fecund)

  #If no offsprings are produced OR if both sexes are not present return NULL
  if(sum(num_offsprings) == 0 | !sex_present){
    message(paste0("Patch ", patch_i, " is extinct"))
    return(NULL)
  }else{


    #Subset females to those that will produce offspring and the number of offsprigns that they will produce
    females <- females[num_offsprings > 0, ]
    num_offsprings <- num_offsprings[num_offsprings > 0]

    #Number of mating events is equal to number of females left after removing 0 births
    num_females <- nrow(females)

    #Create a list with length equal to the number of females
    offsprings <- vector(mode = 'list', length = num_females)


    #Loop through females
    for(f in 1:num_females){

      #Sample a male mate ------
      #compute probabilities of sampling each male given z_i and z_j
      z_i <- females[f, ]$zi
      assort_sigma_i <- females[f, ]$assort_sigma #choosiness of female_i
      z_j <- males$zi
      pij <-(exp((-(z_i - z_j)^2 )/(2*assort_sigma_i)) / sum(exp((-(z_i - z_j)^2 )/(2*assort_sigma_i))))

      #select mate
      male_i <- sample(1:nrow(males), size = 1, prob = pij)

      #parent traits
      pair <- rbind(females[f, ],
                    males[male_i,])

      #Calculate incompatibility ------
      #count the number of loci that are different
      # incomp_n <- pair[,bdmi_names]%>%
      #   apply(2, function(x)abs(diff(x)))%>%
      #   sum
      #offspring_incomp <- 1 - (1 - phi)^incomp_n  #Calculate incomp of offspring


      #Generate offspring phenotypes --------
      offspring_z <- rep(mean(pair$zi), num_offsprings[f]) #offspring z's
      offspring_assort_sigma <- rep(mean(pair$assort_sigma), num_offsprings[f]) #offspring assortativity
      offspring_sex <- sample(c('F','M'), size = num_offsprings[f], replace = TRUE)
      offspring_bdmi <- pair[,bdmi_names]%>%
        apply(2, function(x)sample(x, size = num_offsprings[f], replace = TRUE))


      #BDMI incompatability
      offspring_incomp <- bdmi_I(offspring_bdmi, phi = phi)

      #Convert bdmi to a data.frame
      if(!is.matrix(offspring_bdmi)){
        offspring_bdmi <- offspring_bdmi%>%
          t()%>%
          as.data.frame()
      }else{
        offspring_bdmi <- offspring_bdmi%>%
          as.data.frame()
      }

      #Combine and save
      offsprings[[f]] <- data.frame(zi = offspring_z, sex = offspring_sex,
                                    assort_sigma = offspring_assort_sigma,
                                    incomp = offspring_incomp, patch = patch_i)%>%
        bind_cols(offspring_bdmi)

      #End of female loop
    }
    out <- bind_rows(offsprings)

    return(out)
  }

}


#Simulator  -----
ibm <- function(population, num_gens,
                phi = 0.2, birth = 10,
                theta = c(0, 0.01), sig2_alpha = 0.05, K0 = 50, sig2_k = 0.05,
                migration_p = 0.1, bdmi_mu = 10^-3, eco_mu = 10^-2, assort_mu = 0){

  require(dplyr)

  ###########################
  #1. Setup
  ###########################
  #BDMI names
  bdmi_names <- grep("^bdmi", names(population), value = TRUE)


  #Create output list
  output <- vector('list', length = num_gens + 1)
  output[[1]] <- population #initial condition

  #Stoppers
  g <- 1
  extinction <- nrow(population) == 0

  #Iterate until g = num_gens or extinction occcurs
  while(g <= num_gens && !extinction){

    N <- output[[g]]

    ###########################
    #2. Reproduce by patch
    ###########################
    #Store patch-specific outputs
    patch_results <- vector(mode = 'list', 2)
    #Iterate
    for(p in unique(N$patch)){

      #Resource competition -------------------
      theta_p <- theta[p] #patch optimum
      #Competing individuals
      zi <- N[N$patch == p,]$zi
      #Competitive effects: sum(a_ij)
      c_i <- colSums(exp(-dist(zi, upper = T)^2/(2*sig2_alpha))%>%
                       as.matrix)
      #Abundance of resource_i: (K_i)
      k_i <- K0*exp(-(zi - theta_p)^2/(2*sig2_k))
      #Fitness effect of competition
      comp <- 1 + (c_i/k_i)

      #Print gens + patch
      #print(paste0("Generation:", g, " Patch:", p))
      #print(N)
      #if(g == 2 & p == 2){
      #  return(N)
      #}
      #print("Generation:", g, "Patch:", p)


      #Reproduction -------------------
      patch_results[[p]] <- mating(population = N[N$patch == p,], competition = comp, birth = birth, phi = phi)

    }

    #New generation
    N <- bind_rows(patch_results)

    #Check for extinction
    extinction <- (nrow(N)==0)
    #print(paste0('extinction: ', extinction))


    ###########################
    #3. Mutations
    ###########################
    if(!extinction){
      #eco trait
      N$zi <- rnorm(n = length(N$zi), mean = N$zi, sd = eco_mu)

      #assortativity trait
      N$assort_sigma <- rnorm(n = length(N$assort_sigma), mean = N$assort_sigma, sd = assort_mu)


      #bdmi
      N[,bdmi_names] <- N[,bdmi_names]%>%
        mutate(across(everything(), ~ ifelse(runif(n()) < bdmi_mu, 1 - ., .)))
      #incompatability
      N$incomp <- N[,bdmi_names]%>%
        apply(., 1, function(x)bdmi_I(x, phi = 0.2))


      #location
      N$patch <- ifelse(runif(length(N$patch)) < migration_p, 3 - N$patch, N$patch)



      ###########################
      #4. Update stoppers
      ###########################
      g <- g+1


      ###########################
      #5. Save results
      ###########################
      output[[g]] <- N
    }else{
      #Trim list if extinction occurs before num_gens
      output <- output[1:g]
    }

  }

  return(output)
}


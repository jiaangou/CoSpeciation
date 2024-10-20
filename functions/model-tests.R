#TESTS
library(dplyr)
library(ggplot2)

#Set up initial conditions ----------
N <- 60 #pop size
N0 <- data.frame(zi = runif(N), #starting traits are uniformly distributed
                 sex = sample(c("F", "M"), size = N, replace = TRUE),
                 assort_sigma = 0.05,
                   incomp = 0.1,
                   #patch = '2'
                   patch = sample(1:2, N, replace = TRUE))%>%
  bind_cols(matrix(rbinom(n = 5*N, size = 1, 0.5), ncol = 5)%>%
              as.data.frame()%>%
              setNames(paste0('bdmi',1:5)))
                   #patch = sample(c('1','2'), N, replace = TRUE)



ibm(population = N0, num_gens = 10, birth = 10, eco_mu = 0.1, theta = c(0, 1), migration_p = 0.001)%>%
  bind_rows(.id = 'generations')%>%
  mutate(generations = as.integer(generations))%>%
  mutate(gen = paste0('gen:', generations))%>%
  mutate(patch = as.factor(patch))%>%
  ggplot(aes(x = zi))+
  geom_density(aes(fill = patch), alpha = 0.5)+
  scale_fill_manual(values = c('#386cb0', '#ffff99'))+
  facet_grid(gen~.)+
  theme_bw()
  






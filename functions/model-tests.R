#TESTS
library(dplyr)
library(ggplot2)

#Color palette
col2 <- wesanderson::wes_palette("FrenchDispatch")[1:2]


#Set up initial conditions ----------
N <- 60 #pop size
N0 <- initial_condition(N = N, assort_sigma = 0.05, bdmi_B = 5)


#Demo ---------
sims <- ibm(population = N0, num_gens = 10, birth = 10, eco_mu = 0.1, theta = c(0, 1), migration_p = 0.001)%>%
  bind_rows(.id = 'generations')%>%
  mutate(generations = as.integer(generations))%>%
  mutate(gen = paste0('gen:', generations))%>%
  mutate(patch = as.factor(patch))

#Population-level ---------------
sims %>%
  group_by(patch, generations) %>%
  summarise(zi = mean(zi), n = n(), assort_sigma = mean(assort_sigma)) %>%
  tidyr::pivot_longer(cols = c('zi', 'n', 'assort_sigma'), names_to = 'Variable') %>%
  # Convert Variable to a factor and set the levels to match your labels
  mutate(Variable = factor(Variable, levels = c('assort_sigma', 'n', 'zi'),
                           labels = c(expression(sigma[m]),
                                      expression(N[p]),
                                      expression(bar(z)))))%>%
  ggplot(aes(x = generations, y = value, col = patch)) +
  geom_path(size = 2) +
  scale_color_manual(values = col2) +
  facet_wrap(~Variable, scales = 'free', 
             labeller = as_labeller(variable_labels, label_parsed)) + 
  theme_bw() +
  theme(strip.text = element_text(size = 15))


#BDMIs
sims%>%
  tidyr::pivot_longer(cols = paste0('bdmi', 1:5), names_to = 'bdmi', values_to = 'allele')%>%
  #filter(generations == 3)%>%
  group_by(generations, patch, bdmi, allele, .drop = FALSE)%>%
  summarise(n = n())%>%
  ungroup%>%
  tidyr::complete(generations, patch, bdmi, allele,
                  fill = list(n = 0))%>% #this adds 0 to groups with 0 counts
  mutate(allele = as.factor(allele))%>%
  mutate(patch = plyr::mapvalues(patch, from = c(1,2),to = c("Patch: 1", "Patch: 2")))%>%
  ggplot(aes(x = bdmi, y = n, group = allele))+
  geom_bar(position="dodge", stat="identity", aes(fill = allele))+
  scale_fill_manual(values = col2)+
  labs(x = "BDMI locus", y = "N")+
  facet_grid(generations~patch)+
  theme_bw()+
  theme(strip.text = element_text(size = 13))
  

#BDMI in multidimensional space
# library(vegan)
# bdmi_data <- sims%>%
#   select(patch, generations, paste0('bdmi',1:5))
# 
# PCA <- bdmi_data%>%
#   select(-patch, -generations)%>%
#   rda()%>%
#   scores()
# 
# PCA$sites%>%
#   bind_cols(bdmi_data)%>%
#   ggplot(aes(x = PC1, y = PC2))+
#   geom_point(size = 5, alpha = 0.5, aes(col = generations))+
#   facet_wrap(~patch)
#   geom_jitter(width = 0.01, height = 0.01)





#Individual-level ---------------
sims%>%
  ggplot(aes(x = zi))+
  geom_density(aes(fill = patch), alpha = 0.5)+
  scale_fill_manual(values = col2)+
  facet_grid(gen~.)+
  theme_bw()
  






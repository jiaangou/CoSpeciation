library(dplyr)
library(ggplot2)


#Group BDMI loci in groups of size b
bdmi_groups <- function(genotype, b = 2){
  n_loci <- length(genotype)
  #Check number of loci is divisible by group_size
  if (n_loci %% b != 0) {
    stop("Number of loci is not divisible by group_size")
  }

  #Combine pairs into a single element
  out <- sapply(seq(from = 1, to = n_loci, by = b), function(i)
    paste0(genotype[i:(i + (b-1))], collapse = ''))

  return(out)
}

#Count the number of bdmi_groups that are of specific type (default is counting the number of heterozygote)
bdmi_types_n <- function(b_groups, type = c('10','01')){
  out <- sum(b_groups %in% type)
  return(out)
}

#Example use-------------
#Example parameters
n_loci <- 10
inds <- 10
group_size <- 2
bdmi_loci <- matrix(rbinom(n_loci*inds,1, prob = 0.9), ncol = n_loci)

#Apply to each individual's genome
apply(bdmi_loci, 1, function(x)x%>%
        bdmi_groups()%>% #groups bdmi loci
        bdmi_types_n()) #count the number of heteros


#Illustration-------------
data.frame(pos = 1:n_loci,
           genotype = rbinom(n_loci, 1, 0.8))%>%
  ggplot(aes(x = pos, y = 1))+
  scale_x_continuous(breaks = 1:10)+
  annotate('line', x = c(0.5, 10.5), y = 0.9, size = 3)+
  annotate('line', x = c(0.5, 10.5), y = 1.1, size = 3)+
  annotate('segment',
           x = seq(from = 0, to = n_loci, by = group_size)+0.5,
           xend = seq(from = 0, to = n_loci, by = group_size)+0.5,
           y = 0.9, yend = 1.1,
           size = 2)+
  annotate('text', label = paste0("Group: ", 1:(n_loci/group_size)),
           x = seq(from = 1, to = 10, by = 2), y = 1.15)+
  geom_text(aes(label = genotype), size = 10)+
  lims(y = c(0.7, 1.3))+labs(x = "BDMI locus")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),
        axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank())


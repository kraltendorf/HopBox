# hopbox - optimization of phenotypic data collection
# kayla altendorf
# monday, february 20, 2023

# load packages
library(dplyr)
library(data.table)
library(ggplot2)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(segmented)


# set working directory
setwd("/Users/kayla.altendorf/OneDrive - USDA/Documents/2023/HopBox/Output/")

# read in data
#dat <- read.csv("./Output/hopbox_data.csv")

#### run t-tests between subsamples and full samples #### 

# create a new column for genotype rep
#validation <- dat %>% mutate(genotype_rep = paste(genotype, rep, sep = "_")) %>%
#  rename(density = mean_cone_density, 
#         weight = avg_cone_weight_adj)

# create a vector of traits, sample sizes, iterations, and genotype_reps for iteration purposes
traits <- c("area", "density", "width", "length", "perimeter", "openness", "weight", "dgci")
sample_size <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100) 
iterations <- 500
#genotype_rep <- unique(validation$genotype_rep)

# create lists to record output from the loop
# trait_list <- list()
# genotype_rep_list <- list()
# 
# for (i in 1:length(genotype_rep)) {
#   for (j in 1:length(traits)) {
#     
#     # extract just the genotype_rep that we're evaluating, and the trait
#     full_sample_genotype_rep <- validation[validation$genotype_rep == genotype_rep[i],][, traits[j]]
#     
#     # calculate the mean for the full sample
#     full_sample_mu <- mean(full_sample_genotype_rep, na.rm = T)
#     
#     # create a new dataframe for the results
#     output <- data.frame(genotype = genotype_rep[i], trait = traits[j], sample_size = rep(sample_size, each = iterations), rep = 1:iterations)
#     
#     # iterate through the output dataframe utilizing the sample size and iteration information 
#     for (k in 1:nrow(output)) {
#       sub_sample <- sample(full_sample_genotype_rep, size = output$sample_size[k], replace = F)
#       output$p_value[k] <- t.test(sub_sample, full_sample_genotype_rep, var.equal = TRUE)[3][[1]]
#       output$full_sample_mu[k] <- full_sample_mu
#       output$sub_sample_mean[k] <- mean(sub_sample)
#       output$coeff_var[k] <- sd(sub_sample) / mean(sub_sample)
#       }
#     
#     trait_list[[j]] <- output
#   }
#   genotype_rep_list[[i]] <-  do.call("rbind", trait_list)
#   print(paste("done with genotype ", genotype_rep[i], " at ", Sys.time(), sep = "")) # provide updates to the terminal when a genotype_rep is complete
# }
# 
# output_all <- do.call("rbind", genotype_rep_list)
# write.csv(output_all, "./Output/output_all_02.20.23.csv")
# 
output_all <- read.csv("./output_all_02.20.23.csv")

#### calculate frequency of non-representative sample #### 
significantly_different <- output_all %>% 
  filter(p_value < 0.05) %>% 
  group_by(genotype, trait, sample_size) %>% 
  tally() %>%
  mutate(merge = paste(genotype, trait, sample_size, sep = "_")) %>%
  ungroup() %>%
  dplyr::select(merge, n)

backbone <- output_all %>% 
  mutate(merge = paste(genotype, trait, sample_size, sep = "_")) %>% 
  dplyr::select(merge) %>%
  distinct()

sig <- left_join(backbone, significantly_different, by = "merge")
sig$n[is.na(sig$n)] <- 0 

sig <- sig %>% mutate(n_total = iterations,
               frequency = n / n_total) %>%
  separate(merge, into = c("genotype", "rep", "trait", "sample_size"), sep = "_") 

sig_sum <- sig %>% group_by(trait, sample_size) %>% summarise(mean = mean(frequency)) %>% as.data.frame() %>% mutate(sample_size = as.numeric(sample_size))

# Segmented regression
t.test.matrix<- as_tibble(model.matrix(~0+sig$trait)*as.numeric(sig$sample_size)) %>%
  mutate(trait = sig$trait,
         frequency = sig$frequency) %>%
  rename("area" = `sig$traitarea`,
         "density" = `sig$traitdensity`,
         "dgci" = `sig$traitdgci`,
         "length" = `sig$traitlength`,
         "openness" = `sig$traitopenness`,
         "perimeter" = `sig$traitperimeter`,
         "weight" = `sig$traitweight`,
         "width" = `sig$traitwidth`)

## making a matrix that can be used in segmented for all 8 traits at once
sig.glm<- glm(frequency ~ 0 + trait + area + density + dgci + length + openness + perimeter + weight + width, data = t.test.matrix)
davies.test(sig.glm, ~ area)
davies.test(sig.glm, ~ density)
davies.test(sig.glm, ~ dgci) # significant
davies.test(sig.glm, ~ length) # significant
davies.test(sig.glm, ~ openness)
davies.test(sig.glm, ~ perimeter) # significant
davies.test(sig.glm, ~ weight)
davies.test(sig.glm, ~ width)

sig.seg<-segmented(sig.glm, 
               seg.Z =~ dgci + length +perimeter,
               psi=list(dgci=c(50), 
                        length=c(50), 
                        perimeter=c(50) 
                        ))
slope(sig.seg)
## None of the slope estimates change much ans so I am assuming a linear relationship

sig.lm1<- lm(frequency ~ as.numeric(sample_size) * trait, data = sig)
anova(sig.lm1) # trait does not inlfuence slope and does not interact with sample size

sig.lm2<- lm(frequency ~ as.numeric(sample_size), data = sig)
summary(sig.lm2) 

# plot results
ggplot(sig, aes(x = as.numeric(sample_size), y = frequency, color = trait))  + 
  geom_point(size = 0.1, position = "jitter") +  
  theme_bw() + 
  stat_smooth(method = "lm" , se = FALSE) +
  scale_color_manual('Trait', 
                     labels = c("Area", "Density", "DGCI", "Length", "Openness", "Perimeter", "Weight", "Width"),
                     values = c("#009E73", "#CC79A7", "#E69F00", "#D55E00", "firebrick4", "#56B4E9", "#0072B2", "grey30")) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size=25, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.position = "bottom") + 
  labs(y = "Proportion of Significant T-Tests", x = "Subample Size") 


















#### assess spearman rank correlations between subsample and full sample ####

# first calculate the full sample means
# full_sample_means <- validation %>% group_by(genotype_rep) %>% summarise(area = mean(area, na.rm = T), 
#                                                     density = mean(density, na.rm = T), 
#                                                     width = mean(width, na.rm = T), 
#                                                     length = mean(length, na.rm = T), 
#                                                     perimeter = mean(perimeter, na.rm = T), 
#                                                     openness = mean(openness, na.rm = T), 
#                                                     weight = mean(weight, na.rm = T), 
#                                                     dgci = mean(dgci, na.rm = T)) %>%
#   rename(genotype = genotype_rep)

# create a list for output, with info on traits and sample sizes to utilize for subsetting
output <- list()

for (j in 1:length(traits)) {
  output[[j]] <- data.frame(trait = traits[j], sample_size = rep(sample_size, each = iterations), rep = 1:iterations)
}

output <- do.call("rbind", output)

# utilizing the output_all from the previous step, extract out the subsample means 
for(i in 1:nrow(output)) {
  t <- output_all[output_all$trait == output$trait[i],]
  t_s <- t[t$sample_size == output$sample_size[i],]
  t_s_i <- t_s[t_s$rep == output$rep[i],]
      
  # evaluate correlation between full sample and the subsample
  result <- cor.test(t_s_i$full_sample_mu, t_s_i$sub_sample_mean, method = "spearman")
  
  output$rho[i] <- result$estimate[[1]]
  output$p_value[i] <- result$p.value
  
  # rank changes? 
  test_full <- t_s_i %>% dplyr::select(genotype, full_sample_mu) %>% 
    separate(genotype, into = c("genotype", "rep"), sep = "_") %>% 
    group_by(genotype) %>%
    summarise(mean = mean(full_sample_mu)) %>%
    arrange(- mean) 
  
  test_subsample <- t_s_i %>% dplyr::select(genotype, sub_sample_mean) %>% 
    separate(genotype, into = c("genotype", "rep"), sep = "_") %>% 
    group_by(genotype) %>%
    summarise(mean = mean(sub_sample_mean)) %>%
    arrange(-mean)
  
  # how many times to the genotypes in the top 3 actually change? 
  output$n_intersect[i] <- length(intersect(test_full$genotype[1:3], test_subsample$genotype[1:3]))
}

write.csv(output, "./Output/subsampling_spearman_02.20.23.csv", row.names = F)
subsampling_spearman <- read.csv("./subsampling_spearman_02.20.23.csv")

## making a matrix that can be used in segmented for all 8 traits at once
spearman.matrix<- as_tibble(model.matrix(~0+subsampling_spearman$trait)*as.numeric(subsampling_spearman$sample_size)) %>%
  mutate(trait = subsampling_spearman$trait,
         rho = subsampling_spearman$rho) %>%
  rename("area" = `subsampling_spearman$traitarea`,
         "density" = `subsampling_spearman$traitdensity`,
         "dgci" = `subsampling_spearman$traitdgci`,
         "length" = `subsampling_spearman$traitlength`,
         "openness" = `subsampling_spearman$traitopenness`,
         "perimeter" = `subsampling_spearman$traitperimeter`,
         "weight" = `subsampling_spearman$traitweight`,
         "width" = `subsampling_spearman$traitwidth`)

## glm to feed the segmentation algorithm
spearman.glm<- glm(rho ~ 0 + trait + area + density + dgci + length + openness + perimeter + weight + width, data = spearman.matrix)

## checking to see if a non linear slope is appropriate
davies.test(spearman.glm, ~ area)# significant
davies.test(spearman.glm, ~ density)# significant
davies.test(spearman.glm, ~ dgci) # significant
davies.test(spearman.glm, ~ length) # significant
davies.test(spearman.glm, ~ openness)# significant
davies.test(spearman.glm, ~ perimeter) # significant
davies.test(spearman.glm, ~ weight)# significant
davies.test(spearman.glm, ~ width)# significant

## calculating breakpoints from the glm 
spearman.seg<-segmented(spearman.glm, 
                   seg.Z =~ area + density + dgci + length + openness + perimeter + weight + width,
                   psi=list(area=c(15,40),
                            density=c(15,40),
                            dgci=c(15, 40), 
                            length=c(15, 40),
                            openness=c(15,40),
                            perimeter=c(15, 40),
                            weight=c(15,40),
                            width=c(15, 40)
                   ))



## Slope and breakpoint coefficients
spearman.psi<- spearman.seg$psi 
spearman.slope<- slope(spearman.seg)
write.csv(spearman.psi, "/Users/garett/Desktop/spearman.psi.csv")
write.csv(spearman.slope, "/Users/garett/Desktop/spearman.slope.csv")




## Fitting a regression line using the segmented model
spearman.newdat<- data.frame(expand_grid(trait = c(unique(spearman.matrix$trait)),
                             sample_size = seq(0,100, length.out=10000)))

spearman.newdatMatrix<- as_tibble(model.matrix(~0+spearman.newdat$trait)*as.numeric(spearman.newdat$sample_size)) %>%
  mutate(trait = spearman.newdat$trait) %>%
  rename("area" = `spearman.newdat$traitarea`,
         "density" = `spearman.newdat$traitdensity`,
         "dgci" = `spearman.newdat$traitdgci`,
         "length" = `spearman.newdat$traitlength`,
         "openness" = `spearman.newdat$traitopenness`,
         "perimeter" = `spearman.newdat$traitperimeter`,
         "weight" = `spearman.newdat$traitweight`,
         "width" = `spearman.newdat$traitwidth`)


spearman.fitted<- predict.segmented(spearman.seg, newdata = spearman.newdatMatrix)
spearman.prediction<- data.frame(sample_size = spearman.newdat$sample_size, 
                                  trait = spearman.newdat$trait,
                                  rho = spearman.fitted)



# create a figure showing spearman rank correlation coefficients
ggplot(subsampling_spearman, aes(x = sample_size, y = rho, color = trait))  + 
  geom_point(size = 0.01, 
             position = "jitter",
             alpha=0.5) + 
  theme_bw() +
  geom_line(data=spearman.prediction, 
            mapping=aes(x = sample_size, y = rho, color = trait),
            size=1) +
  scale_color_manual('Trait',
                     labels = c("Area", "Density", "DGCI", "Length", "Openness", "Perimeter", "Weight", "Width"),
                     values = c("#009E73", "#CC79A7", "#E69F00", "#D55E00", "firebrick4", "#56B4E9", "#0072B2", "grey30")) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size=25, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.position = "bottom") + 
  labs(y = "Spearman Rank Coefficient", x = "Subsample Size")

















subsampling_spearman %>% filter(p_value > 0.05) # there are none that are insignificant
head(subsampling_spearman)

zero_wrong <- subsampling_spearman %>% 
  filter(n_intersect == 3) %>% 
  group_by(trait, sample_size) %>%
  tally() %>%
  mutate(frequency = 1 - (n/500))


one_wrong <- subsampling_spearman %>% 
  filter(n_intersect == 2) %>% 
  group_by(trait, sample_size) %>%
  tally() %>%
  mutate(frequency = n / 500)

two_wrong <- subsampling_spearman %>%
  filter(n_intersect == 1) %>%
  group_by(trait, sample_size) %>%
  tally() %>%
  mutate(frequency = n / 500)

three_wrong <- subsampling_spearman %>%
  filter(n_intersect == 0) %>%
  group_by(trait, sample_size) %>%
  tally() %>%
  mutate(frequency = n / 500)

## making a matrix that can be used in segmented for all 8 traits at once
zero_wrong.matrix<- as_tibble(model.matrix(~0+zero_wrong$trait)*as.numeric(zero_wrong$sample_size)) %>%
  mutate(trait = zero_wrong$trait,
         frequency = zero_wrong$frequency) %>%
  rename("area" = `zero_wrong$traitarea`,
         "density" = `zero_wrong$traitdensity`,
         "dgci" = `zero_wrong$traitdgci`,
         "length" = `zero_wrong$traitlength`,
         "openness" = `zero_wrong$traitopenness`,
         "perimeter" = `zero_wrong$traitperimeter`,
         "weight" = `zero_wrong$traitweight`,
         "width" = `zero_wrong$traitwidth`)

## glm to feed the segmentation algorithm
zero_wrong.glm<- glm(frequency ~ 0 + trait + area + density + dgci + length + openness + perimeter + weight + width, data = zero_wrong.matrix)

## checking to see if a non linear slope is appropriate
davies.test(zero_wrong.glm, ~ area)# significant
davies.test(zero_wrong.glm, ~ density)# significant
davies.test(zero_wrong.glm, ~ dgci) # significant
davies.test(zero_wrong.glm, ~ length) # meh
davies.test(zero_wrong.glm, ~ openness)# significant
davies.test(zero_wrong.glm, ~ perimeter) # significant
davies.test(zero_wrong.glm, ~ weight)
davies.test(zero_wrong.glm, ~ width)# significant

## calculating breakpoints from the glm 
zero_wrong.seg<-segmented(zero_wrong.glm, 
                        seg.Z =~ area + density + dgci + length + openness + perimeter + weight + width,
                        psi=list(area=c(20),
                                 density=c(35),
                                 dgci=c(40), 
                                 length=c(25),
                                 openness=c(30),
                                 perimeter=c(28),
                                 weight=c(1),
                                 width=c(30)
                        ))

## Slope and breakpoint coefficients
zero_wrong.psi<- zero_wrong.seg$psi 
zero_wrong.slope<- slope(zero_wrong.seg)
#write.csv(zero_wrong.psi, "/Users/garett/Desktop/zero_wrong.psi.csv")
#write.csv(zero_wrong.slope, "/Users/garett/Desktop/zero_wrong.slope.csv")

## Fitting a regression line using the segmented model
zero_wrong.newdat<- data.frame(expand_grid(trait = c(unique(zero_wrong.matrix$trait)),
                                         sample_size = seq(0,100, length.out=10000)))

zero_wrong.newdatMatrix<- as_tibble(model.matrix(~0+zero_wrong.newdat$trait)*as.numeric(zero_wrong.newdat$sample_size)) %>%
  mutate(trait = zero_wrong.newdat$trait) %>%
  rename("area" = `zero_wrong.newdat$traitarea`,
         "density" = `zero_wrong.newdat$traitdensity`,
         "dgci" = `zero_wrong.newdat$traitdgci`,
         "length" = `zero_wrong.newdat$traitlength`,
         "openness" = `zero_wrong.newdat$traitopenness`,
         "perimeter" = `zero_wrong.newdat$traitperimeter`,
         "weight" = `zero_wrong.newdat$traitweight`,
         "width" = `zero_wrong.newdat$traitwidth`)


zero_wrong.fitted<- predict.segmented(zero_wrong.seg, newdata = zero_wrong.newdatMatrix)
zero_wrong.prediction<- data.frame(sample_size = zero_wrong.newdat$sample_size, 
                                 trait = zero_wrong.newdat$trait,
                                 frequency = zero_wrong.fitted)



# determine the frequency of rank changes in the top 3 spots, which would be about a 20% selection intensity. 
ggplot(zero_wrong, aes(x = sample_size, y = frequency, color = trait))  + 
  geom_point(size = 1, 
             alpha=0.5) + 
  theme_bw() +
  geom_line(data=zero_wrong.prediction, 
            mapping=aes(x = sample_size, 
                        y = frequency, color = trait),
            size=1) +
  scale_color_manual('Trait',
                     labels = c("Area", "Density", "DGCI", "Length", "Openness", "Perimeter", "Weight", "Width"),
                     values = c("#009E73", "#CC79A7", "#E69F00", "#D55E00", "firebrick4", "#56B4E9", "#0072B2", "grey30")) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size=25, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.position = "bottom") + 
  labs(y = expression( bold("Probability ">="1 Erroneous Selections", )), x = "Subsample Size") 



















#### assess the threshold of f-value for determining significant differences among genotypes ####
# create a new list for output
output <- list()
for (j in 1:length(traits)) {
  output[[j]] <- data.frame(trait = traits[j], sample_size = rep(sample_size, each = iterations), rep = 1:iterations)
}
output <- do.call("rbind", output)

# iterate through subsamples, run a model and view f-value from the anova
for (i in 1:nrow(output)) {
  validation_trait <- validation[, c("genotype_rep", "genotype", "rep", output$trait[i])] %>%
    group_by(genotype_rep) %>% slice_sample(n = output$sample_size[i], replace = F)
  
  model <- paste(output$trait[i], " ~ genotype + (1|rep) + (1|rep:genotype)", sep = "")
  output$f_val[i] <- suppressMessages(anova(lmer(model, data = validation_trait))[[5]])
}

write.csv(output, "./Output/subsampling_f_values.csv", row.names = F)
f_val <- read.csv("./subsampling_f_values.csv")

## making a matrix that can be used in segmented for all 8 traits at once
f_val.matrix<- as_tibble(model.matrix(~0+f_val$trait)*as.numeric(f_val$sample_size)) %>%
  mutate(trait = f_val$trait,
         f_val = f_val$f_val) %>%
  rename("area" = `f_val$traitarea`,
         "density" = `f_val$traitdensity`,
         "dgci" = `f_val$traitdgci`,
         "length" = `f_val$traitlength`,
         "openness" = `f_val$traitopenness`,
         "perimeter" = `f_val$traitperimeter`,
         "weight" = `f_val$traitweight`,
         "width" = `f_val$traitwidth`)


## glm to feed the segmentation algorithm
f_val.glm<- glm(f_val ~ 0 + trait + area + density + dgci + length + openness + perimeter + weight + width, data = f_val.matrix)


## checking to see if a non linear slope is appropriate
davies.test(f_val.glm, ~ area)# significant
davies.test(f_val.glm, ~ density)# significant
davies.test(f_val.glm, ~ dgci) 
davies.test(f_val.glm, ~ length) # significant
davies.test(f_val.glm, ~ openness)# significant
davies.test(f_val.glm, ~ perimeter) # significant
davies.test(f_val.glm, ~ weight)# significant
davies.test(f_val.glm, ~ width)# significant


## calculating breakpoints from the glm 
f_val.seg<-segmented(f_val.glm, 
                        seg.Z =~ area + density + dgci + length + openness + perimeter + weight + width,
                        psi=list(area=c(15),
                                 density=c(5),
                                 dgci=c(5), 
                                 length=c(20),
                                 openness=c(20),
                                 perimeter=c(30),
                                 weight=c(15),
                                 width=c(15)
                        ))


## Slope and breakpoint coefficients
f_val.psi<- f_val.seg$psi 
f_val.slope<- slope(f_val.seg)
#write.csv(f_val.psi, "/Users/garett/Desktop/f_val.psi.csv")
#write.csv(f_val.slope, "/Users/garett/Desktop/f_val.slope.csv")


## Fitting a regression line using the segmented model
f_val.newdat<- data.frame(expand_grid(trait = c(unique(f_val.matrix$trait)),
                                           sample_size = seq(0,100, length.out=2000)))

f_val.newdatMatrix<- as_tibble(model.matrix(~0+f_val.newdat$trait)*as.numeric(f_val.newdat$sample_size)) %>%
  mutate(trait = f_val.newdat$trait) %>%
  rename("area" = `f_val.newdat$traitarea`,
         "density" = `f_val.newdat$traitdensity`,
         "dgci" = `f_val.newdat$traitdgci`,
         "length" = `f_val.newdat$traitlength`,
         "openness" = `f_val.newdat$traitopenness`,
         "perimeter" = `f_val.newdat$traitperimeter`,
         "weight" = `f_val.newdat$traitweight`,
         "width" = `f_val.newdat$traitwidth`)


f_val.fitted<- predict.segmented(f_val.seg, newdata = f_val.newdatMatrix)
f_val.prediction<- data.frame(sample_size = f_val.newdat$sample_size, 
                                   trait = f_val.newdat$trait,
                              f_val = f_val.fitted)


# graph of f-values
ggplot(f_val, aes(x = as.numeric(sample_size), y = f_val, color = trait))  + 
  geom_point(size = 0.01, 
             position = "jitter",
             alpha=0.5) + 
  theme_classic() + 
  geom_line(data=f_val.prediction, 
            mapping=aes(x = sample_size, y = f_val, color = trait),
            size=1) +
  geom_hline(yintercept = 2.48, linetype = "dotted") +
  geom_hline(yintercept = 5.93, linetype = "dashed") +
  scale_color_manual('Trait',
                     labels = c("Area", "Density", "DGCI", "Length", "Openness", "Perimeter", "Weight", "Width"),
                     values = c("#009E73", "#CC79A7", "#E69F00", "#D55E00", "firebrick4", "#56B4E9", "#0072B2", "grey30")) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size=25, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.position = "bottom") + 
  ylim(0,50)+
  labs(y = "F-Value", x = "Subsample Size") 



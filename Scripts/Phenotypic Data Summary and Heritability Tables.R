# hopbox - summary of phenotypic data
# kayla altendorf
# monday, january 9 2023

# load packages

# set directory
setwd("/users/kayla.altendorf/OneDrive - USDA/Documents/2023/HopBox/")

# read in data
files <-  list.files("./Output", pattern = c("emmean"), full.names = T)

# set output
output <- list()

for (i in 1:length(files)) {
  dat <- read.csv(files[i]) %>% dplyr::select(1:4) 
  colnames(dat)[4] <- "value"
  dat$trait <- gsub(".csv", "", str_split(files[i],pattern = "_")[[1]][2])
  output[[i]] <- dat
}

output <- do.call("rbind", output) %>% 
  pivot_wider(names_from = trait, values_from = value) %>% 
  rename(weight = avg, 
         cone_density = cone)

colnames(output)

summary <- output %>% group_by(genotype) %>% summarise(mean_length = mean(length, na.rm = T), 
                                            min_length = min(length, na.rm = T), 
                                            max_length = max(length, na.rm = T), 
                                            cv_length = (sd(length, na.rm = T) / mean(length, na.rm = T)),
                                            mean_width = mean(width, na.rm = T), 
                                            min_width = min(width, na.rm = T), 
                                            max_width = max(width, na.rm = T), 
                                            cv_width = (sd(width, na.rm = T) / mean(width, na.rm = T)),
                                            mean_area = mean(area, na.rm = T), 
                                            min_area = min(area, na.rm = T), 
                                            max_area = max(area, na.rm = T), 
                                            cv_area = (sd(area, na.rm = T) / mean(area, na.rm = T)),
                                            mean_perimeter = mean(perimeter, na.rm = T), 
                                            min_perimeter = min(perimeter, na.rm = T), 
                                            max_perimeter = max(perimeter, na.rm = T), 
                                            cv_perimeter = (sd(perimeter, na.rm = T) / mean(perimeter, na.rm = T)),
                                            mean_openness = mean(openness, na.rm = T), 
                                            min_openness = min(openness, na.rm = T), 
                                            max_openness = max(openness, na.rm = T), 
                                            cv_openness = (sd(openness, na.rm = T) / mean(openness, na.rm = T)),
                                            mean_weight = mean(weight, na.rm = T), 
                                            min_weight = min(weight, na.rm = T), 
                                            max_weight = max(weight, na.rm = T), 
                                            cv_weight = (sd(weight, na.rm = T) / mean(weight, na.rm = T)),
                                            mean_dgci = mean(dgci, na.rm = T), 
                                            min_dgci = min(dgci, na.rm = T), 
                                            max_dgci = max(dgci, na.rm = T), 
                                            cv_dgci = (sd(dgci, na.rm = T) / mean(dgci, na.rm = T)),
                                            mean_cone_density = mean(cone_density, na.rm = T), 
                                            min_cone_density = min(cone_density, na.rm = T), 
                                            max_cone_density = max(cone_density, na.rm = T), 
                                            cv_cone_density = (sd(cone_density, na.rm = T) / mean(cone_density, na.rm = T)))
                                            
                                            
write.csv(summary, "./Output/Formatted Tables and Figures/phentoypic_data_summary.csv", row.names = F)

heritability <- data.frame(trait = c("length", "width", "area", "perimeter", "openness", "avg_cone_weight_adj", "dgci", "cone_density"), 
                           heritability = c( 0.2309142, 0.50834, 0.3609637, 0.4026552, 0.590153, 0.4153846, 0.2899855, -0.5164697))

write.csv(heritability, "./Output/Formatted Tables and Figures/heritability.csv", row.names = F)


traits <- c("length", "width", "area", "perimeter", "openness", "avg_cone_weight_adj", "dgci", "cone_density")

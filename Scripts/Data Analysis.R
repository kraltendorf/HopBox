# hopbox data analysis
# kayla altendorf

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
library(cowplot)

# set working directory
setwd("/users/kayla.altendorf/OneDrive - USDA/Documents/2023/HopBox/")

# read in hopbox data and filter out images to just include HB images
files <- list.files("./Data", full.names = TRUE, pattern = "results")
dat <- bind_rows(lapply(files, fread)) %>% filter(grepl("HB", QR_name))

# set pixels per cm
pixels_per_cm <- 597.9674/5

# format data - extract genotype name
dat <- dat %>% separate(File_name, sep =  "_", into = c("HB", "genotype", "rep", "image")) %>%
  dplyr::select(-V1, -HB) %>% 
  mutate(area = Area / (pixels_per_cm^2), 
         perimeter = Perimeter / pixels_per_cm, 
         length = Major_axis / pixels_per_cm, 
         width = Minor_axis / pixels_per_cm, 
         openness = perimeter / length,
         QR_name = str_remove(QR_name, ".JPG"), 
         image = str_remove(image, ".JPG")) %>%
  dplyr::select(-Area, -Perimeter, -Major_axis, -Minor_axis)

# switch B and R for color output - this issue was resolved in an updated version of the image analysis pipeline
colnames(dat)[11] <- "Median_B"
colnames(dat)[13] <- "Median_R"
colnames(dat)[14] <- "std_B"
colnames(dat)[16] <- "std_R"
colnames(dat)[17] <- "Mean_B"
colnames(dat)[19] <- "Mean_R"

# read in and format image weight data
weights <- read_excel("./Data/HB_Image_Weights_092221.xlsx") %>%
  separate(Genotype, sep = "_", into = c("HB", "genotype", "rep")) %>%
  mutate(QR_name = paste("HB", genotype, rep, Image, sep = "_"),
         genotype_rep = paste(genotype, rep, sep = "_")) %>%
  rename(image = Image, 
         weight_g = Weight)

# read in and format dry matter data
dry_matter <- read_excel("./Data/HB_Dry_Matter_092321.xlsx") %>%
  mutate(genotype_rep = paste(genotype, rep, sep = "_")) %>%
  dplyr::select(genotype_rep, `Dry Matter`) %>%
  rename(dry_matter = `Dry Matter`)

# calculate cone weight per image adjusted for dry matter
weights_dry_matter <- left_join(weights, dry_matter, by = "genotype_rep") %>% 
  mutate(weight_g_adj = (weight_g * dry_matter)) %>% 
  dplyr::select(QR_name, weight_g, dry_matter, weight_g_adj)

# determine the number of cones per image
cone_n <- dat %>% group_by(QR_name) %>% tally() %>% rename(cone_n = n)

# calculate average cone weight, adjusted for dry matter
avg_cone_weight <- left_join(weights_dry_matter, cone_n, by = "QR_name") %>% 
  mutate(avg_cone_weight_adj = weight_g_adj / cone_n)

# bind image weight info to data
dat2 <- left_join(dat, cone_n, by = "QR_name")
dat2 <- left_join(dat2, avg_cone_weight, by = "QR_name")

# create a data frame of means for each trait to calculate heritability
dat3 <- dat2 %>% group_by(genotype, rep, image) %>% summarise(mean_length = mean(length, na.rm = T), 
                                                              mean_area = mean(area, na.rm = T), 
                                                              mean_perimeter = mean(perimeter, na.rm = T), 
                                                              mean_openness = mean(openness, na.rm = T), 
                                                              mean_width = mean(width, na.rm = T),
                                                              mean_weight = mean(avg_cone_weight_adj, na.rm = T)) %>%
  mutate(mean_cone_density = (mean_weight / mean_area))

# linear models for each of the traits

#### area ####
# linear model
lm_area <- lmer(sqrt(area) ~ rep + genotype*image + (1|genotype:rep) + (1|genotype:rep:image), data = dat2)
anova_area <- anova(lm_area, type = "II")
write.csv(anova_area, "./Output/anova_area.csv", row.names = T)

# heritability
lm_area_random <- lmer(sqrt(mean_area) ~ (1|genotype) + (1|rep), data = dat3)
summary(lm_area_random)

vg <- (0.081289 - 0.038169) / 2
ve <- 0.038169

H_area <- vg / (vg + ve)
H_area

# emmeans
emmean_area <- as.data.frame(emmeans(lm_area, ~ genotype + rep + image, type = "response")) %>%
  rename(area = response)
write.csv(emmean_area, "./Output/emmean_area.csv", row.names = F)
emmean_area <- read.csv("./Output/emmean_area.csv")


emmean_genotype_area <- emmeans(lm_area, ~ genotype, type = "response")
cld_genotype_area <- as.data.frame(cld(emmean_genotype_area, Letters=letters)) %>%
  rename(area = response)
write.csv(emmean_genotype_area, "./Output/emmean_genotype_area.csv", row.names = F)
emmean_genotype_area <- read.csv("./Output/emmean_genotype_area.csv")
cld_genotype_area <- cld_genotype_area %>% mutate(.group = gsub(" ", "", .group))


# ggplot2
gg_area <- ggplot() +
  geom_boxplot(emmean_area, mapping = aes(x = reorder(genotype, area), y = area), color = "#000000", size = 0.3)  +
  geom_text(cld_genotype_area, mapping = aes(x = reorder(genotype, area),
                                             y = area + 2.5,
                                             label = .group), size = 2.5) +
  theme_bw() + 
  labs(y = bquote(~bold('Cone Area ('*cm^2*')')), x = " ") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10),
        axis.title = element_text(size=13, face = "bold")) 
  #annotate("text", label = "H = 0.65", x = 14.5, y = 1, color = "black", size = 8)

ggsave(filename = "/users/kayla.altendorf/OneDrive - USDA/Documents/2023/HopBox/Output/Formatted Tables and Figures/Area_Boxplot.jpg", 
       width = 5, height = 4, units = "in")


#### perimeter ####
# linear model
lm_perimeter <- lmer(sqrt(perimeter) ~ rep + genotype*image + (1|genotype:rep) + (1|genotype:rep:image), data = dat2)
anova_perimeter <- anova(lm_perimeter, type = "II")
write.csv(anova_perimeter, "./Output/anova_perimeter.csv", row.names = T)

# heritability
lm_perimeter_random  <- lmer(sqrt(mean_perimeter) ~ (1|genotype) + (1|rep), data = dat3)
summary(lm_perimeter_random)

vg <- (0.084524 - 0.035996) / 2
ve <- 0.035996

H_perimeter <- vg / (vg + ve)
H_perimeter

# emmeans
emmean_perimeter <- as.data.frame(emmeans(lm_perimeter, ~ genotype + rep + image, type = "response")) %>%
  rename(perimeter = response)
write.csv(emmean_perimeter, "./Output/emmean_perimeter.csv", row.names = F)


emmean_genotype_perimeter <- emmeans(lm_perimeter, ~ genotype, type = "response")
cld_genotype_perimeter <- as.data.frame(cld(emmean_genotype_perimeter, Letters=letters)) %>%
  rename(perimeter = response)


# ggplot2
ggplot() +
  geom_boxplot(emmean_perimeter, mapping = aes(x = reorder(genotype, perimeter), y = perimeter), color = "#00868B", size = 1)  +
  geom_text(cld_genotype_perimeter, mapping = aes(x = reorder(genotype, perimeter),
                                             y = perimeter + 2.5,
                                             label = .group), size = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25),
        axis.title = element_text(size=25, face = "bold")) + 
  labs(y = bquote("Cone Perimeter (cm)"), x = "Genotype")  + 
  annotate("text", label = "H = 0.59", x = 14.5, y = 7, color = "black", size = 8)




#### length #####
# linear model
lm_length <- lmer(sqrt(length) ~ rep + genotype*image + (1|genotype:rep) + (1|genotype:rep:image), data = dat2)
anova_length <- anova(lm_length, type = "II")
write.csv(anova_length, "./Output/anova_length.csv", row.names = T)

# heritability
lm_length_random  <- lmer(sqrt(mean_length) ~ (1|genotype) + (1|rep), data = dat3)
summary(lm_length_random)

vg <- (0.0164580 - 0.0102831) / 2
ve <- 0.0102831

H_length <- vg / (vg + ve)
H_length

# emmeans
emmean_length <- as.data.frame(emmeans(lm_length, ~ genotype + rep + image, type = "response")) %>%
  rename(length = response)
write.csv(emmean_length, "./Output/emmean_length.csv", row.names = F)


emmean_genotype_length <- emmeans(lm_length, ~ genotype, type = "response")
cld_genotype_length <- as.data.frame(cld(emmean_genotype_length, Letters=letters)) %>%
  rename(length = response)


# ggplot2
ggplot() +
  geom_boxplot(emmean_length, mapping = aes(x = reorder(genotype, length), y = length), color = "#00868B", size = 1)  +
  geom_text(cld_genotype_length, mapping = aes(x = reorder(genotype, length),
                                                  y = length + 1,
                                                  label = .group), size = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25),
        axis.title = element_text(size=25, face = "bold")) + 
  labs(y = bquote("Cone Length (cm)"), x = "Genotype")  + 
  annotate("text", label = "H = 0.59", x = 14.5, y = 2, color = "black", size = 8)


#### width #####
# linear model
lm_width <- lmer(sqrt(width) ~ rep + genotype*image + (1|genotype:rep) + (1|genotype:rep:image), data = dat2)
anova_width <- anova(lm_width, type = "II")
write.csv(anova_width, "./Output/anova_width.csv", row.names = T)

# heritability
lm_width_random  <- lmer(sqrt(mean_width) ~ (1|genotype) + (1|rep), data = dat3)
summary(lm_width_random)

vg <- (0.0078356 - 0.0025541) / 2
ve <- 0.0025541

H_width <- vg / (vg + ve)
H_width

# emmeans
emmean_width <- as.data.frame(emmeans(lm_width, ~ genotype + rep + image, type = "response")) %>%
  rename(width = response)
write.csv(emmean_width, "./Output/emmean_width.csv", row.names = F)


emmean_genotype_width <- emmeans(lm_width, ~ genotype, type = "response")
cld_genotype_width <- as.data.frame(cld(emmean_genotype_width, Letters=letters)) %>%
  rename(width = response)


# ggplot2
ggplot() +
  geom_boxplot(emmean_width, mapping = aes(x = reorder(genotype, width), y = width), color = "black", size = 1)  +
  geom_text(cld_genotype_width, mapping = aes(x = reorder(genotype, width),
                                               y = width + 0.5,
                                               label = .group), size = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25),
        axis.title = element_text(size=25, face = "bold")) + 
  labs(y = bquote("Cone width (cm)"), x = "Genotype")  + 
  annotate("text", label = "H = 0.59", x = 14.5, y = 1.5, color = "black", size = 8)


#### openness ####
# linear model
lm_openness <- lmer(sqrt(openness) ~ rep + genotype*image + (1|genotype:rep) + (1|genotype:rep:image), data = dat2)
anova_openness <- anova(lm_openness, type = "II")
write.csv(anova_openness, "./Output/anova_openness.csv", row.names = T)

# heritability
lm_openness_random  <- lmer(sqrt(mean_openness) ~ (1|genotype) + (1|rep), data = dat3)
summary(lm_openness_random)

vg <- (0.0048898 - 0.0012603) / 2
ve <- 0.0012603

H_openness <- vg / (vg + ve)
H_openness

# emmeans
emmean_openness <- as.data.frame(emmeans(lm_openness, ~ genotype + rep + image, type = "response")) %>%
  rename(openness = response)
write.csv(emmean_openness, "./Output/emmean_openness.csv", row.names = F)
emmean_openness <- read.csv("./Output/emmean_openness.csv")

emmean_genotype_openness <- emmeans(lm_openness, ~ genotype, type = "response")
cld_genotype_openness <- as.data.frame(cld(emmean_genotype_openness, Letters=letters)) %>%
  rename(openness = response)
write.csv(cld_genotype_openness, "./Output/cld_genotype_openness.csv", row.names = F)
cld_genotype_openness <- read.csv("./Output/cld_genotype_openness.csv")
cld_genotype_openness <- cld_genotype_openness %>% mutate(.group = gsub(" ", "", .group))

# ggplot2
gg_openess <- ggplot() +
  geom_boxplot(emmean_openness, mapping = aes(x = reorder(genotype, openness), y = openness), color = "black", size = 0.3)  +
  geom_text(cld_genotype_openness, mapping = aes(x = reorder(genotype, openness),
                                              y = openness + 0.5,
                                              label = .group), size = 2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")) + 
  labs(y = bquote("Cone Openness"), x = "Genotype")  
  #annotate("text", label = "H = 0.59", x = 14.5, y = 1.5, color = "black", size = 8)


#### avg_cone_weight_adj ####
dat4 <- dat2 %>% group_by(genotype, rep, image) %>% summarise(avg_cone_weight_adj = mean(avg_cone_weight_adj, na.rm = T))
View(dat2)
lm_avg_cone_weight_adj <- lmer(sqrt(avg_cone_weight_adj) ~ rep + genotype*image + (1|genotype:rep), data = dat4)
anova_avg_cone_weight_adj <- anova(lm_avg_cone_weight_adj, type = "II")
write.csv(anova_avg_cone_weight_adj, "./Output/anova_avg_cone_weight_adj.csv", row.names = T)

# heritability
lm_weight_random  <- lmer(sqrt(mean_weight) ~ (1|genotype) + (1|rep), data = dat3)
summary(lm_weight_random)

View(dat3)
vg <- (3.680e-03 - 1.520e-03) / 2
ve <- 1.520e-03

H_weight <- vg / (vg + ve)
H_weight


# emmeans
emmean_avg_cone_weight_adj <- as.data.frame(emmeans(lm_avg_cone_weight_adj, ~ genotype + rep + image, type = "response")) %>%
  rename(avg_cone_weight_adj = response)
write.csv(emmean_avg_cone_weight_adj, "./Output/emmean_avg_cone_weight_adj.csv", row.names = F)


emmean_genotype_avg_cone_weight_adj <- emmeans(lm_avg_cone_weight_adj, ~ genotype, type = "response")
cld_genotype_avg_cone_weight_adj <- as.data.frame(cld(emmean_genotype_avg_cone_weight_adj, Letters=letters)) %>%
  rename(avg_cone_weight_adj = response)

# ggplot2
ggplot() +
  geom_boxplot(emmean_avg_cone_weight_adj, mapping = aes(x = reorder(genotype, avg_cone_weight_adj), y = avg_cone_weight_adj), color = "black", size = 1)  +
  geom_text(cld_genotype_avg_cone_weight_adj, mapping = aes(x = reorder(genotype, avg_cone_weight_adj),
                                                 y = avg_cone_weight_adj + 0.1,
                                                 label = .group), size = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25),
        axis.title = element_text(size=25, face = "bold")) + 
  labs(y = bquote("Cone avg_cone_weight_adj"), x = "Genotype")  + 
  annotate("text", label = "H = 0.59", x = 14.5, y = 0.5, color = "black", size = 8)


#### dark green color index ####
# dark green color index
dat5 <- dat2 %>% dplyr::select(QR_name, Label, Mean_R, Mean_G, Mean_B) 
IDs <- dat5[,1:2]
rgb <- as.matrix(t(dat5[,-1:-2]))

hsv <- rgb2hsv(rgb)
hsv <- as.data.frame(t(hsv))
id_hsv <- cbind(IDs, hsv) %>% 
  mutate(h = h*360, 
         QR_name_Label = paste(QR_name, Label, sep = "_"))
                                     

id_hsv_dgci <- id_hsv %>% mutate(dgci = (((h-60)/60) + (1 - s) + (1 - v))/3)

dat6 <- dat2 %>% mutate(QR_name_Label = paste(QR_name, Label, sep = "_"))
dat8 <- dat3 %>% mutate(QR_name = paste("HB", genotype, rep, image, sep = "_")) %>% ungroup() %>% dplyr::select(QR_name, mean_cone_density)
dat7 <- left_join(dat6, id_hsv_dgci %>% dplyr::select(-QR_name, -Label), by = "QR_name_Label")
dat9 <- left_join(dat7, dat8, by = "QR_name")
write.csv(dat9, "./Output/hopbox_data.csv", row.names = F)

mu <- id_hsv_dgci %>% separate(QR_name, sep = "_", into = c("HB", "genotype", "rep", "image")) %>% group_by(genotype) %>% summarise(mean = mean(dgci)) %>% arrange(-mean)

# linear model
lm_dgci <- lmer(dgci ~ rep + genotype*image + (1|genotype:rep) + (1|genotype:rep:image), data = dat7)
anova_dgci <- anova(lm_dgci, type = "II")
write.csv(anova_dgci, "./Output/anova_dgci.csv", row.names = T)

# heritability
dat8 <- dat7 %>% group_by(genotype, rep, image) %>% summarise(mean_dgci = mean(dgci, na.rm = T))

lm_dgci_random  <- lmer(mean_dgci ~ (1|genotype) + (1|rep), data = dat8)
summary(lm_dgci_random)

vg <- (2.222e-04 - 1.223e-04) / 2
ve <- 1.223e-04

H_dgci <- vg / (vg + ve)
H_dgci

# emmeans
emmean_dgci <- as.data.frame(emmeans(lm_dgci, ~ genotype + rep + image))
write.csv(emmean_dgci, "./Output/emmean_dgci.csv", row.names = F)
emmean_dgci <- read.csv("./Output/emmean_dgci.csv")

emmean_genotype_dgci <- emmeans(lm_dgci, ~ genotype)
cld_genotype_dgci <- as.data.frame(cld(emmean_genotype_dgci, Letters=letters))
write.csv(cld_genotype_dgci, "./Output/cld_genotype_dgci.csv")
cld_genotype_dgci <- read.csv("./Output/cld_genotype_dgci.csv")
cld_genotype_dgci <- cld_genotype_dgci %>% mutate(.group = gsub(" ", "", .group))

View(dat7)

# ggplot2
gg_dgci <- ggplot() +
  geom_boxplot(emmean_dgci, mapping = aes(x = reorder(genotype, emmean), y = emmean), color = "black", size = 0.3)  +
  geom_text(cld_genotype_dgci, mapping = aes(x = reorder(genotype, emmean),
                                                 y = emmean + 0.025,
                                                 label = .group), size = 2.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")) + 
  labs(y = bquote("DGCI"), x = " ")   
  #annotate("text", label = "H = 0.59", x = 14.5, y = 0.5, color = "black", size = 8)


#### cone density ####
lm_cone_density <- lmer(mean_cone_density ~ rep + genotype*image + (1|genotype:rep), data = dat3)
anova_cone_density <- anova(lm_cone_density, type = "II")
write.csv(anova_cone_density, "./Output/anova_anova_cone_density.csv", row.names = T)

# heritability ### is negative and need to figure out why
test_dat <- dat3 %>% dplyr::select(genotype, rep, image, mean_cone_density)
hist(sqrt(dat3$mean_cone_density))
lm_cone_density_random  <- lmer(sqrt(mean_cone_density) ~ (1|genotype) + (1|image), data = as.data.frame(test_dat))
summary(lm_cone_density_random)

vg <- (0.0001534 - 0.0004811) / 2
ve <- 0.0004811

H_cone_density <- vg / (vg + ve)
H_cone_density

# emmeans - this doesn't make sense when there's only one value per image
emmean_cone_density <- as.data.frame(emmeans(lm_cone_density, ~ genotype + image + rep, type = "response", df = Inf)) #%>%
  rename(cone_density = emmean)
write.csv(emmean_cone_density, "./Output/emmean_cone_density.csv", row.names = F)

hist(dat3$mean_cone_density)


emmean_genotype_cone_density <- emmeans(lm_cone_density, ~ genotype, type = "response")
cld_genotype_cone_density <- as.data.frame(cld(emmean_genotype_cone_density, Letters=letters)) #%>%
  #rename(cone_density = response)
write.csv(cld_genotype_cone_density, "./Output/cld_genotype_conedensity.csv")

# ggplot2
ggplot() +
  geom_boxplot(emmean_cone_density, mapping = aes(x = reorder(genotype, cone_density), y = cone_density), color = "black", size = 1)  +
  geom_text(cld_genotype_cone_density, mapping = aes(x = reorder(genotype, avg_cone_weight_adj),
                                                            y = avg_cone_weight_adj + 0.1,
                                                            label = .group), size = 10) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 25, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25),
        axis.title = element_text(size=25, face = "bold")) + 
  labs(y = bquote("Cone avg_cone_weight_adj"), x = "Genotype")  + 
  annotate("text", label = "H = 0.59", x = 14.5, y = 0.5, color = "black", size = 8)

grid <- plot_grid(gg_area, gg_dgci, gg_openess, labels = c("A", "B", "C"), nrow = 3)

ggsave("/users/kayla.altendorf/OneDrive - USDA/Documents/2023/HopBox/Output/Formatted Tables and Figures/boxplot_grid.jpg", 
       width = 6, height = 9, units = "in")

### compile all data for output ###
Heritability <- data.frame(Trait = c("Length", "Width", "Area", "Perimeter", "Openness", "Weight", "DGCI"), 
           Heritability = c(H_length, H_width, H_area, H_perimeter, H_openness, H_weight, H_dgci))


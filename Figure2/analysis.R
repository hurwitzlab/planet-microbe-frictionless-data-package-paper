install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
library("ggplot2")
library("dplyr")
library("tidyverse")

#http://ggplot2.tidyverse.org/reference/facet_grid.html
#http://ggplot2.tidyverse.org/reference/geom_boxplot.html

directory <- getwd()
setwd(directory)
custom_pallet = c("#00AFBB", "#f36213", "#E7B800", "#bb00af", "#b40025", "#002fe7", "#ff0136", "#08efff", "#fbbc6a","#ff08ef", "#a18000","#736afb", "#7c3100", "#00676f", "#82007a", "#2f0c65","#0079c9", "#f3db7f" , "#ecf920", "#004b7c", "#b4faff", "#ffe1ff", "#000000", "#bebebe", "#698b69")

data <- read.csv(file = "data.csv",header = TRUE)
data_si <- read.csv(file = "data_silic_acid.csv",header = TRUE)

# Linear model for Redfield ratio
model <- lm(formula = Phosphate ~  Nitrate, data = data)
model
summary(model)

# Correlation of Redfield ratio to Depth to get correlation coefficients:
data <- data %>% mutate(redfield_ratio = Nitrate/Phosphate)
cor.test(data$redfield_ratio, data$Depth, method = "spearman")

# Sililic acid 
model_si <- lm(formula = Phosphate ~ Nitrate, data = data_si)
model_si
summary(model_si) 

#Correlation of Redfield ratio to silici acid to get correlation coefficients:
data_si <- data_si %>% mutate(redfield_ratio = Nitrate/Phosphate)
cor.test(data_si$redfield_ratio, data_si$Silicic.Acid, method = "spearman")

# Individual spearman correlation between phosphate with silicic acid
cor.test(data_si$Phosphate, data_si$Silicic.Acid, method = "spearman")

# Individual spearman correlation between nitrate with silicic acid
cor.test(data_si$Nitrate, data_si$Silicic.Acid, method = "spearman")

# Plots for Figure2:
# Environmental Feature
myColors <- custom_pallet
names(myColors) <- levels(data$Environmental.Feature)
custom_colors <- scale_colour_manual(name = "ENVO Feature", values = myColors)

redfield_feature <- data %>%
  ggplot(aes(x= Nitrate, Phosphate)) +
  theme_light() +
  xlab("Nitrate Concentration (micromolar)") +
  xlim(NA, 45) +
  ylab("Phosphate Concentration (micromolar)") +
  ylim(NA, 3.5) +
  geom_smooth(se = FALSE, method = "glm", color = 'black') +
  geom_point(size= 1, stroke=1, alpha=3/4, aes(colour = factor(Environmental.Feature)), show.legend = FALSE ) +
  custom_colors 
plot(redfield_feature)

path_string = paste(directory,"redfield_feature.pdf",sep="/")

#save plot
ggsave(filename = path_string, 
       plot =  redfield_feature, width = 11, height = 11, units = "cm")

# Biome
myColors <- custom_pallet
names(myColors) <- levels(data$Biome)
custom_colors <- scale_colour_manual(name = "ENVO Biome", values = myColors)

redfield_biome <- data %>%
  ggplot(aes(x= Nitrate, Phosphate)) +
  theme_light() +
  xlab("Nitrate Concentration (micromolar)") +
  xlim(NA, 45) +
  ylab("Phosphate Concentration (micromolar)") +
  ylim(NA, 3.5) +
  geom_smooth(se = FALSE, method = "lm", color = 'black') +
  geom_point(size= 1, stroke=1, alpha=3/4, aes(colour = factor(Biome)), show.legend = FALSE ) +
  custom_colors
plot(redfield_biome)

path_string = paste(directory,"redfield_biome.pdf",sep="/")

#save plot
ggsave(filename = path_string, 
       plot =  redfield_biome, width = 11, height = 11, units = "cm")

# Project
myColors <- custom_pallet
names(myColors) <- levels(data$Project.Name..)
custom_colors <- scale_colour_manual(name = "Project", values = myColors)

redfield_project <- data %>%
  ggplot(aes(x= Nitrate, Phosphate)) +
  theme_light() +
  xlab("Nitrate Concentration (micromolar)") +
  xlim(NA, 45) +
  ylab("Phosphate Concentration (micromolar)") +
  ylim(NA, 3.5) +
  geom_smooth(se = FALSE, method = "lm", color = 'black') +
  geom_point(size= 1, stroke=1, alpha=3/4, aes(colour = factor(Project.Name..)), show.legend = FALSE ) +
  custom_colors
plot(redfield_project)

path_string = paste(directory,"redfield_project.pdf",sep="/")

#save plot
ggsave(filename = path_string, 
       plot =  redfield_project, width = 11, height = 11, units = "cm")

# Silicic acid 
redfield_si_acid <- data_si %>%
  ggplot(aes(x= Nitrate, Phosphate)) +
  geom_point(size=1, stroke=1, alpha=3/4 ) +
  theme_light() +
  xlab("Nitrate Concentration (micromolar)") +
  xlim(NA, 45) +
  ylab("Phosphate Concentration (micromolar)") +
  ylim(NA, 3.5) +
  geom_smooth(se = FALSE, method = "lm", color = 'black') +
  geom_point(size= 1, stroke=1, alpha=3/4, aes(colour = Silicic.Acid), show.legend = TRUE ) +
  scale_colour_gradient(low = "blue", high="red")
plot(redfield_si_acid)

path_string = paste(directory,"redfield_si_acid.pdf",sep="/")

#save plot
ggsave(filename = path_string, 
       plot =  redfield_si_acid, width = 11, height = 11, units = "cm")

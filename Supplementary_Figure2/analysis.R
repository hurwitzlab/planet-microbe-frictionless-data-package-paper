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

# Filter data
data = filter(data, Project.Name != 'OSD' )
data = filter(data, Environmental.Feature != 'None' )
data = filter(data, Environmental.Feature != 'coastal water body' )
data = filter(data, Environmental.Feature != 'marine layer' )
data = filter(data, Environmental.Feature != 'oceanic zone' )
data = filter(data, Environmental.Feature != 'ocean' )
data = filter(data, Environmental.Feature != 'marine photic zone' )
data = filter(data, Environmental.Feature != 'marine wind mixed layer' )

# Temp figure with colors
myColors <- custom_pallet
names(myColors) <- levels(data$Project.Name)
custom_colors <- scale_colour_manual(name = "Project", values = myColors)

# Bin by latitude
data <- data %>% mutate(Lat_bin = case_when(Latitude >= 66.57  & Latitude <= 90 ~ 'polar',
                                            Latitude >= 23.43  & Latitude <= 66.57 ~ 'temperate',
                                            Latitude >= -23.43  & Latitude <= 23.43 ~ 'tropical',
                                            Latitude >= -66.57  & Latitude <= -23.43 ~ 'temperate',
                                            Latitude >= -90  & Latitude <= -66.57 ~ 'polar'))
# Remove polar samples for figure
data = filter(data, Lat_bin != 'polar' )

acidity <- data %>% 
  ggplot(aes(Environmental.Feature, Acidity)) + 
  theme_linedraw() +
  geom_boxplot() +
  geom_jitter(size=1.2, alpha=3/4,width=0, aes(color=factor(Project.Name)), show.legend = FALSE ) +
  custom_colors +
  facet_grid(. ~ Lat_bin, scales = "free", space = "free")
plot(acidity)

path_string = paste(directory,"acidity.pdf",sep="/")

#save plot
ggsave(filename = path_string, 
       plot =  acidity, width = 29.7, height = 10, units = "cm")


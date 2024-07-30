##########################################
## Predator-prey interaction experiments of A. tonsa and A. pseudogonyaulax
## Here:  GDA Excretion/metabolization rate of A. tonsa after feeding on A. pseudogonyaulax
## Published in: 
## All raw-data available on PANGAEA: 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 05.22
## Alfred-Wegener-Institute Helgoland
##########################################

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(
  PMCMRplus,
  pacman,
  ggplot2,
  growthrates,
  scales,
  agricolae,
  stats,
  psych,
  multcomp,
  rstatix,
  devtools,
  BiocManager,
  qvalue,
  data.table,
  dplyr,
  ggforce,
  patchwork,
  scales,
  gridExtra,
  ggforce
)

# Load data
data = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/R-studio/Extra.txt",
  header = TRUE,
  sep = ""
)

# Subset Qualitative Metabolization of GDA by copepods
data1 <-
  data.frame(time = unlist(data[1:5, 1]), data = unlist(data[1:5, 2]))

# General data transformation
data1 <- data1 %>%
  mutate(
    data = as.numeric(gsub(",", ".", data)),
    time = as.numeric(gsub(",", ".", time)),
    time = as.factor(time)
  )

# convert pg/ul of data1 to pg/copepod with 200 copepods per bottle and an extraction volume of 500 uL methanol
data1$data <- data1$data * (500 / 200)

# introduce grouping variable for plot
data1$group <- c("Initial", rep("Experiment", each = 4))

data1$group <-
  factor(data1$group, levels = c("Initial", "Experiment"))

# ggplot2:
dodge <- position_dodge(width = 0.9)

A <- ggplot(data1, aes(x = time, y = data)) +
  
  geom_col(position = position_dodge2(width = 0.4, preserve = "single")) +
  
  ggtitle(expression(paste(
    "GDA quota of copepods after feeding on ",
    italic("A. pseudogonyaulax"),
    ""
  ))) +
  
  theme_classic() +
  
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.background = element_rect(fill = 'white')
  ) +
  
  scale_y_continuous(expression("Toxin content (pg GDA copepod" ^ -1 * ")")) +
  
  xlab("Time (days)")

# Final <-
#   A + ggforce::facet_row( ~ group, scales = "free", space = "free") + theme(strip.text.x =
#                                                                               element_blank())

# ggsave(
#   "Qual_Meta_2.png",
#   plot = Final,
#   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/AP2",
#   dpi = 300,
#   width = 18,
#   height = 10,
#   units = "cm"
# )

# Garbage collection: call after large objects have been removed
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()

rm(list = ls())

.rs.restartR()

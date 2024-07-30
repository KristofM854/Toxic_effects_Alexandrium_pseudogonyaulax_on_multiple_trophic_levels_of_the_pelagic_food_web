##########################################
## Predator-prey interaction experiments of copepods and A. pseudogonyaulax
## Here: Analysis of copepodamide toxin induction experiment
## Published in: 
## All raw-data available on PANGAEA: 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 03.24
## Alfred-Wegener-Institute Helgoland
##########################################

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Detach all loaded packages
invisible(lapply(
  paste0('package:', names(sessionInfo()$otherPkgs)),
  detach,
  character.only = TRUE,
  unload = TRUE
))

# Load Windows Fonts and define Times font 
pacman::p_load(extrafont, NCmisc)
# 
# packages <-
#   list.functions.in.file("C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP2-Helgoland\\AP2\\R-files\\GDA.R",
#                          alphabetic = TRUE) # set to your filepath
# summary(packages)
loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

# Check for used packages in the file
# list.functions.in.file("C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP2-Helgoland\\AP2\\GDA.R")

# Load data; replace commas; change strain format
data = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/R-studio/copepodamide_toxin_induction.txt",
  header = TRUE,
  sep = ""
)

# introduce grouping variable 
pacman::p_load(tidyr, stringr, dplyr)

# experiment grouping variable - Toxin induction was once tested with an extract 
# from C. finmarchicus from E. Selander and once with a self-made A. tonsa extract
data <- data %>% mutate(Exp = as.factor(c(
  rep("C. finmarchicus", length.out = 52),
  rep("A. tonsa", length.out = nrow(data) - 52)
)))

data <- data %>% 
  mutate(group = as.factor(str_split(Sample, "_", simplify = TRUE)[, 1])) %>%
  filter(group != "GDA")

# aggregate data
data_agg <- data %>%
  group_by(group, Exp) %>%
  dplyr::summarise(GDA_cell_mean = mean(GDA_cell, na.rm = TRUE),
                   GDA_cell_SD = sd(GDA_cell, na.rm = TRUE))

# As point +/- CI for manuscript
pacman::p_load(broom)

# Calculate confidence intervals
confidence_intervals <- data %>%
  group_by(group, Exp) %>%
  do(tidy(t.test(.$GDA_cell, conf.level = 0.95)))

confidence_intervals <- confidence_intervals %>%
  mutate(group = as.character(group),
         group = case_when(
           !group %in% c("Control", "Cop") ~ paste0(group, " nM"),
           TRUE ~ group
         ))

pacman::p_load(forcats, ggplot2)

confidence_intervals$group <- factor(confidence_intervals$group)

confidence_intervals <- confidence_intervals %>%
  mutate(concentration = ifelse(group %in% c("Control", "Cop"), 0, as.numeric(
    gsub("nM", "", as.character(group))
  )))

confidence_intervals$concentration <- as.factor(sprintf("%.1f", confidence_intervals$concentration))

confidence_intervals$Exp <- factor(confidence_intervals$Exp, levels = c("C. finmarchicus", "A. tonsa"))

limits <- aes(ymin = conf.low, ymax = conf.high)

# P1 <- ggplot(confidence_intervals %>% filter(group != "Cop" & Exp != "A. tonsa"), aes(x = concentration, y = estimate)) +
#   geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
#   theme_classic() +
#   scale_x_discrete("Copepodamide<br>concentration (nM)") +
#   scale_y_continuous("Toxin content (pg GDA cell<sup> -1</sup>)", limits = c(7, 16))  +
#   theme(
#     legend.title = element_blank(),
#     strip.background = element_blank(),
#     legend.key.height = unit(1, "cm"),
#     strip.text = element_text(),
#     axis.title.y = element_markdown(),
#     axis.title.x = element_markdown(), 
#     legend.position = "top",
#     panel.background = element_rect(fill = 'white'),
#     axis.text.x = element_markdown(vjust = 0.5, hjust = 0.5, align_heights = T, angle = 45),
#     legend.text = element_text(),
#     plot.caption = element_textbox_simple(),
#     plot.subtitle = element_markdown()
#   ) +
#   ggtitle("")
# 
# P2 <- ggplot(confidence_intervals %>% filter(Exp == "A. tonsa"), aes(x = group, y = estimate)) +
#   geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
#   theme_classic() +
#   scale_x_discrete("Copepodamide<br>concentration (nM)", labels = c("0", "<i>A. tonsa<br></i> extract")) +
#   scale_y_continuous("Toxin content (pg GDA cell<sup> -1</sup>)", limits = c(7, 16))  +
#   theme(
#     legend.title = element_blank(),
#     strip.background = element_blank(),
#     legend.key.height = unit(1, "cm"),
#     strip.text = element_text(),
#     axis.title.y = element_markdown(),
#     axis.title.x = element_markdown(), 
#     legend.position = "top",
#     panel.background = element_rect(fill = 'white'),
#     axis.text.x = element_markdown(vjust = 0.5, hjust = 0.5, align_heights = T, angle = 45),
#     legend.text = element_text(),
#     plot.caption = element_textbox_simple(),
#     plot.subtitle = element_markdown()
#   ) +
#   ggtitle("")
# Filtered data
filtered_data <- confidence_intervals %>% filter(!(group == "Cop" & Exp == "C. finmarchicus")) %>% arrange(Exp)

filtered_data$group <- c("0.1", "0.2", "0.5", "1.0", "2.0", "5.0","0", "0", "crude<br>extract")

filtered_data$group <- factor(filtered_data$group, levels = c("0", "0.1", "0.2", "0.5", "1.0", "2.0", "5.0", "crude<br>extract"))

pacman::p_load(ggplot2, ggthemes, ggtext)

P_combined <- ggplot(filtered_data, aes(x = group, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_classic() +
  scale_x_discrete("Copepodamide<br>concentration (nM)") +
  scale_y_continuous("Toxin content<br>(pg GDA cell<sup> -1</sup>)", limits = c(7, 16))  +
  facet_grid(
    ~ Exp,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(
      Exp = function(x) {
        ifelse(x == "C. finmarchicus", "<b>a)</b>",
               ifelse(x == "A. tonsa", "<b>b)</b>", x))
               }
      )
    ) +
      theme(
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.key.height = unit(1, "cm"),
        strip.text = element_markdown(hjust = 0),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown(), 
        legend.position = "top",
        panel.background = element_rect(fill = 'white'),
        axis.text.x = element_markdown(vjust = 0.5, hjust = 0.5, align_heights = TRUE, angle = 45),
        legend.text = element_text(),
        plot.caption = element_textbox_simple(),
        plot.subtitle = element_markdown()
      ) +
      ggtitle("")
    
# 
# ggsave(
#   "GDA_induction_copepodamides.png",
#   P_combined,
#   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/AP2",
#   dpi = 300,
#   width = 2.5,
#   height = 3,
#   units = "in",
#   scaling = 1
# )

# Statistical analysis
# Perform ANOVA and Tukey's HSD test for each level of Exp
res.aov1 <- aov(GDA_cell ~ group, data = data %>% filter(Exp == "C. finmarchicus" & group != "Cop"))
summary(res.aov1)
tukey1 <- TukeyHSD(res.aov1)

res.aov2 <- aov(GDA_cell ~ group, data = data %>% filter(Exp == "A. tonsa"))
summary(res.aov2)
tukey2 <- as.data.frame(TukeyHSD(res.aov2)$group)

# No significant differences in both experiments!!

# Garbage collection: call after large objects have been removed 
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up 
dev.off()

rm(list = ls())

.rs.restartR()




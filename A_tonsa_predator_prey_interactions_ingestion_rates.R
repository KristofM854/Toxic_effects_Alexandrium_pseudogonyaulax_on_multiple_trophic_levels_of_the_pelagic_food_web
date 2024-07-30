##########################################
## Predator-prey interaction experiments of A. tonsa and A. pseudogonyaulax
## Here: Calculation and analysis of ingestion rates
## Published in: 
## All raw-data available on PANGAEA: 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 05.22
## Alfred-Wegener-Institute Helgoland
##########################################

# INSTALL AND LOAD PACKAGES ################################
# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Load Windows Fonts and define Times font
pacman::p_load(extrafont)

loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

# Detach all loaded packages
invisible(lapply(
  paste0('package:', names(sessionInfo()$otherPkgs)),
  detach,
  character.only = TRUE,
  unload = TRUE
))

##############
# Load data
data = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/R-studio/Ingestion_rate.txt",
  header = TRUE,
  sep = ""
)

# Add grouping variable (12 treatments with 4-8 replicates)
data1 <-
  data.frame(data1 = unlist(data), group = as.factor(rep(c(1:12), each =
                                                           length(data$L2D2))))
# replace commas with points; introduce strain column
pacman::p_load(stringr)

data1$data1 <- as.numeric(gsub(",", ".", data1$data1))
data1$strain <- str_sub(rownames(data1), 1, 4)

# change strain information for bialgal treatment (A. pseudogonyaulax and Rhodomonas salina)
for (i in seq(1, length(data1$group), 1)) {
  if (data1$group[i] == 4 |
      data1$group[i] == 8 |
      data1$group[i] == 12) {
    data1$strain[i] = "L4B1 + R. salina"
  }
}

# convert ingestion rate (cells/hour) to (carbon/hour); raw data in excel-file "POC/PON"; data represents mean of n = 3-5 samples
# Order: L2D2, L4B1, L4B9, L4-B1+Rho
# All N4 (01.06): C [ng] pro Zelle: 4.183 , 4.726, 4.267
# All C4 (27.5 / 02.06): C [ng] pro Zelle: 3.977 , 4.654, 3.956,
# All copepods (30.5): C [ng] pro Zelle: 4.610, 4.654, 3.543
# L4B1_Rho: A.p.: 4.726; Rho: 0.03923
# Include correction factor in one row in data1
data1$poc <-
  rep(
    c(
      4.183,
      4.726,
      4.267,
      4.726,
      3.977,
      4.654,
      3.956,
      4.726,
      4.610,
      4.654,
      3.543,
      4.726
    ),
    each = length(data$L2D2)
  )
data1$data <- data1$data1 * data1$poc

# subgrous of N4/C4/Cop: Order is always L2D2, L4B1, L4B9, L4B1/Rho and N4/C4/Cop
pacman::p_load(tidyr, dplyr, rstatix, forcats)

N4 <-
  data1 %>% subset(group == 1 |
                     group == 2 |
                     group == 3 |
                     group == 4) %>% convert_as_factor(group, strain) %>% mutate(strain = fct_relevel(strain, c("L2D2", "L4B1", "L4B9", "L4B1 + R. salina")))
C4 <-
  data1 %>% subset(group == 5 |
                     group == 6 |
                     group == 7 |
                     group == 8) %>% convert_as_factor(group, strain) %>% mutate(strain = fct_relevel(strain, c("L2D2", "L4B1", "L4B9", "L4B1 + R. salina")))
Cop <-
  data1 %>% subset(group == 9 |
                     group == 10 |
                     group == 11 |
                     group == 12) %>% convert_as_factor(group, strain) %>% mutate(strain = fct_relevel(strain, c("L2D2", "L4B1", "L4B9", "L4B1 + R. salina")))


data_all <- full_join(full_join(N4, C4), Cop)

data_all <- data_all %>% group_by(strain, group) %>% mutate(replicate = rep(1:length(data))) %>% drop_na(data1)

# # Data export for pangaea
# write.table(data_all, "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/PANGAEA/Ingestion_rate.txt", sep = "\t", row.names = FALSE)

# Check for outliers and remove them #####
# outliers_mad: Median absolute deviation --> Only accetable for n>5
# dixon.test: Dixon outlier test. Also applicable for smaller datasets
pacman::p_load(Routliers, outliers)

# Initiate empty vector to store outlier results
outliers <- c()
outliers2 <- c()

for (each_group in unique(N4$group)) {
  MAD1 <-
    outliers_mad( # Detecting univariate outliers using the robust median absolute deviation
      subset(N4, N4$group == each_group)$data,
      b = 1.4826,
      # constant depending on the underlying data distribution --> value equals an assumed normal distribution
      threshold = 3,
      # number of MAD considered as a threshold to consider a value an outlier
      na.rm = TRUE
    )
  
  dix1 <-
    if (dixon.test(subset(N4, N4$group == each_group)$data)$p.val < 0.05) {
      # only extract outliers from dixon.test if p-value < 0.05
      parse_number(dixon.test(subset(N4, N4$group == each_group)$data)$alternative) # extract outliers
    }
  if (!is.null(MAD1)) {
    # bind MAD results
    MAD1_out <- as.vector(MAD1$outliers)
    outliers <- c(outliers, MAD1_out)
  }
  if (!is.null(dix1)) {
    # bind dixon results
    outliers2 <- c(outliers2, dix1)
  }
}

# remove outliers from dataset
N4_2 <- N4[!(N4$data %in% outliers2), ]
N4 <- N4[!(N4$data %in% outliers), ]

# Initiate empty vector to store outlier results
outliers <- c()
outliers2 <- c()

for (each_group in unique(C4$group)) {
  MAD1 <-
    outliers_mad(
      subset(C4, C4$group == each_group)$data,
      b = 1.4826,
      threshold = 3,
      na.rm = TRUE
    )
  dix1 <-
    if (dixon.test(subset(C4, C4$group == each_group)$data)$p.val < 0.05) {
      parse_number(dixon.test(subset(C4, C4$group == each_group)$data)$alternative)
    }
  if (!is.null(MAD1)) {
    MAD1_out <- as.vector(MAD1$outliers)
    outliers <- c(outliers, MAD1_out)
  }
  if (!is.null(dix1)) {
    outliers2 <- c(outliers2, dix1)
  }
}

# remove outliers from dataset
C4_2 <- C4[!(C4$data %in% outliers2), ]
C4 <- C4[!(C4$data %in% outliers), ]

# Initiate empty vector to store outlier results
outliers <- c()
outliers2 <- c()

for (each_group in unique(Cop$group)) {
  MAD1 <-
    outliers_mad(
      subset(Cop, Cop$group == each_group)$data,
      b = 1.4826,
      threshold = 3,
      na.rm = TRUE
    )
  dix1 <-
    if (dixon.test(subset(Cop, Cop$group == each_group)$data)$p.val < 0.05) {
      parse_number(dixon.test(subset(Cop, Cop$group == each_group)$data)$alternative)
    }
  if (!is.null(MAD1)) {
    MAD1_out <- as.vector(MAD1$outliers)
    outliers <- c(outliers, MAD1_out)
  }
  if (!is.null(dix1)) {
    outliers2 <- c(outliers2, dix1)
  }
}

# remove outliers from dataset
Cop2 <- Cop[!(Cop$data %in% outliers2), ]
Cop <- Cop[!(Cop$data %in% outliers), ]

# Check for equal variances by bartlett test or by levenes test
# Ho: All population variances are equal
levene_test(data ~ group, data = data1)
levene_test(sqrt(data) ~ group, data = data1)

# p > 0.05 so the Ho hypothesis gets not rejected and thus the population variances are equal (after square root transformation!)

# Check wether the data is normally distributed by Shapiro Wilk (more robust for small sample sizes)
# Ho: The data is normally distributed
shapiro <- data1 %>%
  group_by(group) %>%
  na.omit() %>%
  dplyr::summarize(
    statistic = shapiro.test(sqrt(data))$statistic,
    p.data = shapiro.test(sqrt(data))$p.data,
    n = length(data),
    var = var(sqrt(data))
  )

# normal two-way ANOVA for each sub
res.aov1 <- aov(sqrt(data) ~ group, data = N4)
res.aov2 <- aov(sqrt(data) ~ group, data = C4)
res.aov3 <- aov(sqrt(data) ~ group, data = Cop)

# Statistically significant differences between the treatments --> PostHoc Test
# PostHoc Tukey test
tukey.test1 <- TukeyHSD(res.aov1)
tukey.test2 <- TukeyHSD(res.aov2)
tukey.test3 <- TukeyHSD(res.aov3)

# ggplot ####
pacman::p_load(ggplot2, ggtext, ggthemes)

# join the data sets of each A. tonsa lifestage to construct a combined ggplot
data_all <- full_join(full_join(N4, C4), Cop)

# Prepare ggplot parameters
dodge <- position_dodge(width = 0.9)

data_all$group2 <-
  factor(c(rep(c("N4", "C4"), each = 32), rep(c("Cop"), each = 31)), levels = c("N4", "C4", "Cop"))

labels <-
  c("strain A", "strain B", "strain C", "strain B + <i> R. salina </i>")

labels2 = c("N<sub>4</sub>-<br>nauplii", "C<sub>4</sub>-<br>copepodites", "adult<br>copepods")

# add vertical lines to separate A. tonsa lifestages
addline_format <- function(x, ...) {
  gsub('\\s', '\n', x)
}

vline_data <-
  data.frame(x = seq(1.5, 2.5, 1), y = rep(250, length(seq(1.5, 2.5, 1))))

# As point +/- CI for manuscript
pacman::p_load(broom)

# Calculate confidence intervals
confidence_intervals <- data_all %>% 
  rename(life_stage = group2) %>%
  group_by(strain, life_stage) %>%
  do(tidy(t.test(.$data, conf.level = 0.95)))

# View the resulting dataframe with confidence intervals
confidence_intervals

limits <- aes(ymin = conf.low, ymax = conf.high)

colors <- c(colorblind_pal()(4))[c(4, 2, 3)]

P_all_manuscript_ingestion <-
  ggplot(confidence_intervals %>% filter(strain != "L4B1 + R. salina"),
         aes(y = estimate, x = life_stage)) +

  geom_pointrange(aes(ymin = conf.low, ymax = conf.high, col = strain),
                  position = position_dodge(width = 0.8, preserve = "total")) +

  scale_x_discrete(NULL, labels = labels2) +

  scale_color_manual(values = colors, labels = labels) +

  scale_y_continuous(expression("Ingestion rate (ng C copepod" ^ -1 * " h" ^ -1 * ")")) +
  
  theme_classic() +

  theme(
    legend.title = element_blank(),
    strip.background = element_blank(),
    legend.key.height = unit(0.5, "cm"),
    strip.text = element_text(),
    axis.title = element_text(),
    legend.position = "bottom",
    panel.background = element_rect(fill = 'white'),
    axis.text.x = element_markdown(angle = 45, vjust = 0.5, hjust = 0.5, align_heights = T),
    legend.text = element_markdown(),
    plot.caption = element_textbox_simple(hjust = 0),
    plot.subtitle = element_markdown()
  ) +

  coord_cartesian(ylim = c(0, max(confidence_intervals$estimate * 1.25))) +

  geom_segment(
    data = vline_data,
    aes(
      x = x,
      xend = x,
      y = 0,
      yend = y
    ),
    linetype = "dotted",
    color = "black"
  ) +
  labs(subtitle = "<b>a)</b>")

# # save locally to combine with figure in another script
saveRDS(P_all_manuscript_ingestion, "P_all_manuscript_ingestion.rds")

# # export word document with package citations used in ingestion rate R-files
# pacman::p_load(grateful)
# cite_packages(out.format = "docx", out.dir = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/", pkgs = "Session", out.file = "ingestion_rate_packages")

# Garbage collection: call after large objects have been removed 
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up 
dev.off()

rm(list = ls())

.rs.restartR()
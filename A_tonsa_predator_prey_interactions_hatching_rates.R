# Kristof Möller 05/2022
# Experimental Design:
# A. tonsa eggs were subjected to the supernatant of A. pseudogonyaulax cultures for 48h; egg hatching ratio was determined by counting hatched nauplii
# AWI Helgoland

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Detach all loaded packages
invisible(lapply(
  paste0('package:', names(sessionInfo()$otherPkgs)),
  detach,
  character.only = TRUE,
  unload = TRUE
))

# Load Windows Fonts and define Times font
pacman::p_load(extrafont)

loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

## Load data, unlist dataframe, change , to . and introduce grouping factor

data = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/R-studio/Hatching_rate.txt",
  header = FALSE,
  sep = ""
)

data1 <-
  data.frame(data = unlist(data[2:7]),
             group = rep(as.factor(data$V1), each = 1, by = length(data)))

data1$data <- as.numeric(gsub(",", ".", data1$data))

# Check for outliers and remove them #####
# outliers_mad: Median absolute deviation --> Only accetable for n>5
# dixon.test: Dixon outlier test. Also applicable for smaller datasets
pacman::p_load(Routliers, outliers)

outliers <- c()
outliers2 <- c()

for (each_treatment in unique(data1$group)) {
  MAD1 <-
    outliers_mad(
      subset(data1, data1$group == each_treatment)$data,
      b = 1.4826,
      # constant depending on the underlying data distribution. Value equals an assumed normal distribution
      threshold = 3,
      # number of MAD considered as a threshold to consider a value an outlier
      na.rm = TRUE
    )
  dix1 <-
    if (dixon.test(subset(data1, data1$group == each_treatment)$data)$p.val < 0.05) {
      # only extract outliers from dixon.test if p-value < 0.05
      parse_number(dixon.test(subset(data1, data1$group == each_treatment)$data)$alternative) # extract outlier
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

# Both outlier tests find no outliers!!

# Check for equal variances by levenes test
# Ho: All population variances are equal
pacman::p_load(rstatix)

levene_test(data ~ group, data = data1)
# p >> 0.05 so the Ho hypothesis gets not rejected and thus the population variances are equal

# Check whether the data is normally distributed by Shapiro Wilk (more robust for small sample sizes)
# Ho: The data is normally distributed
pacman::p_load(tidyr, dplyr)

test <- data1 %>%
  group_by(group) %>%
  dplyr::summarize(statistic = shapiro.test(data)$statistic,
                   p.data = shapiro.test(data)$p.value)

# all normally distributed!

# normal two-way ANOVA
res.aov <-
  anova_test(data = data1, dv = data, between = group) # differences in egg hatching between each treatment (strains and control)

res.aov2 <- aov(data = data1, data ~ group)

# Statistically significant differences between the treatments --> PostHoc Test

# PostHoc Tukey test

tukey.test <-
  as.data.frame(TukeyHSD(aov(data = data1, data ~ group))$group)

# extract the treatments that were compared from the tukey.test
pacman::p_load(stringr)

tukey.test$group1 <- str_sub(rownames(tukey.test), 1, 5)
tukey.test$group2 <- str_sub(rownames(tukey.test), 7, 13)

# change number format and change values from 0.000 to 0.001 for visualization as p < 0.001 in the graph

tukey.test$`p adj` <- sprintf("%.3f", tukey.test$`p adj`)
tukey.test$`p adj` <-
  ifelse(tukey.test$`p adj` < 0.001, 0.001, tukey.test$`p adj`)

# Calculate cohens d effect size
effsize <- cohens_d(data1, data ~ group)
effsize1 <-
  data.frame(data = c(NA, round(effsize$effsize[1:3], digits = 1)),
             group = c(NA, effsize$group2[1:3])) # Include NA for control to match the treatments of the data1 dataframe

# colorblind colors but skip black (first color)
pacman::p_load(ggthemes)

colors <- c(colorblind_pal()(4))[c(1, 4, 2, 3)]

data1$group2 <-
  ifelse(data1$group == "Control", "Control", "A. pseudogonyaulax")

pacman::p_load(ggplot2, ggtext, ggprism)

labels <-
  c("control",
    "strain A",
    "strain B",
    "strain C")

# ggplot paper

# As point +/- CI for manuscript
pacman::p_load(broom)

# Calculate confidence intervals
confidence_intervals <- data1 %>%
  group_by(group) %>%
  do(tidy(t.test(.$data, conf.level = 0.95)))

# View the resulting dataframe with confidence intervals
confidence_intervals
pacman::p_load(glue, ggtext)

F_stat <-  sprintf("%.2f",{res.aov$F[1]})

anova_subtitle <- glue("ANOVA, *F*<sub>{res.aov$DFn[1]},{res.aov$DFd[1]}</sub> = {F_stat}, *p* < 0.001")

hatching_rate <- ggplot(confidence_intervals, aes(y = estimate, x = group)) +
  geom_pointrange(
    aes(col = group, ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = 0.8, preserve = "total"),
    show.legend = F
  ) +
  add_pvalue(
    data = tukey.test[1:3, ],
    label = "< {`p adj`}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 85,
    label.size = 2.75,
    step.increase = 0.125,
    bracket.shorten = 0
  )  +
  scale_color_manual(values = colors, labels = labels) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_markdown(),
    axis.title = element_markdown(),
    axis.text.x = element_markdown(
      color = colors, angle = 45, vjust = 0.5
    ),
    axis.text.y = element_markdown(),
    legend.text = element_markdown(lineheight = 1.5),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    legend.title = element_blank(),
    plot.subtitle = element_markdown(size = 10)
  ) +
  scale_y_continuous("Hatching rate success (%)", limits = c(40, 100)) +
  scale_x_discrete(labels = labels) +
  xlab("") +
  labs(subtitle = anova_subtitle)

ggsave(
  "Hatching_rate_CI.png",
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/AP2", 
  hatching_rate, 
  width = 3.5, height = 4, units = "in", dpi = 300, scaling = 1.35
)

#
# # export word document with package citations used in ingestion rate and GDA R-files
# pacman::p_load(grateful)
# cite_packages(out.format = "docx", out.dir = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/", pkgs = "Session", out.file = "hatching_rate_packages")

# Garbage collection: call after large objects have been removed
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()

rm(list = ls())

.rs.restartR()

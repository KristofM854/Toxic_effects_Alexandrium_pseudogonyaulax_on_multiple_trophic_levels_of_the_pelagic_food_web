##########################################
## Predator-prey interaction experiments of A. tonsa and A. pseudogonyaulax
## Here: Calculation and analysis of goniodomin cell content
## Published in: 
## All raw-data available on PANGAEA: 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 05.22
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
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/R-studio/GDA.txt",
  header = TRUE,
  sep = ""
)

data$data <- as.numeric(gsub(",",".",data$data))
data$Strain_date <- gsub("_","-",data$Strain_date)

## calculate molar ratio of GDA / C
## pg/cell has to be converted to mol/cell first! M(GDA) = 768.941 g/mol; M(C) = 12.011
data$data_norm <-
  ((data$data * 10 ^ -12) / (768.941)) / ((data$C_pg.cell * 10 ^ -12) / (12.011))

# counts <- data %>% group_by(Strain_date) %>% dplyr::summarise(counts = dplyr::n())
# write.table(counts, "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/counts_GDA.txt", sep = "\t", row.names = FALSE)

##subgrous of L2D2;L4B1;L4B9 and convert the idenfication key (Strain_date) to a factor
pacman::p_load(tidyr, dplyr, rstatix)

L2D2 <-
  data %>% subset(
    Strain_date == "L2-D2-C" |
      Strain_date == "L2-D2-N4" |
      Strain_date == "L2-D2-C4" |
      Strain_date == "L2-D2-Cop"
  ) %>% convert_as_factor(Strain_date)

L4B1 <-
  data %>% subset(
    Strain_date == "L4-B1-C" |
      Strain_date == "L4-B1-N4" |
      Strain_date == "L4-B1-C4" |
      Strain_date == "L4-B1-Cop" |
      Strain_date == "L4-B1-RAG"
  ) %>% convert_as_factor(Strain_date)

L4B9 <-
  data %>% subset(
    Strain_date == "L4-B9-C" |
      Strain_date == "L4-B9-N4" |
      Strain_date == "L4-B9-C4" |
      Strain_date == "L4-B9-Cop"
  ) %>% convert_as_factor(Strain_date)

data_all <- full_join(full_join(L2D2, L4B1), L4B9)

data_all <- data_all %>% mutate(strain = stringr::str_sub(data_all$Strain_date, 1, 5),
                                life_stage = stringr::str_sub(data_all$Strain_date, 7, 9))

data_all <- data_all %>% group_by(strain, life_stage) %>% mutate(replicate = rep(1:length(data)), poc = C_pg.cell / 1000) 

## Data export to pangaea
# write.table(data_all, "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/PANGAEA/GDA.txt", sep = "\t", row.names = FALSE)

## Check for potential outliers and remove them #####
pacman::p_load(Routliers, outliers, readr)
outliers <- c()
outliers2 <- c()

for (each_lifestage in unique(L2D2$Strain_date)) {
  dix1 <-
    if (dixon.test(subset(L2D2, L2D2$Strain_date == each_lifestage)$data_norm)$p.val < 0.05) {
      parse_number(dixon.test(
        subset(L2D2, L2D2$Strain_date == each_lifestage)$data_norm
      )$alternative)
    }
  MAD1 <-
    outliers_mad(
      subset(L2D2, L2D2$Strain_date == each_lifestage)$data_norm,
      b = 1.4826,
      threshold = 3,
      na.rm = TRUE
    )
  if (!is.null(MAD1)) {
    MAD1_out <- as.vector(MAD1$outliers)
    outliers <- c(outliers, MAD1_out)
  }
  if (!is.null(dix1)) {
    outliers2 <- c(outliers2, dix1)
  }
}

L2D2 <- L2D2[!(L2D2$data_norm %in% outliers),]

outliers <- c()
outliers2 <- c()

for (each_lifestage in unique(L4B1$Strain_date)) {
  dix1 <-
    if (dixon.test(subset(L4B1, L4B1$Strain_date == each_lifestage)$data_norm)$p.val < 0.05) {
      parse_number(dixon.test(
        subset(L4B1, L4B1$Strain_date == each_lifestage)$data_norm
      )$alternative)
    }
  MAD1 <-
    outliers_mad(
      subset(L4B1, L4B1$Strain_date == each_lifestage)$data_norm,
      b = 1.4826,
      threshold = 3,
      na.rm = TRUE
    )
  if (!is.null(MAD1)) {
    MAD1_out <- as.vector(MAD1$outliers)
    outliers <- c(outliers, MAD1_out)
  }
  if (!is.null(dix1)) {
    outliers2 <- c(outliers2, dix1)
  }
}

L4B1 <- L4B1[!(L4B1$data_norm %in% outliers),]

outliers <- c()
outliers2 <- c()

for (each_lifestage in unique(L4B9$Strain_date)) {
  dix1 <-
    if (dixon.test(subset(L4B9, L4B9$Strain_date == each_lifestage)$data_norm)$p.val < 0.05) {
      parse_number(dixon.test(
        subset(L4B9, L4B9$Strain_date == each_lifestage)$data_norm
      )$alternative)
    }
  MAD1 <-
    outliers_mad(
      subset(L4B9, L4B9$Strain_date == each_lifestage)$data_norm,
      b = 1.4826,
      threshold = 3,
      na.rm = TRUE
    )
  if (!is.null(MAD1)) {
    MAD1_out <- as.vector(MAD1$outliers)
    outliers <- c(outliers, MAD1_out)
  }
  if (!is.null(dix1)) {
    outliers2 <- c(outliers2, dix1)
  }
}

L4B9 <- L4B9[!(L4B9$data_norm %in% outliers),]

# Visually inspect boxplots of each strain #####
pacman::p_load(ggplot2, ggpubr)
# 
# box1 <- ggplot(L2D2, aes(x = Strain_date, y = data_norm)) +
#   geom_boxplot() +
#   ggtitle("L2D2") +
#   theme(
#     plot.title = element_blank(),
#     legend.title = element_blank(),
#     strip.background = element_blank(),
#     axis.text.x = element_text(angle = 90),
#     strip.text = element_text(size = 14),
#     axis.title = element_text(size = 14),
#     panel.background = element_rect(fill = 'white')
#   ) +
#   ylab("GDA:C")
# 
# box2 <- ggplot(L4B1, aes(x = Strain_date, y = data_norm)) +
#   geom_boxplot() +
#   ggtitle("L2D2") +
#   theme(
#     plot.title = element_blank(),
#     legend.title = element_blank(),
#     strip.background = element_blank(),
#     axis.text.x = element_text(angle = 90),
#     strip.text = element_text(size = 14),
#     axis.title = element_text(size = 14),
#     panel.background = element_rect(fill = 'white')
#   ) +
#   ylab("GDA:C")
# 
# box3 <- ggplot(L4B9, aes(x = Strain_date, y = data_norm)) +
#   geom_boxplot() +
#   ggtitle("L2D2") +
#   theme(
#     plot.title = element_blank(),
#     legend.title = element_blank(),
#     strip.background = element_blank(),
#     axis.text.x = element_text(angle = 90),
#     strip.text = element_text(size = 14),
#     axis.title = element_text(size = 14),
#     panel.background = element_rect(fill = 'white')
#   ) +
#   ylab("GDA:C")
# 
# ggarrange(box1, box2, box3, ncol = 1)

####
# Try typical data transformations to reduces skewness of data
# Test for equal variances and normality of each subgroup #####
pacman::p_load(stats)
L2D2 %>% dplyr::summarize(normal_lev = levene_test(data = L2D2, data_norm ~ Strain_date)$p,
                          log_lev = levene_test(data = L2D2, log(data_norm) ~ Strain_date)$p,
                          rezi_lev = levene_test(data = L2D2, 1/data_norm ~ Strain_date)$p)
## --> reziprok transformation leads to equal variances as determined by the levene_test!

L4B1 %>% dplyr::summarize(normal_lev = levene_test(data = L4B1, data_norm ~ Strain_date)$p,
                          log_lev = levene_test(data = L4B1, log(data_norm) ~ Strain_date)$p,
                          rezi_lev = levene_test(data = L4B1, 1/data_norm ~ Strain_date)$p)
## --> log transformation leads to equal variances as determined by the levene_test!

L4B9 %>% dplyr::summarize(normal_lev = levene_test(data = L4B9, data_norm ~ Strain_date)$p,
                          log_lev = levene_test(data = L4B9, log(data_norm) ~ Strain_date)$p,
                          rezi_lev = levene_test(data = L4B9, 1/data_norm ~ Strain_date)$p)
## --> log transformation leads to equal variances as determined by the levene_test!

L2D2 %>%
  group_by(Strain_date) %>%
  dplyr::summarize(
    statistic = shapiro.test(log(data))$statistic,
    p.data = shapiro.test(log(data))$p.value,
    n = length(data),
    var = var(log(data))
  )

L4B1 %>%
  group_by(Strain_date) %>%
  dplyr::summarize(
    statistic = shapiro.test(data)$statistic,
    p.data = shapiro.test(data)$p.value,
    n = length(data),
    var = var(data)
  )

L4B9 %>%
  group_by(Strain_date) %>%
  dplyr::summarize(
    statistic = shapiro.test(data)$statistic,
    p.data = shapiro.test(data)$p.value,
    n = length(data),
    var = var(data)
  )
#######
## anova test of each subgroup, but with transformed data as determined before 

L2D2$data_norm2 <- 1/L2D2$data_norm
L4B1$data_norm2 <- log(L4B1$data_norm)
L4B9$data_norm2 <- log(L4B9$data_norm)

res.aov1 <- anova_test(L2D2, dv = data_norm2, between = Strain_date)
res.aov2 <- anova_test(L4B1 %>% filter(Strain_date != "L4-B1-RAG"), dv = data_norm2, between = Strain_date)
res.aov3 <- anova_test(L4B9, dv = data_norm2, between = Strain_date)

# cohens d effect sizes

L2D2$data_norm2 <- 1 / L2D2$data_norm
L4B1$data_norm2 <- log(L4B1$data_norm)
L4B9$data_norm2 <- log(L4B9$data_norm)

cohens_d(data_norm2 ~ Strain_date, data = L2D2, hedges.correction = T)
cohens_d(data_norm2 ~ Strain_date, data = L4B1, hedges.correction = T)
cohens_d(data_norm2 ~ Strain_date, data = L4B9, hedges.correction = T)

## Tukeys test for each strain, as requirements are fulfilled after data transformation #####
pacman::p_load(stringr)

res.aov1a <- aov(data = L2D2, 1/data_norm ~ Strain_date)
res.aov2a <- aov(data = L4B1 %>% filter(Strain_date != "L4-B1-RAG"), log(data_norm) ~ Strain_date)
res.aov3a <- aov(data = L4B9, log(data_norm) ~ Strain_date)

combinations <- str_split(rownames(as.data.frame(TukeyHSD(res.aov1a)[[1]])), pattern = "-")
combinations2 <- c()
combinations3 <- c()

for (each_combination in seq(1, length(combinations), 1)) {
  combinations2 = rbind(combinations2, combinations[[each_combination]][3])
  combinations3 = rbind(combinations3, combinations[[each_combination]][6])
}

WX_L2D2 <-
  data.frame(
    p = sprintf("%.3f", as.data.frame(TukeyHSD(res.aov1a)[[1]])$'p adj'),
    group1 = combinations2,
    group2 = combinations3
  )
WX_L2D2$p2 <- ifelse(WX_L2D2$p < 0.001, 0.001, WX_L2D2$p)

combinations <-
  str_split(rownames(as.data.frame(TukeyHSD(res.aov2a)[[1]])), pattern = "-")
combinations2 <- c()
combinations3 <- c()

for (each_combination in seq(1, length(combinations), 1)) {
  combinations2 = rbind(combinations2, combinations[[each_combination]][3])
  combinations3 = rbind(combinations3, combinations[[each_combination]][6])
}

WX_L4B1 <-
  data.frame(
    p = sprintf("%.3f", as.data.frame(TukeyHSD(res.aov2a)[[1]])$'p adj'),
    group1 = combinations2,
    group2 = combinations3
  )
WX_L4B1$p2 <- ifelse(WX_L4B1$p < 0.001, 0.001, WX_L4B1$p)

combinations <-
  str_split(rownames(as.data.frame(TukeyHSD(res.aov3a)[[1]])), pattern = "-")
combinations2 <- c()
combinations3 <- c()

for (each_combination in seq(1, length(combinations), 1)) {
  combinations2 = rbind(combinations2, combinations[[each_combination]][3])
  combinations3 = rbind(combinations3, combinations[[each_combination]][6])
}

WX_L4B9 <-
  data.frame(
    p = sprintf("%.3f", as.data.frame(TukeyHSD(res.aov3a)[[1]])$'p adj'),
    group1 = combinations2,
    group2 = combinations3
  )
WX_L4B9$p2 <- ifelse(WX_L4B9$p < 0.001, 0.001, WX_L4B9$p)

L4B1$Strain_date <- gsub("RAG", "Gaze", x = L4B1$Strain_date)

L2D2$strain <- substr(L2D2$Strain_date, 1, 5)
L2D2$life_stage <-
  factor(substr(L2D2$Strain_date, 7, 9), levels = c("C", "N4", "C4", "Cop"))

L4B1$strain <- substr(L4B1$Strain_date, 1, 5)
L4B1$life_stage <-
  factor(substr(L4B1$Strain_date, 7, 10),
         levels = c("C", "N4", "C4", "Cop", "Gaze"))

L4B9$strain <- substr(L4B9$Strain_date, 1, 5)
L4B9$life_stage <-
  factor(substr(L4B9$Strain_date, 7, 9), levels = c("C", "N4", "C4", "Cop"))

scientific_10 <- function(x) {
  parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

labels = c(
  "<i>A. pseudogonyaulax</i><br>control",
  "N<sub>4</sub>-nauplii",
  "C<sub>4</sub>-copepodites",
  "adult copepods",
  "gaze"
)

labels2 = c(
  "<i>A. pseudogonyaulax</i><br>control",
  "N<sub>4</sub>-<br>nauplii",
  "C<sub>4</sub>-<br>copepodites",
  "adult<br>copepods"
)

pacman::p_load(ggplot2, ggthemes, ggtext, ggprism)

dodge <- position_dodge(width = 0.9)

data_all <- full_join(full_join(L2D2, L4B1), L4B9)

vline_data <- data.frame(x = seq(1.5, 4.5, 1), y = rep(2*10^-4, length(seq(1.5, 4.5, 1))))

# colorblind colors but skip black (first color)
colors <- c(colorblind_pal()(4))[c(4,2,3)]

# For paper:
anova_subtitle <- "<b>b)</b> ANOVA, *F*<sub>3,22</sub> = 40-1343, *p* < 0.001"

vline_data <-
  data.frame(x = seq(1.5, 3.5, 1), y = rep(45, length(seq(1.5, 3.5, 1))))

labels = c("no grazer", "N<sub>4</sub>-<br>nauplii", "C<sub>4</sub>-<br>copepodites", "adult<br>copepods")

# As point +/- CI for manuscript
pacman::p_load(broom)

# Calculate confidence intervals
confidence_intervals <- data_all %>%
  group_by(strain, life_stage) %>%
  do(tidy(t.test(.$data, conf.level = 0.95)))

# View the resulting dataframe with confidence intervals
confidence_intervals

limits <- aes(ymin = conf.low, ymax = conf.high)

P_all_manuscript <-
  ggplot(subset(confidence_intervals, confidence_intervals$life_stage != "Gaze"),
         aes(y = estimate, x = life_stage)) +

  geom_pointrange(aes(ymin = conf.low, ymax = conf.high, col = strain),
                  position = position_dodge(width = 0.8, preserve = "total")) +

  add_pvalue(
    data = WX_L4B9[c(3, 1, 2), ],
    label = "p < {p2}",
    tip.length = 0.05,
    fontface = "bold",
    y.position = 45,
    label.size = 2.75,
    step.increase = 0.1,
    bracket.shorten = 0
  ) +

   scale_x_discrete(NULL, labels = labels) +

  scale_color_manual(values = colors, labels = c("strain A", "strain B", "strain C")) +

  scale_y_continuous(expression("Toxin content (pg GDA cell" ^ -1 * ")"), breaks = seq(0, 50, 10)) +

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

  coord_cartesian(ylim = c(0, max(data$data * 1.25))) +

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
  labs(subtitle = anova_subtitle)

# ggsave("GDA_manuscript.png", P_all_manuscript, path="C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/AP2",dpi=300, width=20, height=15, units="cm")

# combine ingestion rates and GDA figure in one for manuscript
pacman::p_load(ggpubr, patchwork)

ingestion_rate <- readRDS("P_all_manuscript_ingestion.rds")

P_combined <- ingestion_rate + P_all_manuscript + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))
ggsave(
  "P_combined2.png",
  P_combined,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/AP2",
  dpi = 300,
  width = 7,
  height = 4,
  units = "in"
)

# export word document with package citations used in ingestion rate and GDA R-files
# pacman::p_load(grateful)
# cite_packages(
#   out.format = "docx",
#   out.dir = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/",
#   pkgs = "Session",
#   out.file = "GDA_packages"
# )

# Garbage collection: call after large objects have been removed 
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up 
dev.off()

rm(list = ls())

.rs.restartR()





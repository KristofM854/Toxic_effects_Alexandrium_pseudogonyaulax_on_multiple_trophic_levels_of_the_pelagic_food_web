##########################################
## Gill cell line assays with cell-free supernatans of A. pseudogonyaulax and purified goniodomins
## Here: filtering of data, Dose-response-curves and plotting of data, EC50 calculation
## Published in: 
## All raw-data available on PANGAEA: 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 05-06.23
## University of Vienna 
##########################################

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Detach all loaded packages
# invisible(lapply(
#   paste0('package:', names(sessionInfo()$otherPkgs)),
#   detach,
#   character.only = TRUE,
#   unload = TRUE
# ))

# Find used packages

#### PACKAGES ####
## Check which packages you actually used for documentation and tidying purposes
pacman::p_load(NCmisc)

# packages <-
#   list.functions.in.file("C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP3\\Ap_all_new.R",
#                          alphabetic = TRUE) # set to your filepath
# summary(packages)

# Load Windows Fonts and define Times font
pacman::p_load(extrafont)

loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

# Load Goniodomins CTB data (metabolic activity) ####
GD_CTB = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/GDs_CTB.txt",
  header = TRUE,
  sep = ""
)

pacman::p_load(rstatix, tidyr, dplyr, stringr)

# General data transformations ####
GD_CTB <-
  GD_CTB %>% convert_as_factor(plate)

GD_CTB <- GD_CTB %>%
  pivot_longer(
    cols = grep("[0-9]", colnames(GD_CTB)),
    names_to = "treat",
    values_to = "data"
  )

GD_CTB <-
  GD_CTB %>%  drop_na(data) %>% 
  mutate(data = as.numeric(sprintf("%.2f", data))) %>% 
  convert_as_factor(treat) %>%
  mutate(treat = str_sub(treat, start = 1, end = -3)) %>%
  arrange(plate, conc_pg_ul) %>%
  mutate(replicate = rep(1:3, length.out = nrow(GD_CTB)))

# Dixon & MAD (median absolute deviation) outlier tests ####
pacman::p_load(Routliers, outliers)

# Aggregate CTB data as the dose response model function requires mean data as input ####
GD_CTB_agg <- GD_CTB %>% 
  aggregate(data ~ conc_pg_ul + treat, FUN = "mean")

GD_CTB_agg$treat <-
  factor(GD_CTB_agg$treat, levels = c("GDA", "GDA_s", "GDA_GDAsa"))

GD_CTB_agg$sd <-
  GD_CTB %>% aggregate(data ~ conc_pg_ul + treat, FUN = sd) %>% 
  pull(data)

# Model dose response curves ####
model_GD <- list()  # Create an empty list to store the models

pacman::p_load(drc, ggplot2)

for (each_treat in GD_CTB_agg$treat) {
  model_GD[[each_treat]] <- drm(
    data = filter(GD_CTB_agg, treat == each_treat),
    data ~ conc_pg_ul,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )
}

list2env(model_GD, envir = globalenv())

# # predictions and confidence intervals.
model.fits_CTB_GD <- expand.grid(conc = seq(1, 500, length = 1000))

# Create an empty list to store model predictions and confidence intervals
model_predictions <- list()

for (each_treat in GD_CTB_agg$treat) {
  model <- drm(
    data = filter(GD_CTB_agg, treat == each_treat),
    data ~ conc_pg_ul,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )
  
  # Create a data frame for predictions and confidence intervals
  model_fits <- expand.grid(conc = seq(1, 500, length = 1000))
  pm <-
    predict(model, newdata = model_fits, interval = "confidence")
  model_fits$p <- pm[, 1]
  model_fits$pmin <- pm[, 2]
  model_fits$pmax <- pm[, 3]
  model_fits$assay <- "CTB"
  model_fits$treat <- as.factor(each_treat)
  
  # Store the model fits in the list
  model_predictions[[each_treat]] <- model_fits
}

list2env(model_predictions, envir = globalenv())

model.fits_CTB_GD_all <- rbind(GDA, GDA_s, GDA_GDAsa)

# extract EC-50 values ####
EC50_GDA <-
  toString(format(model_GD[[1]]$coefficients[[4]], digits = 2))
EC50_GDA_sa <-
  toString(format(model_GD[[3]]$coefficients[[4]], digits = 2))
EC50_GDA_GDA_sa <-
  toString(format(model_GD[[2]]$coefficients[[4]], digits = 2))

# Goniodomins LDH (lytic activity - lactate dehydrogenase assay)
# Load LDH data ####
GD_LDH = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/GDs_LDH.txt",
  header = TRUE,
  sep = ""
)

# General data transformations ####
GD_LDH <-
  GD_LDH %>% convert_as_factor(plate)

GD_LDH <- GD_LDH %>%
  pivot_longer(
    cols = grep("[0-9]", colnames(GD_LDH)),
    names_to = "treat",
    values_to = "data"
  )

GD_LDH <-
  GD_LDH %>%  
  drop_na(data) %>% 
  mutate(data = as.numeric(sprintf("%.2f", data))) %>% 
  convert_as_factor(treat) %>%
  mutate(treat = str_sub(treat, start = 1, end = -3)) %>%
  arrange(plate, conc_pg_ul) %>%
  mutate(replicate = rep(1:3, length.out = nrow(GD_LDH %>% drop_na(data))))

# Aggregate LDH data for dose response curves 
GD_LDH_agg <- GD_LDH %>%
  # filter(!(data %in% outliers_dix)) %>%
  aggregate(data ~ conc_pg_ul + treat, FUN = "mean")

GD_LDH_agg$treat <-
  factor(GD_LDH_agg$treat, levels = c("GDA", "GDA_s", "GDA_GDAsa"))

GD_LDH_agg$sd <-
  GD_LDH %>% 
  aggregate(data ~ conc_pg_ul + treat, FUN = sd) %>% 
  pull(data)

# Model dose response curves ####

model_GD <- list()  # Create an empty list to store the models

for (each_treat in GD_LDH_agg$treat) {
  model_GD[[each_treat]] <- drm(
    data = filter(GD_LDH_agg, treat == each_treat),
    data ~ conc_pg_ul,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )
}

list2env(model_GD, envir = globalenv())

# # predictions and confidence intervals.
model.fits_LDH_GD <- expand.grid(conc = seq(1, 500, length = 1000))

# Create an empty list to store model predictions and confidence intervals
model_predictions <- list()

for (each_treat in GD_LDH_agg$treat) {
  model <- drm(
    data = filter(GD_LDH_agg, treat == each_treat),
    data ~ conc_pg_ul,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )
  
  # Create a data frame for predictions and confidence intervals
  model_fits <- expand.grid(conc = seq(1, 500, length = 1000))
  pm <-
    predict(model, newdata = model_fits, interval = "confidence")
  model_fits$p <- pm[, 1]
  model_fits$pmin <- pm[, 2]
  model_fits$pmax <- pm[, 3]
  model_fits$assay <- "LDH"
  model_fits$treat <- as.factor(each_treat)
  
  # Store the model fits in the list
  model_predictions[[each_treat]] <- model_fits
}

list2env(model_predictions, envir = globalenv())

model.fits_LDH_GD_all <- rbind(GDA, GDA_s, GDA_GDAsa)

GD_all <- rbind(GD_CTB_agg, GD_LDH_agg)

GD_all$assay <-
  as.factor(rep(c("CTB", "LDH"), each = length(GD_all$data) / 2))

for (i in seq(1, length(GD_all$data), 1)) {
  if (GD_all$assay[i] == "LDH") {
    GD_all$data[i] = GD_all$data[i] * 10
  }
}

model.fits_LDH_GD_all <- rbind(GDA, GDA_s, GDA_GDAsa)

# colorblind colors but skip black (first color)
pacman::p_load(ggthemes, ggtext)
colors <- c(colorblind_pal()(4))[c(4, 2, 3)]

labels = c("GDA", "GDA-sa", "GDA + GDA-sa")

P_GD_all <-
  ggplot(data = GD_all, aes(
    x = conc_pg_ul,
    y = data,
    col = treat,
    linetype = assay
  )) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown(),
    legend.text = element_markdown(lineheight = 1.5),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    plot.subtitle = element_markdown()
  ) +
  scale_linetype_manual(
    labels = c("viability", "lytic activity"),
    values = c("solid", "longdash")
  ) +
  scale_color_manual(labels = labels,
                     values = colors) +
  geom_line(data = model.fits_LDH_GD_all, aes(x = conc, y = p * 10), linewidth = 1) +
  geom_line(data = model.fits_CTB_GD_all, aes(x = conc, y = p), linewidth = 1) +
  ggtitle("DRC L2D2 all replicates; filtered data; all plates") +
  ylab("Gill cell viability (%)") +
  xlab("Goniodomin conc. (pg <i>&mu;</i>L<sup>-1</sup>)") +
  scale_x_log10(limits = c(50, 500)) +
  scale_y_continuous(sec.axis = sec_axis(trans =  ~ . / 10 ,
                                         name = expression("Lytic activity (%)"))) +
  guides(linetype = "none") +
  labs(subtitle = "<b>b)</b>")

# for graphical abstract
P_GD_all_ga <-
  ggplot(data = GD_all, aes(
    x = conc_pg_ul,
    y = data,
    col = treat,
    linetype = assay
  )) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown(),
    legend.text = element_markdown(lineheight = 1.5),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    plot.subtitle = element_text(),
    legend.box.margin = margin(-10, 0, 0, 0)
  ) +
  scale_linetype_manual(
    labels = c("viability", "lytic activity"),
    values = c("solid", "longdash")
  ) +
  scale_color_manual(labels = labels,
                     values = colors) +
  geom_line(data = model.fits_LDH_GD_all, aes(x = conc, y = p * 10), linewidth = 1) +
  geom_line(data = model.fits_CTB_GD_all, aes(x = conc, y = p), linewidth = 1) +
  ggtitle("DRC L2D2 all replicates; filtered data; all plates") +
  ylab("Gill cell viability (%)") +
  xlab("Goniodomin conc. (pg <i>&mu;</i>L<sup>-1</sup>)") +
  scale_x_log10(limits = c(50, 500)) +
  scale_y_continuous(sec.axis = sec_axis(trans =  ~ . / 10 ,
                                         name = expression("Lytic activity (%)"))) +
  guides(linetype = "none") +
  labs(subtitle = expression(bold("b" %->% "marginal effect")))

# Graded EC50 values from goniodomins:
EC50_GDA <-
  toString(format(model_GD[[1]]$coefficients[[4]], digits = 2))
EC50_GDA_sa <-
  toString(format(model_GD[[3]]$coefficients[[4]], digits = 2))
EC50_GDA_GDA_sa <-
  toString(format(model_GD[[2]]$coefficients[[4]], digits = 2))

# CTB L2D2
L2D2_CTB = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/L2_D2_raw_CTB.txt",
  header = TRUE,
  sep = ""
)

## Dixon outlier tests
outliers <- c()
outliers_dix <- c()

for (each_treat in unique(L2D2_CTB$treat)) {
  for (each_cell_count in unique(L2D2_CTB$cell_count)) {
    for (each_plate in unique(L2D2_CTB$plate)) {
      outlier_sub <- L2D2_CTB %>% ungroup() %>%
        filter(cell_count == each_cell_count, plate == each_plate, treat == each_treat) %>%
        drop_na(data)
      
      if (nrow(outlier_sub) > 5) {
        MAD1 <- outliers_mad(outlier_sub$data, b = 1.4826, threshold = 3, na.rm = TRUE)
        dix <- dixon.test(outlier_sub$data, opposite = TRUE, two.sided = TRUE)
        
        if (!is.null(MAD1$outliers)) {
          outliers <- c(outliers, MAD1$outliers)
        }
        
        if (dix$p.value < 0.05) {
          dix.out <- as.numeric(sub(".*?([0-9.]+).*", "\\1", dix$alternative))
          outliers_dix <- c(outliers_dix, dix.out)
        }
      }
    }
  }
}

L2D2_CTB <-
  L2D2_CTB %>% 
  filter(!data %in% outliers_dix) %>% 
  convert_as_factor(treat, plate) %>% 
  mutate(strain = "L2-D2") %>%
  filter(
    !(plate == 6 & cell_count == 500) &
      !(plate == 6 & cell_count == 250) &
      !(plate == 8 & cell_count == 250) &
      !(plate == 8 & cell_count == 500)
  ) 

L2D2_CTB_agg <- L2D2_CTB %>%
  aggregate(data ~ cell_count, FUN = "mean") %>% 
  mutate(strain = "L2-D2")

L2D2_CTB_agg$sd <- L2D2_CTB %>%
  aggregate(data ~ cell_count, FUN = sd) %>%
  pull(data)

model_L2D2 <-
  drm(
    data = L2D2_CTB_agg,
    data ~ cell_count,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )

# # predictions and confidence intervals.
model.fits_CTB_2D2 <-
  expand.grid(conc = seq(1, 6000, length = 1000))
# new data with predictions
pm <-
  predict(model_L2D2, newdata = model.fits_CTB_2D2, interval = "confidence")
model.fits_CTB_2D2$p <- pm[, 1]
model.fits_CTB_2D2$pmin <- pm[, 2]
model.fits_CTB_2D2$pmax <- pm[, 3]

# CTB L4-B1
L4B1_CTB = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/L4_B1_raw_CTB.txt",
  header = T,
  sep = ""
)

L4B1_CTB <-
  L4B1_CTB %>% convert_as_factor(plate)

L4B1_CTB <- L4B1_CTB %>%
  pivot_longer(
    cols = matches("^[A-C][0-9]"),
    names_to = "treat",
    values_to = "data"
  )

# Temporarily split dataset as 250-2000 cells/mL are the same for all, but unconcentrated
# supernatants have different cell density equivalents. Makes pivoting easier
L4B1_CTB_1 <- L4B1_CTB %>% filter(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b) %>% 
  dplyr::rename(cell_count = cell_counts_c)

L4B1_CTB_2 <- L4B1_CTB %>% filter(!c(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000))

L4B1_CTB_2 <- L4B1_CTB_2 %>% mutate(cell_count = c(cell_counts_a[1:3], cell_counts_b[1:3], cell_counts_c[1:3])) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b, -cell_counts_c)  
  
L4B1_CTB <- full_join(L4B1_CTB_1, L4B1_CTB_2) %>%
  mutate(treat = str_sub(treat, 1, 1))

## Dixon outlier tests
outliers <- c()
outliers_dix <- c()

for (each_treat in unique(L4B1_CTB$treat)) {
  for (each_cell_count in unique(L4B1_CTB$cell_count)) {
    for (each_plate in unique(L4B1_CTB$plate)) {
      outlier_sub <- L4B1_CTB %>% ungroup() %>%
        filter(cell_count == each_cell_count, plate == each_plate, treat == each_treat) %>%
        drop_na(data)
      
      if (nrow(outlier_sub) > 5) {
        MAD1 <- outliers_mad(outlier_sub$data, b = 1.4826, threshold = 3, na.rm = TRUE)
        dix <- dixon.test(outlier_sub$data, opposite = TRUE, two.sided = TRUE)
        
        if (!is.null(MAD1$outliers)) {
          outliers <- c(outliers, MAD1$outliers)
        }
        
        if (dix$p.value < 0.05) {
          dix.out <- as.numeric(sub(".*?([0-9.]+).*", "\\1", dix$alternative))
          outliers_dix <- c(outliers_dix, dix.out)
        }
      }
    }
  }
}

L4B1_CTB <-
  L4B1_CTB %>% 
  filter(!data %in% outliers_dix) %>% 
  convert_as_factor(treat, plate) %>% 
  mutate(strain = "L4-B1") 

L4B1_CTB_agg <- L4B1_CTB %>%
  aggregate(data ~ cell_count, FUN = "mean") %>% 
  mutate(strain = "L4-B1")

L4B1_CTB_agg$sd <-
  L4B1_CTB %>% 
  aggregate(data ~ cell_count, FUN = sd) %>% 
  pull(data)

model_L4B1 <-
  drm(
    data = L4B1_CTB_agg,
    data ~ cell_count,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )

# # predictions and confidence intervals.
model.fits_CTB_4B1 <-
  expand.grid(conc = seq(1, 6000, length = 1000))
# new data with predictions
pm <-
  predict(model_L4B1, newdata = model.fits_CTB_4B1, interval = "confidence")
model.fits_CTB_4B1$p <- pm[, 1]
model.fits_CTB_4B1$pmin <- pm[, 2]
model.fits_CTB_4B1$pmax <- pm[, 3]

# CTB L4-B9
L4B9_CTB = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/L4_B9_raw_CTB.txt",
  header = T,
  sep = ""
)

L4B9_CTB <- L4B9_CTB %>%
  pivot_longer(
    cols = matches("^[A-C][0-9]"),
    names_to = "treat",
    values_to = "data"
  )

# Temporarily split dataset as 250-2000 cells/mL are the same for all, but unconcentrated
# supernatants have different cell density equivalents. Makes pivoting easier
L4B9_CTB_1 <- L4B9_CTB %>% filter(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b) %>% 
  dplyr::rename(cell_count = cell_counts_c)

L4B9_CTB_2 <- L4B9_CTB %>% filter(!c(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000))

L4B9_CTB_2 <- L4B9_CTB_2 %>% mutate(cell_count = c(cell_counts_a[1:3], cell_counts_b[1:3], cell_counts_c[1:3])) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b, -cell_counts_c)  

L4B9_CTB <- full_join(L4B9_CTB_1, L4B9_CTB_2) %>%
  mutate(treat = str_sub(treat, 1, 1))


## Dixon outlier tests
outliers <- c()
outliers_dix <- c()

for (each_treat in unique(L4B9_CTB$treat)) {
  for (each_cell_count in unique(L4B9_CTB$cell_count)) {
    for (each_plate in unique(L4B9_CTB$plate)) {
      outlier_sub <- L4B9_CTB %>% ungroup() %>%
        filter(cell_count == each_cell_count, plate == each_plate, treat == each_treat) %>%
        drop_na(data)
      
      if (nrow(outlier_sub) > 5) {
        MAD1 <- outliers_mad(outlier_sub$data, b = 1.4826, threshold = 3, na.rm = TRUE)
        dix <- dixon.test(outlier_sub$data, opposite = TRUE, two.sided = TRUE)
        
        if (!is.null(MAD1$outliers)) {
          outliers <- c(outliers, MAD1$outliers)
        }
        
        if (dix$p.value < 0.05) {
          dix.out <- as.numeric(sub(".*?([0-9.]+).*", "\\1", dix$alternative))
          outliers_dix <- c(outliers_dix, dix.out)
        }
      }
    }
  }
}

L4B9_CTB <-
  L4B9_CTB %>% 
  filter(!data %in% outliers_dix) %>% 
  convert_as_factor(treat, plate) %>% 
  mutate(strain = "L4-B9")

L4B9_CTB_agg <- L4B9_CTB %>%
  aggregate(data ~ cell_count, FUN = "mean") %>% mutate(strain = "L4-B9")

L4B9_CTB_agg$sd <-
  L4B9_CTB %>% aggregate(data ~ cell_count, FUN = sd) %>% pull(data)

model_L4B9 <-
  drm(
    data = L4B9_CTB_agg,
    data ~ cell_count,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )

# # predictions and confidence intervals.
model.fits_CTB_4B9 <-
  expand.grid(conc = seq(1, 6000, length = 1000))
# new data with predictions
pm <-
  predict(model_L4B9, newdata = model.fits_CTB_4B9, interval = "confidence")
model.fits_CTB_4B9$p <- pm[, 1]
model.fits_CTB_4B9$pmin <- pm[, 2]
model.fits_CTB_4B9$pmax <- pm[, 3]

# all CTB plots together
CTB_all <- rbind(L2D2_CTB_agg, L4B1_CTB_agg, L4B9_CTB_agg)

model.fits_CTB <-
  rbind(model.fits_CTB_2D2, model.fits_CTB_4B1, model.fits_CTB_4B9)
model.fits_CTB$strain <-
  rep(c("L2-D2", "L4-B1", "L4-B9"), each = 1000)

EC50_L2D2 <-
  toString(format(model_L2D2$coefficients[[4]], digits = 2))
EC50_L4B1 <-
  toString(format(model_L4B1$coefficients[[4]], digits = 2))
EC50_L4B9 <-
  toString(format(model_L4B9$coefficients[[4]], digits = 2))

## LDH all strains
# LDH L2D2
L2D2_LDH = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/L2_D2_raw_LDH.txt",
  header = TRUE,
  sep = ""
)
outliers <- c()
outliers_dix <- c()

for (each_treat in unique(L2D2_LDH$treat)) {
  for (each_cell_count in unique(L2D2_LDH$cell_count)) {
    for (each_plate in unique(L2D2_LDH$plate)) {
      outlier_sub <- L2D2_LDH %>%
        ungroup() %>%
        filter(cell_count == each_cell_count, plate == each_plate, treat == each_treat) %>%
        drop_na(data)
      
      if (nrow(outlier_sub) > 5) {
        MAD1 <- outliers_mad(outlier_sub$data, b = 1.4826, threshold = 3, na.rm = TRUE)
        dix <- dixon.test(outlier_sub$data, opposite = TRUE, two.sided = TRUE)
        
        if (!is.null(MAD1$outliers)) {
          outliers <- c(outliers, MAD1$outliers)
        }
        
        if (dix$p.value < 0.05) {
          dix.out <- as.numeric(sub(".*?([0-9.]+).*", "\\1", dix$alternative))
          outliers_dix <- c(outliers_dix, dix.out)
        }
      }
    }
  }
}

L2D2_LDH <-
  L2D2_LDH %>% 
  filter(!data %in% outliers_dix) %>% 
  convert_as_factor(treat, plate) %>% 
  mutate(strain = "L2-D2") %>% filter(
    !(plate == 6 & cell_count == 500) &
      !(plate == 6 & cell_count == 250) &
      !(plate == 8 & cell_count == 250) &
      !(plate == 8 & cell_count == 500)
  )

L2D2_LDH_agg <- L2D2_LDH %>%
  aggregate(data ~ cell_count, FUN = "mean") %>% 
  mutate(strain = "L2-D2")

L2D2_LDH_agg$sd <- L2D2_LDH %>% 
  aggregate(data ~ cell_count, FUN = sd) %>% 
  pull(data)

model_L2D2 <-
  drm(
    data = L2D2_LDH_agg,
    data ~ cell_count,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )

# # predictions and confidence intervals.
model.fits_LDH_2D2 <-
  expand.grid(conc = seq(1, 6000, length = 1000))
# new data with predictions
pm <-
  predict(model_L2D2, newdata = model.fits_LDH_2D2, interval = "confidence")
model.fits_LDH_2D2$p <- pm[, 1]
model.fits_LDH_2D2$pmin <- pm[, 2]
model.fits_LDH_2D2$pmax <- pm[, 3]

# LDH L4-B1
L4B1_LDH = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/L4_B1_raw_LDH.txt",
  header = T,
  sep = ""
)

L4B1_LDH <-
  L4B1_LDH %>% convert_as_factor(plate)

L4B1_LDH <- L4B1_LDH %>%
  pivot_longer(
    cols = matches("^[A-C][0-9]"),
    names_to = "treat",
    values_to = "data"
  )

# Temporarily split dataset as 250-2000 cells/mL are the same for all, but unconcentrated
# supernatants have different cell density equivalents. Makes pivoting easier
L4B1_LDH_1 <- L4B1_LDH %>% filter(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b) %>% 
  dplyr::rename(cell_count = cell_counts_c)

L4B1_LDH_2 <- L4B1_LDH %>% filter(!c(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000))

L4B1_LDH_2 <- L4B1_LDH_2 %>% mutate(cell_count = c(cell_counts_a[1:3], cell_counts_b[1:3], cell_counts_c[1:3])) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b, -cell_counts_c)  

L4B1_LDH <- full_join(L4B1_LDH_1, L4B1_LDH_2) %>%
  mutate(treat = str_sub(treat, 1, 1))

## Dixon outlier tests
outliers <- c()
outliers_dix <- c()
all_outliers <- data.frame()

for (each_treat in unique(L4B1_LDH$treat)) {
  for (each_cell_count in unique(L4B1_LDH$cell_count)) {
    for (each_plate in unique(L4B1_LDH$plate)) {
      outlier_sub <- L4B1_LDH %>%
        ungroup() %>%
        filter(cell_count == each_cell_count, plate == each_plate, treat == each_treat) %>%
        drop_na(data)
      
      if (nrow(unique(outlier_sub)) > 5) {
        MAD1 <- outliers_mad(outlier_sub$data, b = 1.4826, threshold = 3, na.rm = TRUE)
        dix <- dixon.test(outlier_sub$data, opposite = TRUE, two.sided = TRUE)
        
        if (!is.null(MAD1$outliers)) {
          outliers <- c(outliers, MAD1$outliers)
        }
        
        if (dix$p.value < 0.05) {
          dix.out <- as.numeric(sub(".*?([0-9.]+).*", "\\1", dix$alternative))
          outliers_dix <- c(outliers_dix, dix.out)
          outlier_sub <- outlier_sub %>%
            mutate(dix = ifelse(data == dix.out, dix.out, NA))
          all_outliers <- rbind(all_outliers, outlier_sub)
        }
      }
    }
  }
}


L4B1_LDH <-
  L4B1_LDH %>% 
  filter(!data %in% outliers_dix) %>% 
  convert_as_factor(treat, plate) %>% 
  mutate(strain = "L2-D2")

L4B1_LDH_agg <- L4B1_LDH %>%
  aggregate(data ~ cell_count, FUN = "mean") %>% 
  mutate(strain = "L4-B1")

L4B1_LDH_agg$sd <-
  L4B1_LDH %>% 
  aggregate(data ~ cell_count, FUN = sd) %>% 
  pull(data)

model_L4B1 <-
  drm(
    data = L4B1_LDH_agg,
    data ~ cell_count,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )

# # predictions and confidence intervals.
model.fits_LDH_4B1 <-
  expand.grid(conc = seq(1, 6000, length = 1000))
# new data with predictions
pm <-
  predict(model_L4B1, newdata = model.fits_LDH_4B1, interval = "confidence")
model.fits_LDH_4B1$p <- pm[, 1]
model.fits_LDH_4B1$pmin <- pm[, 2]
model.fits_LDH_4B1$pmax <- pm[, 3]

# LDH L4-B9
L4B9_LDH = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/L4_B9_raw_LDH.txt",
  header = T,
  sep = ""
)

L4B9_LDH <-
  L4B9_LDH %>% convert_as_factor(plate)

L4B9_LDH <- L4B9_LDH %>%
  pivot_longer(
    cols = matches("^[A-C][0-9]"),
    names_to = "treat",
    values_to = "data"
  )

# Temporarily split dataset as 250-2000 cells/mL are the same for all, but unconcentrated
# supernatants have different cell density equivalents. Makes pivoting easier
L4B9_LDH_1 <- L4B9_LDH %>% filter(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b) %>% 
  dplyr::rename(cell_count = cell_counts_c)

L4B9_LDH_2 <- L4B9_LDH %>% filter(!c(cell_counts_a < 4000 & cell_counts_b < 4000 & cell_counts_c < 4000))

L4B9_LDH_2 <- L4B9_LDH_2 %>% mutate(cell_count = c(cell_counts_a[1:3], cell_counts_b[1:3], cell_counts_c[1:3])) %>%
  dplyr::select(-cell_counts_a, -cell_counts_b, -cell_counts_c)  

L4B9_LDH <- full_join(L4B9_LDH_1, L4B9_LDH_2) %>%
  mutate(treat = str_sub(treat, 1, 1))

## Dixon outlier tests
outliers <- c()
outliers_dix <- c()
all_outliers <- data.frame()

for (each_treat in unique(L4B9_LDH$treat)) {
  for (each_cell_count in unique(L4B9_LDH$cell_count)) {
    for (each_plate in unique(L4B9_LDH$plate)) {
      outlier_sub <- L4B9_LDH %>%
        ungroup() %>%
        filter(cell_count == each_cell_count, plate == each_plate, treat == each_treat) %>%
        drop_na(data)
      
      if (nrow(unique(outlier_sub)) > 5) {
        MAD1 <- outliers_mad(outlier_sub$data, b = 1.4826, threshold = 3, na.rm = TRUE)
        dix <- dixon.test(outlier_sub$data, opposite = TRUE, two.sided = TRUE)
        
        if (!is.null(MAD1$outliers)) {
          outliers <- c(outliers, MAD1$outliers)
        }
        
        if (dix$p.value < 0.05) {
          dix.out <- as.numeric(sub(".*?([0-9.]+).*", "\\1", dix$alternative))
          outliers_dix <- c(outliers_dix, dix.out)
          outlier_sub <- outlier_sub %>%
            mutate(dix = ifelse(data == dix.out, dix.out, NA))
          all_outliers <- rbind(all_outliers, outlier_sub)
        }
      }
    }
  }
}

L4B9_LDH <-
  L4B9_LDH %>% 
  filter(!data %in% outliers_dix) %>% 
  convert_as_factor(treat, plate) %>% 
  mutate(strain = "L2-D2")

L4B9_LDH_agg <- L4B9_LDH %>%
  aggregate(data ~ cell_count, FUN = "mean") %>% mutate(strain = "L4-B9")

L4B9_LDH_agg$sd <-
  L4B9_LDH %>% 
  aggregate(data ~ cell_count, FUN = sd) %>% 
  pull(data)

model_L4B9 <-
  drm(
    data = L4B9_LDH_agg,
    data ~ cell_count,
    fct = LL.4(),
    na.action = na.omit,
    lowerl = c(-Inf, 0, -Inf, -Inf),
    upperl = c(Inf, Inf, 100, Inf)
  )

# # predictions and confidence intervals.
model.fits_LDH_4B9 <-
  expand.grid(conc = seq(1, 6000, length = 1000))
# new data with predictions
pm <-
  predict(model_L4B9, newdata = model.fits_LDH_4B9, interval = "confidence")
model.fits_LDH_4B9$p <- pm[, 1]
model.fits_LDH_4B9$pmin <- pm[, 2]
model.fits_LDH_4B9$pmax <- pm[, 3]

# all LDH plots together
LDH_all <- rbind(L2D2_LDH_agg, L4B1_LDH_agg, L4B9_LDH_agg)

model.fits_LDH <-
  rbind(model.fits_LDH_2D2, model.fits_LDH_4B1, model.fits_LDH_4B9)
model.fits_LDH$strain <-
  rep(c("L2-D2", "L4-B1", "L4-B9"), each = 1000)
model.fits_LDH[, c(2:4)] <- model.fits_LDH[, c(2:4)] * 3

EC50_L2D2 <-
  toString(format(model_L2D2$coefficients[[4]], digits = 2))
EC50_L4B1 <-
  toString(format(model_L4B1$coefficients[[4]], digits = 2))
EC50_L4B9 <-
  toString(format(model_L4B9$coefficients[[4]], digits = 2))

CTB_LDH_all <-
  rbind(CTB_all, LDH_all) %>% mutate(assay = c(
    rep("CTB", length.out = length(CTB_all$data)),
    rep("LDH", length.out = length(LDH_all$data))
  ))

CTB_LDH_all2 <-
  aggregate(CTB_LDH_all, data ~ strain + cell_count + assay, FUN = "mean")

for (i in seq(1, length(CTB_LDH_all$data), 1)) {
  if (CTB_LDH_all$assay[i] == "LDH") {
    CTB_LDH_all$data[i] = CTB_LDH_all$data[i] * 3
  }
}

# colorblind colors but skip black (first color)
colors <- c(colorblind_pal()(4))[c(4, 2, 3)]

labels2 <-
  c(
    "strain A",
    "strain B",
    "strain C"
  )

model.fits_LDH$assay <- "LDH"
model.fits_CTB$assay <- "CTB"

P3a <-
  ggplot(data = CTB_LDH_all, aes(
    x = cell_count,
    y = data,
    col = strain,
    linetype = assay
  )) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(),
    axis.title = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_markdown(),
    axis.text.y = element_text(),
    legend.text = element_markdown(lineheight = 1.5),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    plot.subtitle = element_markdown()
  ) +
  scale_color_manual(values = colors, labels = labels2) +
  geom_line(data = model.fits_LDH, aes(x = conc, y = p), linewidth = 1) +
  geom_line(data = model.fits_CTB, aes(x = conc, y = p), linewidth = 1) +
  ylab("Gill cell viability (%)") +
  xlab("*A. pseudogonyaulax*<br>cell density (cells mL<sup>-1</sup>)") +
  scale_x_log10(limits = c(100, 10000)) +
  scale_y_continuous(sec.axis = sec_axis(trans =  ~ . / 3 ,
                                         name = expression("Lytic activity (%)"))) +
  guides(linetype = "none") +
  labs(subtitle = "<b>a)</b>")

# for graphical abstract
P3a_ga <-
  ggplot(data = CTB_LDH_all, aes(
    x = cell_count,
    y = data,
    col = strain,
    linetype = assay
  )) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(),
    axis.title = element_markdown(),
    axis.title.x = element_markdown(),
    axis.text.x = element_markdown(),
    axis.text.y = element_text(),
    legend.text = element_markdown(lineheight = 1.5),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    plot.subtitle = element_text(),
    legend.box.margin = margin(-20, 0, 0, 0)
  ) +
  scale_color_manual(values = colors, labels = labels2) +
  geom_line(data = model.fits_LDH, aes(x = conc, y = p), linewidth = 1) +
  geom_line(data = model.fits_CTB, aes(x = conc, y = p), linewidth = 1) +
  ylab("Gill cell viability (%)") +
  xlab("*A. pseudogonyaulax*<br>cell density (cells mL<sup>-1</sup>)") +
  scale_x_log10(limits = c(100, 10000)) +
  scale_y_continuous(sec.axis = sec_axis(trans =  ~ . / 3 ,
                                         name = expression("Lytic activity (%)"))) +
  guides(linetype = "none") +
  labs(subtitle = expression(bold("a+b" %->% "lysis / death")))

# ggsave(
#   "DRC_LDH_CTB.png",
#   P3a,
#   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
#   dpi = 300,
#   width = 20,
#   height = 15,
#   units = "cm"
# )

CTB_all_raw <- full_join(L2D2_CTB, full_join(L4B1_CTB, L4B9_CTB)) %>% arrange(plate, cell_count, strain)

# introduce replicate
CTB_all_raw <- CTB_all_raw %>%
  group_by(plate, cell_count, treat, strain) %>%
  drop_na(data) %>%
  dplyr::summarise(replicate = 1:n(), data = data)


# # Data export for pangaea (CTB)
write.table(
  CTB_all_raw,
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/PANGAEA/CTB_all_raw.txt",
  sep = "\t",
  row.names = FALSE
)
write.table(
  GD_CTB,
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/PANGAEA/CTB_all_raw_GD.txt",
  sep = "\t",
  row.names = FALSE
)

LDH_all_raw <- full_join(L2D2_LDH, full_join(L4B1_LDH, L4B9_LDH))

LDH_all_raw <- LDH_all_raw %>%
  group_by(plate, cell_count, treat, strain) %>% 
  drop_na(data) %>%
  dplyr::summarise(replicate = 1:n(), data = data)

# # Data export for pangaea (LDH)
write.table(
  LDH_all_raw ,
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/PANGAEA/LDH_all_raw.txt",
  sep = "\t",
  row.names = FALSE
)
write.table(
  GD_LDH ,
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/PANGAEA/LDH_all_raw_GD.txt",
  sep = "\t",
  row.names = FALSE
)

# combine overall figure of goniodomins and supernatants in one figure
pacman::p_load(ggpubr, patchwork)

# P_combined <- P3a + P_GD_all + plot_layout(ncol = 2) &
#   theme(legend.position = "bottom",
#         legend.box.margin = margin(-10, -10, -10, -10))

P_combined <- P3a + P_GD_all + plot_layout(ncol = 1) &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))

ggsave(
  "DRC_all.png",
  P_combined,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
  dpi = 300,
  width = 3.5,
  height = 5,
  units = "in"
)

# # Rhodomonas assay A. monilatum and GDs
# Rho_assay = read.csv(
#   "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/Rho_assay.txt",
#   header = T,
#   sep = ""
# )
# 
# Rho_assay <-
#   Rho_assay %>% convert_as_factor(treat, time)
# 
# Rho_assay <-
#   Rho_assay %>% drop_na(data)
# 
# Rho_assay_agg <-
#   aggregate(data = Rho_assay, data ~ conc_pg_ul + treat + time, FUN = "mean")
# Rho_assay_agg$sd <-
#   aggregate(data = Rho_assay, data ~ conc_pg_ul + treat + time, FUN =SD)$data
# 
# # predictions and confidence intervals.
# model.fits_Rho <-
#   expand.grid(conc = c(seq(0.001, 1, length = 1000), seq(1, 250, length = 1000)))
# 
# # Create an empty list to store model predictions and models
# model_predictions <- list()
# model_Rho <- list()
# EC50_Rho <- data.frame()
# 
# pacman::p_load(drc)
# 
# for (each_treat in Rho_assay_agg$treat) {
#   for (each_time in unique(Rho_assay_agg$time)) {
#     if (filter(Rho_assay_agg, treat == each_treat &
#                time == each_time) %>% pull(data) %>% length() > 0)
#       model <- drm(
#         data = filter(Rho_assay_agg, treat == each_treat &
#                         time == each_time),
#         data ~ conc_pg_ul,
#         fct = LL.4(),
#         na.action = na.omit,
#         lowerl = c(-Inf, 0, -Inf, -Inf),
#         upperl = c(Inf, Inf, 100, Inf)
#       )
#     
#     # Create a data frame for predictions and confidence intervals
#     model_fits <-
#       expand.grid(conc = c(seq(0.001, 1, length = 1000), seq(1, 250, length = 1000)))
#     pm <-
#       predict(model, newdata = model_fits, interval = "confidence")
#     model_fits$p <- pm[, 1]
#     model_fits$pmin <- pm[, 2]
#     model_fits$pmax <- pm[, 3]
#     model_fits$treat <- as.factor(each_treat)
#     model_fits$time <- as.factor(each_time)
#     # Store the model fits in the list
#     model_predictions[[paste0(each_treat, "_", each_time)]] <-
#       model_fits
#     # extract EC50 values and store in dataframe
#     model_coef <- data.frame(treat = each_treat, time = each_time, EC50 = format(model$coefficients[[4]], digits = 2))
#     EC50_Rho <- rbind(EC50_Rho, model_coef)
#   }
# }
# 
# EC50_Rho <- unique(EC50_Rho)
# 
# list2env(model_predictions, envir = globalenv())
# 
# Rho_3 <- rbind(GDA_3, A_monilatum_3)
# Rho_24 <- rbind(GDA_24, A_monilatum_24)
# 
# labels = c(expression(italic("A. monilatum")), "GDA")
# colors <- c(colorblind_pal()(4))[c(2, 4)]
# 
# all_Rho <- rbind(Rho_3, Rho_24, GDA_72)
# 
# Rho_all <-
#   ggplot(data = Rho_assay_agg, aes(x = conc_pg_ul, y = data, col = treat)) +
#   geom_pointrange(aes(ymin = data - sd, ymax = data + sd), show.legend = T) +
#   theme_classic() +
#   facet_wrap( ~ time, scales = "free_y") +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank(),
#     axis.title = element_text(size = 12),
#     axis.title.y = element_markdown(size = 12, lineheight = 1.5),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_markdown(size = 12),
#     legend.text = element_text(size = 12),
#     panel.background = element_rect(fill = 'white'),
#     legend.position = "top",
#     legend.key.size = unit(1, "cm"),
#     legend.key.height = unit(1, "cm"),
#     legend.key.width = unit(2, "cm"),
#     legend.text.align = 0.5,
#     legend.title = element_blank()
#   ) +
#   geom_line(data = all_Rho, aes(x = conc, y = p), show.legend = F) +
#   scale_y_continuous(
#     limits = c(0, 125),
#     labels = seq(0, 125, by = 25),
#     breaks = seq(0, 125, by = 25)
#   ) +
#   ylab("Intact *R. salina*<br>(%) of control") +
#   scale_color_manual(labels = labels, values = colors) +
#   scale_x_continuous(
#     trans = 'log10',
#     limits = c(0.001, 1000),
#     labels = scales::number_format()
#   ) +
#   xlab(expression(paste(
#     "Goniodomin conc.", " (pg ", mu, "L" ^ -1 * ")"
#   ))) 
# 
# # ggsave(
# #   "Rho_all.png",
# #   Rho_all,
# #   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
# #   dpi = 300,
# #   width = 16,
# #   height = 10,
# #   units = "cm"
# # )
# 
# # colorblind colors but skip black (first color)
# colors <- c(colorblind_pal()(4))[c(4, 2, 3)]
# labels = c("<i>A. monilatum</i>", "GDA")
# labels2 = c("3 h", "24 h", "72 h")
# # Rho_assay_agg$combined_legend <- interaction(Rho_assay_agg$time, Rho_assay_agg$treat)
# # all_Rho$combined_legend <- interaction(all_Rho$time, all_Rho$treat)
# 
# pacman::p_load("glue")
# 
# Rho_all2 <- ggplot(data = Rho_assay_agg, aes(x = conc_pg_ul, y = data, shape = time, col = treat)) +
#   geom_pointrange(aes(ymin = data - sd, ymax = data + sd), show.legend = F) +
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank(),
#     axis.title = element_markdown(size = 30),
#     axis.title.y = element_markdown(size = 30, lineheight = 1.5, margin = margin(r = 10)),
#     axis.title.x = element_markdown(size = 30),
#     axis.text.x = element_markdown(size = 24),
#     axis.text.y = element_markdown(size = 24, margin = margin(r = 10)),
#     legend.text = element_markdown(size = 30),
#     panel.background = element_rect(fill = 'white'),
#     legend.position = "top",
#     # legend.box = "vertical",
#     legend.key.size = unit(3, "cm"),
#     legend.key.height = unit(1, "cm"),
#     legend.key.width = unit(1, "cm"),
#     legend.text.align = 0.5,
#     legend.title = element_blank(),
#     axis.ticks.length = unit(0.25, "cm"),
#     legend.spacing.x = unit(1, "cm")
#   ) +
#   geom_line(data = all_Rho, aes(x = conc, y = p), show.legend = T, linewidth = 2) +
#   scale_y_continuous(
#     limits = c(0, 125),
#     labels = seq(0, 125, by = 25),
#     breaks = seq(0, 125, by = 25)
#   ) +
#   ylab("intact *R. salina*<br>(%) of control") +
#   scale_color_manual(labels = labels, values = colors[1:2]) +
#   scale_shape_manual(name = "Time",
#                      labels = NULL,
#                      breaks = NULL,
#                      values = c(1, 2, 3)) +
#   scale_x_continuous(
#     trans = 'log10',
#     limits = c(0.001, 1000),
#     labels = scales::number_format()
#   ) +  
#   xlab("Goniodomin conc. (pg &mu;L<sup>-1</sup>)") 
# 
# # +
# #   guides(
# #     col = guide_legend(title = "Treat", direction = "horizontal", ncol = 2, order = 1, byrow = T)
# #     # ,
# #     # shape = guide_legend(title = "Time", direction = "horizontal", ncol = 3, byrow = T, order = 2)
# #   )
# 
# 
# # ggsave(
# #   "Rho_all2.png",
# #   Rho_all2,
# #   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
# #   dpi = 300,
# #   width = 22.5,
# #   height = 15,
# #   units = "cm"
# # )
# 
# Rho_assay <- Rho_assay %>% group_by(treat, time, conc_pg_ul) %>% mutate(replicate = rep(1:length(data)))
# 
# # # Data export for pangaea (Rho assays)
# # write.table(Rho_assay %>% drop_na(data), "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/PANGAEA/Rho_assay.txt", sep = "\t", row.names = FALSE)
# 
# # Rhodomonas assay Alexandrium pseudogonyaulax - done by Francesco (Master student)
# Rho_assay2 = read.csv(
#   "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/Rho_assay_Ap.txt",
#   header = T,
#   sep = ""
# )
# 
# pacman::p_load(tidyverse, dplyr, rstatix)
# 
# Rho_assay2 <-
#   Rho_assay2 %>% convert_as_factor(treat, replicate, strain) %>%
#   mutate(log_conc_Ap = log10(conc_Ap))
# 
# # Convert A. pseudogonyaulax cell counts to GDA equivalents by the concentration found in the supernatant
# # to plot GDA and A. monilatum / A. pseudogonyaulax bioassay results in one figure 
# # divided by 1000 to get ng/mL from pg/mL, which matches pg/uL from the GDs 
# 
# Rho_assay2 <- Rho_assay2 %>% mutate(conc_pg_ul = conc_Ap * 0.74 / 1000)
# 
# # Rhodomonas mean cell counts for normalization corresponding to bioassays with L2-D2 / L4-B1 / L4-B9 see corresponding Excel.´
# Rho_control <- c(259.7777778, 199.3333333, 201.4444444)
# 
# Rho_assay2 <-
#   Rho_assay2 %>% mutate(Rho_control = ifelse(strain == "L2-D2", Rho_control[1],
#                                             ifelse(strain == "L4-B9", Rho_control[2],
#                                                    ifelse(strain == "L4-B1", Rho_control[3], NA))))
# 
# Rho_assay2 <- Rho_assay2 %>% mutate(data = counts/Rho_control*100)
# 
# # # Data export for pangaea (Rho assays)
# # write.table(Rho_assay2, "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/PANGAEA/Rho_assay_Ap.txt", sep = "\t", row.names = FALSE)
# 
# Rho_assay2_agg <-
#   aggregate(data = Rho_assay2, data ~ conc_pg_ul + treat + strain, FUN = "mean")
# Rho_assay2_agg$sd <-
#   aggregate(data = Rho_assay2, data ~ conc_pg_ul + treat + strain, FUN =SD)$data
# 
# # predictions and confidence intervals.
# model.fits_Rho <-
#   expand.grid(conc = c(seq(0, 3, length = 1000)))
# 
# # Create an empty list to store model predictions and models
# model_predictions <- list()
# model_Rho <- list()
# EC50_Rho <- data.frame()
# 
# pacman::p_load(drc)
# 
# for (each_strain in Rho_assay2_agg$strain) {
#       model <- drm(
#         data = filter(Rho_assay2_agg, strain == each_strain),
#         data ~ conc_pg_ul,
#         fct = LL.4(),
#         na.action = na.omit,
#         lowerl = c(-Inf, 0, -Inf, -Inf),
#         upperl = c(Inf, Inf, 100, Inf)
#       )
#     
#     # Create a data frame for predictions and confidence intervals
#     model_fits <-
#       expand.grid(conc = seq(0, 3, length = 1000))
#     pm <-
#       predict(model, newdata = model_fits, interval = "confidence")
#     model_fits$p <- pm[, 1]
#     model_fits$pmin <- pm[, 2]
#     model_fits$pmax <- pm[, 3]
#     model_fits$strain <- as.factor(each_strain)
#     # Store the model fits in the list
#     model_predictions[[paste0(each_strain)]] <-
#       model_fits
#     # extract EC50 values and store in dataframe
#     model_coef <- data.frame(strain = each_strain, EC50 = format(model$coefficients[[4]], digits = 2))
#     EC50_Rho <- rbind(EC50_Rho, model_coef)
#   }
# 
# EC50_Rho <- unique(EC50_Rho)
# 
# list2env(model_predictions, envir = globalenv())
# 
# Rho_fit_all_Ap <- rbind(`L2-D2`, `L4-B1`, `L4-B9`)
# 
# colors <- c(colorblind_pal()(4))[c(1, 2, 4)]
# 
# labels <- unique(Rho_assay2_agg$strain)
# 
# pacman::p_load(ggthemes, ggtext, ggplot2)
# 
# Rho_all_Ap <-
#   ggplot(data = Rho_assay2_agg, aes(x = conc_pg_ul, y = data, col = strain)) +
#   geom_pointrange(aes(ymin = data - sd, ymax = data + sd), show.legend = T) +
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank(),
#     axis.title.x = element_markdown(size = 12),
#     axis.title.y = element_markdown(size = 12, lineheight = 1.5),
#     axis.text.x = element_markdown(size = 12),
#     axis.text.y = element_markdown(size = 12),
#     legend.text = element_text(size = 12),
#     panel.background = element_rect(fill = 'white'),
#     legend.position = "top",
#     legend.key.size = unit(1, "cm"),
#     legend.key.height = unit(1, "cm"),
#     legend.key.width = unit(2, "cm"),
#     legend.text.align = 0.5,
#     legend.title = element_blank()
#   ) +
#   geom_line(data = Rho_fit_all_Ap, aes(x = conc, y = p, col = strain), show.legend = F) +
#   scale_y_continuous(
#     limits = c(0, 125),
#     labels = seq(0, 125, by = 25),
#     breaks = seq(0, 125, by = 25)
#   ) +
#   ylab("intact *R. salina*<br>(%) of control") +
#   scale_color_manual(labels = labels, values = colors) +
#   scale_x_continuous(
#     limits = c(0, 3),
#     labels = scales::number_format()
#   ) +  
#   xlab("Goniodomin conc. (pg &mu;L<sup>-1</sup>)") 
# 
# Rho_all_Ap
# 
# # ggsave(
# #   "Rho_all_Ap.png",
# #   Rho_all_Ap,
# #   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
# #   dpi = 300,
# #   width = 16,
# #   height = 10,
# #   units = "cm"
# # )
# 
# 
# Rho_assay2 <- Rho_assay2 %>% group_by(treat, conc_pg_ul) %>% mutate(replicate = rep(1:length(data)))
# 
# # # Data export for pangaea (Rho assays)
# # write.table(Rho_assay2, "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/PANGAEA/Rho_assay_Ap.txt", sep = "\t", row.names = FALSE)
# 
# # Combine all R. salina bioassays in one figure 
# 
# # Match data and colnames of both aggregated dataframes first
# Rho_assay2_agg <- Rho_assay2_agg %>% mutate(treat = strain, time = 24) %>% dplyr::select(-strain) 
# 
# # Combine bioassay results
# Rho_assay_agg_all <- rbind(Rho_assay_agg, Rho_assay2_agg)
# 
# # Combine modelled dose response curves (DRC)
# 
# # Match data and colnames of both aggregated dataframes first
# Rho_fit_all_Ap <- Rho_fit_all_Ap %>% mutate(treat = strain, time = 24) %>% dplyr::select(-strain) 
# 
# # Combine DRC results
# Rho_fit_all <- rbind(all_Rho, Rho_fit_all_Ap)
# 
# colors <- c(colorblind_pal()(6))
# 
# labels <- unique(Rho_assay_agg_all$treat)
# 
# pacman::p_load(ggthemes, ggtext, ggplot2)
# 
# Rho_all <-
#   ggplot(data = Rho_assay_agg_all, aes(x = conc_pg_ul, y = data, col = treat)) +
#   geom_pointrange(aes(ymin = data - sd, ymax = data + sd), show.legend = T) +
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank(),
#     axis.title.x = element_markdown(size = 12),
#     axis.title.y = element_markdown(size = 12, lineheight = 1.5),
#     axis.text.x = element_markdown(size = 12),
#     axis.text.y = element_markdown(size = 12),
#     legend.text = element_text(size = 12),
#     panel.background = element_rect(fill = 'white'),
#     legend.position = "top",
#     legend.key.size = unit(1, "cm"),
#     legend.key.height = unit(1, "cm"),
#     legend.key.width = unit(2, "cm"),
#     legend.text.align = 0.5,
#     legend.title = element_blank()
#   ) +
#   geom_line(data = Rho_fit_all, aes(x = conc, y = p, col = treat), show.legend = F) +
#   scale_y_continuous(
#     limits = c(0, 125),
#     labels = seq(0, 125, by = 25),
#     breaks = seq(0, 125, by = 25)
#   ) +
#   facet_wrap(~time) +
#   ylab("intact *R. salina*<br>(%) of control") +
#   scale_color_manual(labels = labels, values = colors, guide = guide_legend(nrow = 2)) +
#   scale_x_continuous(
#     trans = 'log10',
#     limits = c(0.001, 1000),
#     labels = scales::number_format()
#   ) +  
#   xlab("Goniodomin conc. (pg &mu;L<sup>-1</sup>)") 
# 
# Rho_all
# 
# ggsave(
#   "Rho_all_new.png",
#   Rho_all,
#   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
#   dpi = 300,
#   width = 16,
#   height = 10,
#   units = "cm"
# )

# Rhodomonas salina bioassays of Alexandrium monilatum, Alexandrium pseudogonyaulax and purified goniodomins (GDA, GDA-sa, GDB, GDA+GDA-sa) combined
# Rhodomonas assay Alexandrium pseudogonyaulax - done by Francesco (Master student)
Rho_assay_all = read.csv(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/raw/Rho_assay_all.txt",
  header = T,
  sep = ""
)

pacman::p_load(tidyverse, dplyr, rstatix)

# Combine L4B9 and A. monilatum datasets
Rho_assay_all <- Rho_assay_all %>%
  mutate(treat = case_when(
    treat == "A_monilatum_old" ~ "A_monilatum_a",
    treat == "A_monilatum" ~ "A_monilatum_b",
    treat == "L4B9" ~ "L4-B9a",
    treat == "L4-B9" ~ "L4-B9b",
    TRUE ~ treat
  ))

Rho_assay_all <-
  Rho_assay_all %>% 
  convert_as_factor(treat, time) %>% 
  arrange(treat) %>% 
  mutate(conc_pg_ul2  = ifelse(!grepl("GD", treat), conc_pg_ul / 100, conc_pg_ul))
                                                                                 
Rho_assay_all_agg <-
  aggregate(data = Rho_assay_all, data ~ conc_pg_ul2 + treat + time, FUN = "mean")

# Calculate confidence intervals
pacman::p_load(broom)
confidence_intervals <- Rho_assay_all %>%
  drop_na(data) %>%
  group_by(conc_pg_ul2, time, treat) %>%
  do(tidy(t.test(.$data, conf.level = 0.95)))

# Create an empty list to store model predictions and models
model_predictions <- list()
model_Rho <- list()
EC50_Rho <- data.frame()

pacman::p_load(drc)

for (each_treat in unique(Rho_assay_all_agg$treat)) {
  for (each_time in unique(Rho_assay_all_agg$time)) {
    if (!grepl("GD", each_treat)) {
      if (nrow(Rho_assay_all_agg %>% filter(treat == each_treat &
                                            time == each_time)) > 0) {
        model <- drm(
          data = Rho_assay_all_agg %>% filter(treat == each_treat &
                                                time == each_time),
          data ~ conc_pg_ul2,
          fct = LL.4(),
          na.action = na.omit,
          lowerl = c(-Inf, 0, -Inf, -Inf),
          upperl = c(Inf, Inf, 100, Inf)
        )
        
        # Create a data frame for predictions and confidence intervals
        model_fits <-
          expand.grid(conc = c(seq(0, 1, length = 10000), seq(1, 50, length = 10000)))
        pm <-
          predict(model, newdata = model_fits, interval = "confidence")
        model_fits$p <- pm[, 1]
        model_fits$pmin <- pm[, 2]
        model_fits$pmax <- pm[, 3]
        model_fits$treat <- as.factor(each_treat)
        model_fits$time <- as.factor(each_time)
        # Store the model fits in the list
        model_predictions[[paste0(each_treat, each_time)]] <-
          model_fits
        # extract EC50 values and store in dataframe
        model_coef <-
          data.frame(
            treat = each_treat,
            time = each_time,
            EC50 = format(model$coefficients[[4]], digits = 2)
          )
        EC50_Rho <- rbind(EC50_Rho, model_coef)
      }
    }
  }
}

EC50_Rho <- unique(EC50_Rho)

Rho_fit_all <- do.call(rbind, model_predictions)

colors <- c(colorblind_pal()(6))

labels <-
  unique(levels(factor(
    Rho_assay_all_agg %>% filter(!grepl("GD", treat)) %>% pull(treat)
  )))

labels = c(
  "A_monilatum_a" = "<i>A. monilatum</i> I",
  "A_monilatum_b" = "<i>A. monilatum</i> II",
  "L2-D2" = "strain A",
  "L4-B1" = "strain B",
  "L4-B9a" = "strain C I",
  "L4-B9b" = "strain C II"
)

pacman::p_load(ggthemes, ggtext, ggplot2, ggtext)

Rho_all <-
  ggplot(data = confidence_intervals %>% filter(!grepl("GD", treat) & time == 24), aes(x = conc_pg_ul2, y = estimate, col = treat)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = T, size = 0.25) +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(lineheight = 1.5),
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown(),
    legend.text = element_markdown(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    legend.title = element_blank(),
    plot.subtitle = element_markdown()
  ) +
  geom_line(data = Rho_fit_all %>% filter(time == 24), aes(x = conc, y = p, col = treat), show.legend = F) +
  scale_y_continuous(
    limits = c(0, 125),
    labels = seq(0, 125, by = 25),
    breaks = seq(0, 125, by = 25)
  ) +
  ylab("Intact *R. salina* (% of control)") +
  scale_color_manual(labels = labels, values = colors) +
  scale_x_continuous(
    trans = 'log10',
    limits = c(0.1, 200),
    labels = scales::number_format()
  ) +  
  xlab("Goniodomin conc. (pg <i>&mu;</i>L<sup>-1</sup>)")  + guides(
    color = guide_legend(ncol = 3)) +
  labs(subtitle = "<b>a)</b>")

# for graphical abstract
Rho_all_ga <-
  ggplot(data = confidence_intervals %>% filter(!grepl("GD", treat) & time == 24), aes(x = conc_pg_ul2, y = estimate, col = treat)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = T, size = 0.25) +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(lineheight = 1.5),
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown(),
    legend.text = element_markdown(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    legend.title = element_blank(),
    plot.subtitle = element_text(),
    legend.box.margin = margin(-20, 0, 0, 0)
  ) +
  geom_line(data = Rho_fit_all %>% filter(time == 24), aes(x = conc, y = p, col = treat), show.legend = F) +
  scale_y_continuous(
    limits = c(0, 125),
    labels = seq(0, 125, by = 25),
    breaks = seq(0, 125, by = 25)
  ) +
  ylab("Intact *R. salina* (% of control)") +
  scale_color_manual(labels = labels, values = colors) +
  scale_x_continuous(
    trans = 'log10',
    limits = c(0.1, 100),
    labels = scales::number_format()
  ) +  
  xlab("Goniodomin conc. (pg <i>&mu;</i>L<sup>-1</sup>)")  + guides(
    color = guide_legend(ncol = 3)) +
  labs(subtitle = expression(bold("a+b" %->% "lysis / death")))

labels2 <- unique( Rho_assay_all_agg %>% filter(grepl("GD", treat)) %>% pull(treat))

Rho_all_GD <- ggplot(data = confidence_intervals %>% filter(grepl("GD", treat) & time == 24), aes(x = conc_pg_ul2, y = estimate, col = treat)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = T, size = 0.25, position = position_jitter(width = 0.1)) +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(lineheight = 1.5),
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown(),
    legend.text = element_text(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    legend.title = element_blank(),
    plot.subtitle = element_markdown()
  ) +
  scale_y_continuous(
    limits = c(0, 125),
    labels = seq(0, 125, by = 25),
    breaks = seq(0, 125, by = 25)
  ) +
  ylab("Intact *R. salina* (% of control)") +
  scale_color_manual(labels = labels2, values = colors) +
  scale_x_continuous(
    trans = 'log10',
    limits = c(0.1, 200),
    labels = scales::number_format()
  ) +  
  xlab("Goniodomin conc. (pg <i>&mu;</i>L<sup>-1</sup>)") + guides(
      color = guide_legend(ncol = 3)) +
  labs(subtitle = "<b>b)</b>")

# for graphical abstract
Rho_all_GD_ga <- ggplot(data = confidence_intervals %>% filter(grepl("GD", treat) & time == 24), aes(x = conc_pg_ul2, y = estimate, col = treat)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), show.legend = T, size = 0.25, position = position_jitter(width = 0.1)) +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(lineheight = 1.5),
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown(),
    legend.text = element_text(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    legend.title = element_blank(),
    plot.subtitle = element_text(),
    legend.box.margin = margin(-10, -10, 0, -10)
  ) +
  scale_y_continuous(
    limits = c(0, 125),
    labels = seq(0, 125, by = 25),
    breaks = seq(0, 125, by = 25)
  ) +
  ylab("Intact *R. salina* (% of control)") +
  scale_color_manual(labels = labels2, values = colors) +
  scale_x_continuous(
    trans = 'log10',
    limits = c(0.1, 100),
    labels = scales::number_format()
  ) +  
  xlab("Goniodomin conc. (pg <i>&mu;</i>L<sup>-1</sup>)") + guides(
    color = guide_legend(ncol = 3)) +
  labs(subtitle = expression(bold("b" %->% "no effect")))

pacman::p_load(ggpubr, patchwork)

Rho_all2 <- Rho_all + Rho_all_GD + plot_layout(ncol = 1, axis_titles = "collect") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))  

pacman::p_load(cowplot)

ggsave(
  "Rho_all_test.png",
  Rho_all2,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
  dpi = 300,
  width = 3.5,
  height = 5,
  units = "in"
)

# graphical abstract plot: 
graph_abstract <- P3a_ga + Rho_all_ga + P_GD_all_ga + Rho_all_GD_ga + plot_layout(ncol = 2, axis_titles = "collect_y") &
  theme(legend.position = "bottom",
        plot.margin = margin(t = 0, r = 20, b = 0, l = 20, unit = "pt"))


graph_abstract

ggsave(
  "graph_abstract.png",
  graph_abstract,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/figures",
  dpi = 300,
  width = 8,
  height = 5,
  units = "in"
)

# export word document with all packages used 

# install_github("Pakillo/grateful")
# cite_packages(out.format = "docx", out.dir = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP3/", out.file = "DRC_packages", pkgs = "Session")

# Garbage collection: call after large objects have been removed 
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up 
dev.off()

rm(list = ls())

.rs.restartR()
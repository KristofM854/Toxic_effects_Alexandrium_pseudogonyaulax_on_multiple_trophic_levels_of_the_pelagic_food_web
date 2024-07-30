##########################################
## Analysis of predator-prey interaction experiments of P. kofoidii and various Alexandrium species (A. pseudogonyaulax, A. catenella, A. limii)
## Published in: 
## All raw-data available on PANGAEA: 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 11.22-03.23
## Alfred-Wegener-Institute Bremerhaven
##########################################

# INSTALL AND LOAD PACKAGES ################################
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

##############
# Load data
pacman::p_load(readr)

data = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP2-Helgoland\\R-studio\\PK_mix.txt",
  delim = "\t",
  col_names = F
)

# Calculation of ingestion rates of PK/Alex5 or PK/A.t. and PK/A.p. ######
# Pk= Polykrikos kofoidii, Alex5 = A. catenella strain; A.t. = A. taylorii/A. limii, A.p. = A. pseudogonyaulax

# General data transformations
data <- as.data.frame(t(data[2:19, ]))
colnames(data) <- data[1, ]
data <- data[-1, ]

# introduce replicate and algae combination grouping variable
pacman::p_load(tidyr, dplyr)

data <- data %>%
  mutate(replicate = as.factor(c(rep(1:6, each = 1), rep(1:3, each = 1))),
         algae = as.factor(rep(c(
           "A5/Ap", "A5/Ap", "At/Ap"
         ), each = 3)))

# transform wide into long format
data1 <-
  data %>% pivot_longer(cols = c(-algae, -replicate),
                        names_to = "treat",
                        values_to = "data")

# extract data needed for calculation of ingestion rates 
ingestion_rates <-
  subset(
    data1,
    data1$treat == "Aufnahme_Ap_1h" |
      data1$treat == "Aufnahme_Alex5/At_1h" |
      data1$treat == "Aufnahme_Ap_3h" |
      data1$treat == "Aufnahme_Alex5/At_3h" |
      data1$treat == "Aufnahme_Ap_6h" |
      data1$treat == "Aufnahme_Alex5/At_6h"
  )

pacman::p_load(stringr)

# change data to numeric
ingestion_rates$data <- as.numeric(ingestion_rates$data)

# extract time and change to numeric
ingestion_rates$time <-
  as.numeric(str_sub(ingestion_rates$treat, -2, -2))

# combine algae and treat column into single column for statistical analysis
ingestion_rates <-
  ingestion_rates %>% arrange(treat) %>% mutate(treat2 = paste(algae, treat, sep = "-"))

ingestion_rates$treat2 <- str_sub(ingestion_rates$treat2, 0, -4)

# statistical analysis
# A. taylorii and A. pseudogonyaulax
pacman::p_load(rstatix)

stat1 <- ingestion_rates %>% filter(algae == "At/Ap")

stat1 %>% group_by(time) %>% kruskal_test(data ~ treat)

# A. taylorii and A. pseudogonyaulax
stat2 <- ingestion_rates %>% filter(algae == "A5/Ap")

stat2 %>% group_by(time) %>% kruskal_test(data ~ treat)

# # Data export for pangaea
write.table(
  ingestion_rates,
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/PANGAEA/PK_bialgal.txt",
  sep = "\t",
  row.names = FALSE
)

# calculate means and SD of ingestion rate export of ratios and plotting
ingestion_rates_mean <-
  aggregate(ingestion_rates, data ~ treat + algae, FUN = "mean")

# Calculate confidence intervals
pacman::p_load(broom)
confidence_intervals <- ingestion_rates %>%
  group_by(treat, algae) %>%
  do(tidy(t.test(.$data, conf.level = 0.95)))

confidence_intervals$time <- str_sub(confidence_intervals$treat, -2, -2)

confidence_intervals$a <- confidence_intervals$estimate - confidence_intervals$conf.low
confidence_intervals$b <- confidence_intervals$conf.high - confidence_intervals$estimate 

ingestion_rates_mean$sd <-
  aggregate(ingestion_rates, data ~ treat + algae, FUN = sd)$data

# # Reintroduce time variable
# ingestion_rates_mean$time <-
#   rep(c(1, 3), each = 1, length.out = 8)
# 
# # Introduce shape factor variable for plotting
# ingestion_rates_mean$shape <-
#   as.factor(rep(c(1:2), each = 2, length.out = 8))

# Prepare ggplot parameters
pacman::p_load(ggplot2, ggpubr)
limits <- aes(ymin = conf.low, ymax = conf.high)

dodge <- position_dodge(width = 0.6)

labels = expression(paste(atop(italic("A. catenella"), " non-lytic/PSP")),
                    paste(atop(
                      italic("A. pseudogonyaulax"), " lytic/GDA"
                    )),
                    paste(atop(italic("A. taylorii"), " GDA/lytic?")),
                    paste(atop(
                      italic("A. pseudogonyaulax"), " lytic/GDA"
                    )))

# calculate mean and SD of ingestion ratios of Alex5/Ap and At/Ap
ratios <- ingestion_rates_mean %>% 
  mutate(
    algae = rep(c("Alex5", "Ap", "At", "Ap"), each = 3), # prey item
    time = rep(c(1, 3, 6), each = 1, length.out = 12),
    exp = rep(c(1, 2), each = 6) # experiment grouping variable
  )

# Custom function to calculate ratios 
calculate_divided_data <- function(data) {
  divided_data <- data[1] / data[2]
  return(divided_data)
}

ratios$data <- as.numeric(sprintf("%.2f", ratios$data))
ratios$sd <- as.numeric(sprintf("%.2f", ratios$sd))

# Calculate ingestion ratios
ratios_results <- ratios %>%
  group_by(exp, time) %>%
  dplyr::summarise(
    combinations = combn(unique(algae), 2, paste, collapse = "/"),
    divided_data = calculate_divided_data(data),
    .groups = "drop"
  )

# Cell densities and ingestion rates of P. kofoidii feeding on monoalgal cultures of Alex5 or A. pseudogonyaulax #####
# Load data
data = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP2-Helgoland\\R-studio\\PK_sum.txt",
  delim = "\t",
  col_names = T
)

# introduce algae grouping factor
data$treat2 <-
  as.factor(rep(c("A5", "Ap", "Ap"), each = 3, length.out = nrow(data)))

# introduce time variable; each time point consists of 9 samples in total
data$time <- rep(c(3, 6, 9, 22.5, 27, 48, 72),
                 each = 9,
                 length.out = nrow(data))

# introduce time column of hunger control
data$time2 <-
  c(rep(c(3, 6, 9, 22.5, 27, 48, 72), each = 3), rep(NA, each = length(data$time) -
                                                       length(rep(
                                                         c(3, 6, 9, 22.5, 27, 48, 72), each = 3
                                                       ))))

# General data transformations for following statistical analysis
data1 <-
  data %>% drop_na(hunger) %>% dplyr::select(hunger, time2) %>% mutate(treat2 = "hunger")

colnames(data1) <- c("data", "time", "treat2")

data_aov <-
  full_join(data1,
            data %>% dplyr::select(data, treat2, time) %>% drop_na(data))

# introduce replicate column for repeated measures ANOVA and remove t = 9h, as there is only data for A. pseudogonyaulax
# and not for Alex 5 --> Error in repeated measures ANOVA
pacman::p_load(rstatix)
data_aov <- data_aov %>% arrange(treat2) %>%
  mutate(replicate = c(
    rep(1:3, length.out = 18),
    rep(1:6, length.out = 42),
    rep(1:3, length.out = 21)
  )) %>%
  convert_as_factor(time, treat2, replicate) %>%
  filter(time != 9)

# Statistical analysis #####
# log-transform cell counts first
data_aov$data2 <- log(data_aov$data)

# arrange time and introduce replicate grouping variable for rANOVA
data_aov <- data_aov %>% arrange(time)
data_aov$replicate2 <- rep(1:12, length.out = 72)

# # check distribution of the groups first:
# data_aov %>%
#   ggplot(aes(x = log(data), fill = treat2)) +
#   geom_histogram(
#     binwidth = 0.2,
#     position = "dodge",
#     color = "black",
#     alpha = 0.7
#   ) +
#   labs(title = "Histogram of Log-transformed Data by treat2",
#        x = "Log-transformed Data",
#        y = "Frequency") +
#   theme_minimal()

# perform repeated measures ANOVA
res.aov <-
  anova_test(
    data = data_aov,
    dv = data2,
    wid = replicate2,
    within = time,
    between = treat2
  )

# Display ANOVA table
get_anova_table(res.aov)

# Perform one-way ANOVAs and pairwise comparisons at each time point after significant interactions in repeated measures ANOVA
one.way <- data_aov %>%
  group_by(time) %>%
  anova_test(dv = data2, wid = replicate, within = treat2) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")

# Display adjusted p-values for one-way ANOVA
one.way

# Perform pairwise comparisons
pwc <- data_aov %>%
  group_by(time) %>%
  pairwise_t_test(data2 ~ treat2, paired = F,
                  p.adjust.method = "bonferroni") %>% filter(p.adj.signif != "ns" &
                                                               as.numeric(as.character(time)) > 47)

# Display pairwise comparisons
pwc

# # Data export for pangaea
# write.table(
#   data_aov,
#   "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/PANGAEA/PK_monoalgal.txt",
#   sep = "\t",
#   row.names = FALSE
# )

# Calculate P. kofoidii as cells / mL. Total volume was 2.5 mL in the experiments!
data$data <- data$data / 2.5
data$hunger <- data$hunger / 2.5

# calculate mean and sd of treatments and hunger control over time for plotting of mean Pk cell densities
data2 <- aggregate(data, data ~ treat2 + time, FUN = "mean")

# Calculate confidence intervals
pacman::p_load(broom)
confidence_intervals <- data %>%
  drop_na(data) %>%
  group_by(treat2, time) %>%
  do(tidy(t.test(.$data, conf.level = 0.95))) %>%
  rename(treat = treat2)

data2$sd <- aggregate(data, data ~ treat2 + time, FUN = "sd")$data

data3 <- aggregate(data, hunger ~ time2, FUN = "mean")

pacman::p_load(broom)
confidence_intervals2 <- data %>%
  drop_na(hunger) %>%
  group_by(time2) %>%
  do(tidy(t.test(.$hunger, conf.level = 0.95))) %>%
  rename(time = time2) %>%
  mutate(treat = "hunger")

data3$sd <- aggregate(data, hunger ~ time2, FUN = "sd")$hunger

data3$treat2 <- "hunger"
colnames(data3) <- c("time", "data", "sd", "treat2")

# full join both aggregated dataframes
data4 <- full_join(data2, data3)

confidence_intervals_all <- full_join(confidence_intervals, confidence_intervals2)

# prepare labels for ggplot figure
labels = c("<i>A. catenella</i>",
           "strain B",
           "hunger<br>control")

# colorblind colors but skip black (first color)
pacman::p_load(ggthemes, ggtext)
colors <- c(colorblind_pal()(4))[c(4, 2, 3)]

# prepare limits
limits <- aes(ymin = conf.low, ymax = conf.high)

dodge <- position_dodge(width = 0.6)

# Define the desired absolute width for the legend
absolute_legend_width <- 10

P1 <- ggplot(confidence_intervals_all, aes(x = time, y = estimate, col = treat)) +
  geom_pointrange(limits, position = position_jitter(width = 0.2), size = 0.5) +
  geom_line() +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = 'transparent'),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(1, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    axis.title.y = element_markdown(),
    legend.text = element_markdown(lineheight = 1.5),
    legend.margin = margin(r = absolute_legend_width),
    plot.subtitle = element_markdown()
  ) +
  scale_x_continuous(expression("Time (hours)")) +
  scale_color_manual(labels = labels,
                     values = colors) +
  scale_y_continuous("Cell density <i>P. kofoidii<br></i>(cells mL<sup>-1</sup>)")+
  labs(subtitle = "<b>a)</b>")

# # export figure
# ggsave(
#   "Mean_sum_PK.png",
#   P1,
#   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/AP2",
#   dpi = 300,
#   width = 18,
#   height = 10,
#   units = "cm"
# )

# Calculate ingestion rates of P. kofoidii feeding on Alex 5 or A. pseudogonyaulax
data$ingestion_rate <- data$Futterzellen_pro_P.k. / data$time

# calculate mean and SD of ingestion rates
pacman::p_load(broom)
confidence_intervals3 <- data[1:51, ] %>%
  drop_na(ingestion_rate) %>%
  group_by(time, treat2) %>%
  do(tidy(t.test(.$ingestion_rate, conf.level = 0.95))) %>%
  rename(treat = treat2) 

# introduce replicate column for repeated measures ANOVA
data <-
  data %>% arrange(treat2) %>% mutate(replicate = c(rep(1:3, length.out =
                                                          21), rep(1:6, length.out = 42)))
# # check distribution of the groups first:
# data_aov2 %>%
#   ggplot(aes(x = log(ingestion_rate), fill = treat2)) +
#   geom_histogram(
#     binwidth = 0.2,
#     position = "dodge",
#     color = "black",
#     alpha = 0.7
#   ) +
#   labs(title = "Histogram of Log-transformed Data by treat2",
#        x = "Log-transformed Data",
#        y = "Frequency") +
#   theme_minimal()

# Repeated measures ANOVA of ingestion rates
data_aov2 <-
  data %>% drop_na(ingestion_rate) %>% filter(time != 9) %>% convert_as_factor(replicate, treat2, time)

# # export ingestion rates as table
# data5$data <- paste(data5$ingestion_rate, data5$sd, sep = " +/- ")
# data5$treat2 <-
#   ifelse(data5$treat2 == "A5",
#          "Alex 5",
#          ifelse(data5$treat2 == "Ap", "A. pseudogonyaulax", NA))
# 
# tab <-
#   data5 %>% arrange(time) %>% flextable() %>% autofit() %>% save_as_docx(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/Ingestion_rate_PK_A5_AP.docx")

# add small constant to all ingestion rates to treat zeros for log-transformation
data_aov2$ingestion_rate <- data_aov2$ingestion_rate + 0.0001
data_aov2$ingestion_rate2 <- log(data_aov2$ingestion_rate)

data_aov2 <-
  data_aov2 %>% arrange(time) %>% mutate(replicate2 = rep(1:9, length.out = length(data_aov2$replicate)))

# perform repeated measures ANOVA
res.aov2 <-
  anova_test(
    data = data_aov2,
    dv = ingestion_rate2,
    wid = replicate2,
    within = time,
    between = treat2
  )

# Display results of rANOVA
get_anova_table(res.aov2)

# Perform one-way ANOVAs and pairwise comparisons at each time point after significant interactions in repeated measures ANOVA
one.way2 <- data_aov2 %>%
  group_by(time) %>%
  anova_test(dv = ingestion_rate2, wid = replicate, within = treat2) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")

# Display adjusted p-values for one-way ANOVA
one.way2

# Perform pairwise comparisons
pwc2 <- data_aov2 %>%
  group_by(time) %>%
  pairwise_t_test(ingestion_rate ~ treat2,
                  paired = F,
                  p.adjust.method = "BH")

# Display pairwise comparisons
pwc2

# # Data export for pangaea
# write.table(
#   data_aov2,
#   "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/PANGAEA/PK_monoalgal_ingestion.txt",
#   sep = "\t",
#   row.names = FALSE
# )

# prepare ggplot parameters
limits <- aes(ymin = conf.low, ymax = conf.high)

P2 <-
  ggplot(confidence_intervals3,
         aes(
           x = time,
           y = estimate,
           col = treat
         )) +
  geom_pointrange(limits, size = 0.5) +
  geom_line() +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.text.align = 0.5,
    legend.text = element_markdown(lineheight = 1.5),
    axis.title.y = element_markdown(),
    legend.margin = margin(r = absolute_legend_width),
    plot.subtitle = element_markdown()
  ) +
  ylab("Ingestion rates <br><i>P. kofoidii </i>(cells  h<sup>-1</sup>)") +
  xlab(expression("Time (hours)")) +
  scale_color_manual(
    labels = c("<i>A. catenella</i>", "strain B"),
    values = colors[1:2]
  ) +
  labs(subtitle = "<b>b)</b>")

# combine both plots and export as figure
pacman::p_load(patchwork)
g <- P1 + P2 + plot_layout(ncol = 1, axis_titles = "collect_x")

ggsave(
  "PK_1.png",
  g,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/AP2",
  dpi = 300,
  width = 3.5,
  height = 4.5,
  units = "in"
)

# export word document with package citations used in ingestion rate and GDA R-files
# pacman::p_load(grateful)
# cite_packages(
#   out.format = "docx",
#   out.dir = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/",
#   pkgs = "Session",
#   out.file = "PK_packages"
# )

# Garbage collection: call after large objects have been removed
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()
rm(list = ls())

# restart R
.rs.restartR()

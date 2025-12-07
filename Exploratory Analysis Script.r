#  Setup & Load 
library(tidyverse)
library(ggplot2)
library(corrplot)
library(naniar)      
library(reshape2)

# Open Source Dataset via Kaggle: 
# https://www.kaggle.com/datasets/shreyasvedpathak/pcos-dataset 
pcos <- read_csv("PCOS_data.csv")

# Quick checks
glimpse(pcos)
summary(pcos)

#  Missing Values Check 
missing_vals <- pcos %>% summarise(across(everything(), ~ sum(is.na(.))))

# Visual missingness
x <- skim(pcos)
x <- x %>% select(skim_variable, numeric.hist)

vis_miss(pcos) +
  ggtitle("Missing Values by Variable and Observation")

missing_df <- pcos %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0) %>%           
  arrange(n_missing)                 

ggplot(missing_df, aes(x = n_missing, y = variable)) +
  geom_col(fill = "#C9941D") +       
  scale_x_continuous(expand = c(0,0)) +
  labs(
    title = "Missing Values by Column",
    x = "Count of Missing",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),     
    panel.grid.major.x = element_line(      
      color = "grey70",
      linetype = "dotted"
    ),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 13)
  )


#  Data Types & Conversions
pcos <- pcos %>%
  rename(pcos = `PCOS (Y/N)`) 

continuous_vars <- pcos %>%
  select(where(is.numeric))

summary(continuous_vars)

# Compare PCOS vs non-PCOS groups for key variables
ggplot(pcos, aes(x = factor(pcos), fill = factor(pcos))) +
  geom_bar() +
  scale_fill_manual(
    values = c(
      "0" = "#C9941D",  # gold
      "1" = "#68AFD7"   # blue
    )
  ) +
  labs(
    title = "Class Balance: PCOS vs Non-PCOS",
    x = "PCOS (Y/N)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    
    # match dashed gridlines
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey80",
      linetype = "dashed"
    ),
    
    plot.title = element_text(size = 18, face = "bold")
  )

pcos %>%
  group_by(pcos) %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, median = median, sd = sd),
                   na.rm = TRUE))

#  Histograms / Distributions 
ggplot(pcos, aes(x = BMI)) +
  geom_histogram(bins = 30, fill = "gold", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of BMI", x = "BMI", y = "Count")

# Outlier Detection / Graphics 
pcos_num <- pcos %>% select(where(is.numeric))

# Compute z-scores for each numeric variable for outliers
pcos_z <- pcos_num %>% 
  mutate(across(everything(), ~ as.numeric(scale(.))))

pcos_outliers <- pcos_z %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(
    cols = -row_id,
    names_to = "variable",
    values_to = "z"
  ) %>%
  mutate(
    is_outlier = abs(z) > 3
  )

outlier_summary <- pcos_outliers %>%
  group_by(variable) %>%
  summarise(
    n_outliers = sum(is_outlier, na.rm = TRUE),
    pct_outliers = mean(is_outlier, na.rm = TRUE) * 100
  ) %>%
  arrange(desc(n_outliers))

# Boxplot of BMI by PCOS status
ggplot(pcos, aes(x = factor(pcos), y = BMI, fill = factor(pcos))) +
  geom_boxplot(
    color = "grey30",
    outlier.shape = 16,
    outlier.size = 2
  ) +
  scale_fill_manual(
    values = c(
      "0" = "#C9941D",  # gold
      "1" = "#68AFD7"   # blue
    )
  ) +
  labs(
    title = "BMI by PCOS Status",
    x = "PCOS (Y/N)",
    y = "BMI"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    
    # match dotted horizontal gridlines
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey80",
      linetype = "dashed"
    ),
    
    plot.title = element_text(size = 18, face = "bold")
  )



#  Correlation Matrix (Continuous Features) 
corr_mat <- cor(continuous_vars, use = "pairwise.complete.obs")
corrplot(corr_mat, method = "color", type = "upper", tl.cex = 0.7)

# Melt correlation for easier filtering 
melted_corr <- melt(corr_mat)
melted_corr %>%
  filter(abs(value) > 0.6 & Var1 != Var2) %>%
  arrange(desc(abs(value)))

#  Initial Insights & Flagged Variables 
if (all(c("LH", "FSH") %in% names(pcos))) {
  pcos <- pcos %>%
    mutate(lh_fsh_ratio = LH / FSH)
  ggplot(pcos, aes(x = pcos, y = lh_fsh_ratio)) +
    geom_boxplot(fill = "lightgreen") +
    labs(title = "LH/FSH Ratio by PCOS Status.")
}

# Differences in means test (t-test) for BMI
t.test(BMI ~ pcos, data = pcos)



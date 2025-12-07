# Libraries loaded
library(tidyverse)
library(vip)
library(xgboost)
library(cluster)
library(factoextra)
library(corrplot)
library(skimr)

set.seed(123)

# Load & clean data 

#If need to reload data:
# pcos <- read_csv("PCOS_data.csv")

pcos <- pcos %>%
  select(-c(`Sl. No`, `Patient File No.`, `...45`)) %>%
  mutate(`PCOS (Y/N)` = factor(`pcos`, levels=c(0,1), labels=c("No","Yes"))) %>%
  select(-pcos)

skim(pcos)

# Descriptive Stats & Visualizations
summary(select(pcos, where(is.numeric)))

ggplot(pcos, aes(x=BMI, fill=`PCOS (Y/N)`)) +
  geom_histogram(bins=30, alpha=0.6, position="identity")

pcos_num <- pcos %>% select(where(is.numeric))
corr <- cor(pcos_num, use="pairwise.complete.obs")
corrplot(corr, type="lower", tl.cex=.6)


pcos_num_measurements <- pcos %>% select(BMI, `Waist(inch)`, `Hip(inch)`, `Weight (Kg)`)
corr2 <- cor(pcos_num_measurements, use="pairwise.complete.obs")
corrplot(corr2, type="lower", tl.cex=.6)

# Feature Creation
# BMI CATs derived via CDC:
# Source: https://www.cdc.gov/bmi/adult-calculator/bmi-categories.html
pcos <- pcos %>%
  mutate(BMI_cat = case_when(
    BMI < 18.5 ~ "Underweight",
    BMI < 25   ~ "Normal",
    BMI < 30   ~ "Overweight",
    TRUE       ~ "Obese"
  ))

# PCA : Dimension Reduction
pcos_num <- pcos %>%
  select(where(is.numeric)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm=TRUE), .)))

pcos_scaled <- scale(pcos_num)
pc <- prcomp(pcos_scaled, center=TRUE, scale.=TRUE)

# scree plot
fviz_eig(pc)      
# biplot
fviz_pca_biplot(pc, repel=TRUE)

cum_var <- summary(pc)$importance[3,]
keep <- which(cum_var >= 0.50)[1]
pcos_pca <- as.data.frame(pc$x[,1:keep])
pcos_pca$PCOS <- pcos$`PCOS (Y/N)`

# Clustering on PCA 
set.seed(123)
km <- kmeans(pcos_pca %>% select(-PCOS), centers=3, nstart=20)
fviz_cluster(km, data=pcos_pca %>% select(-PCOS))

pcos$cluster <- factor(km$cluster)
table(pcos$cluster, pcos$pcos)

pcos %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric),
                   list(mean = mean, median = median, sd = sd),
                   na.rm = TRUE))


# Logistic Regression 
logit_model <- glm(`PCOS (Y/N)` ~ ., data=pcos %>% select(-BMI_cat), family="binomial")
summary(logit_model)
logit_summary <- broom::tidy(logit_model)
write.csv(logit_summary, "plots/logit_summary.csv", row.names=FALSE)

# XGBoost
target <- as.numeric(pcos$`PCOS (Y/N)`) - 1
dtrain <- xgb.DMatrix(data=as.matrix(pcos_num), label=target)

xgb_model <- xgboost(data=dtrain,
                     max_depth=3,
                     nrounds=50,
                     objective="binary:logistic",
                     verbose=0)

vip(xgb_model, num_features=10)


## PLOTS SECTION: 
# Create a folder for plots
dir.create("plots", showWarnings = FALSE)

# Descriptive Stats Plots
library(ggplot2)
p1 <- ggplot(pcos, aes(x=BMI, fill=`PCOS (Y/N)`)) +
  geom_histogram(bins=30, alpha=0.6, position="identity") +
  labs(title="BMI Distribution by PCOS status")
ggsave("plots/bmi_distribution.png", p1, width=6, height=4)

p2 <- ggplot(pcos, aes(x=`Age (yrs)`, y=BMI, color=`PCOS (Y/N)`)) +
  geom_point(alpha=0.6) +
  labs(title="Age vs BMI by PCOS status")
ggsave("plots/age_bmi_scatter.png", p2, width=6, height=4)

# Correlation heatmap
pcos_num <- pcos %>% select(where(is.numeric))
corr <- cor(pcos_num, use="pairwise.complete.obs")
png("plots/correlation_heatmap.png", width=800, height=600)
corrplot(corr, type="lower", tl.cex=.6)
dev.off()

# PCA Plots
fviz_eig(pc)
ggsave("plots/pca_scree.png", width=6, height=4)

fviz_pca_biplot(pc, repel=TRUE)
ggsave("plots/pca_biplot.png", width=6, height=4)

# Clustering Plots
factoextra::fviz_eig(pc)
factoextra::fviz_pca_biplot(pc, repel=TRUE)
#fviz_cluster(km, data=pcos_pca %>% select(-PCOS)) #Run prior in script
ggsave("plots/cluster_plot.png", width=6, height=4)

# XGBoost Plot
png("plots/xgb_vip.png", width=800, height=600)
vip(xgb_model, num_features=10)
dev.off()


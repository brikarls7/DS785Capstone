# DS785Capstone
Detection of PCOS Through Clinical and Lifestyle Data Using Machine Learning Techniques

This repository serves as my codebase for the work completed for my capstone project. The three main scripts were developed for the stages in which the project was completed: Initial Analysis, Feature Engineering, and Final Modeling. 

Project Abstract:
Polycystic Ovary Syndrome (PCOS) is one of the most prevalent endocrine disorders among women of reproductive age, yet it remains severely under diagnosed. This project applies machine learning techniques to predict PCOS using a combination of clinical, lifestyle and hormonal features. The primary objective is to evaluate and compare the performance of three models—Logistic Regression, K-Nearest Neighbors (KNN), and XGBoost—in identifying patients at risk for PCOS. Data preprocessing included removal of non-predictive identifiers, median imputation for missing values, stratified train-test splitting, and feature engineering to create BMI categories and LH/FSH ratios for clinical interpretability. Each model was trained and validated using five-fold cross-validation, optimizing for ROC AUC as the primary performance metric. Results indicate that hormonal markers such as Follicle Numbers, AMH and the LH/FSH ratio were the most influential predictors across models. The findings highlight the potential for integrating machine learning into women’s health diagnostics, offering a scalable and data-driven approach to early PCOS detection. 


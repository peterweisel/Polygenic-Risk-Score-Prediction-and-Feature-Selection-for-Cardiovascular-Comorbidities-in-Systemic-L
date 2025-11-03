#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Peter Joseph Weisel
Affiliation: Forskargruppen reumatologi @ Institutionen för medicinska vetenskaper
Description: Lasso regression for feature selection
Date: 2025-10-24
Developed using https://medium.com/@agrawalsam1997/feature-selection-using-lasso-regression-10f49c973f08
"""

# load libraries
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split, GridSearchCV, KFold


def run_lasso_feature_selection(
    # perform feature selection on rsID/patient matrix using Lass regression
    # input matrix:
    input_csv: str,
    output_dir: str,
    # cardiovascular comorbidity:
    target_col: str = "PE",
    # define for GridSearchCV:
    alpha_min: float = 0.00001,
    alpha_max: float = 10,
    alpha_step: float = 500,
    # determines number of important features:
    importance_threshold: float = 0.001
):

    print(f"\n Processing: {os.path.basename(input_csv)}")
    df = pd.read_csv(input_csv)
    print(f"Shape: {df.shape}")

    # fill NaNs with 0 effect alleles
    df = df.fillna(0)

    # target column (cardiovascular comorbidity)
    if target_col not in df.columns:
        print(f"Using last column: {df.columns[-1]}")
        target_col = df.columns[-1]

    # drop target (only keep SNPs)
    X = df.drop(target_col, axis=1).values
    y = df[target_col].values
    # split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # select best hyperparameter
    params = {"alpha": np.arange(alpha_min, alpha_max, alpha_step)}
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    lasso = Lasso(max_iter=10000)

    lasso_cv = GridSearchCV(lasso, param_grid=params, cv=kf)
    lasso_cv.fit(X_train, y_train)
    best_alpha = lasso_cv.best_params_["alpha"]
    print(f"Best alpha: {best_alpha}")

    # final Lasso model trained on subset
    lasso_final = Lasso(alpha=best_alpha, max_iter=10000)
    lasso_final.fit(X_train, y_train)
    coef = np.abs(lasso_final.coef_)
    names = df.drop(target_col, axis=1).columns

    # plot (plot feature importance)
    plt.figure(figsize=(12, 6))
    plt.bar(names, coef)
    plt.xticks(rotation=90)
    plt.grid(True)
    plt.title(f"Lasso Feature Importance – {os.path.basename(input_csv)}")
    plt.xlabel("Features")
    plt.ylabel("Importance")
    plt.tight_layout()

    # save image to Desktop
    basename = os.path.splitext(os.path.basename(input_csv))[0]
    plot_path = os.path.join(output_dir, f"{basename}_feature_importance.png")
    plt.savefig(plot_path)
    plt.close()
    print(f"Plot saved: {plot_path}")

    # save important features (> threshold) to Desktop
    feature_importance = pd.DataFrame({"Feature": names, "Importance": coef})
    feature_importance = (
        feature_importance[feature_importance["Importance"] > importance_threshold]
        .sort_values(by="Importance", ascending=False)
        .reset_index(drop=True)
    )

    csv_path = os.path.join(output_dir, f"{basename}_feature_importance.csv")
    feature_importance.to_csv(csv_path, index=False)
    print(f"CSV saved: {csv_path}")

# run on all cardiovascular comorbidities
if __name__ == "__main__":
    base_dir = "/Users/peterweisel/Library/CloudStorage/OneDrive-Uppsalauniversitet/PWeisel_Research_Training/_plink_dataset/feature_selection/directionality"
    output_dir = os.path.join(base_dir, "female_cohort")

    os.makedirs(output_dir, exist_ok=True)

    files = [
        "MI_direction_df_f.csv",
        "PE_direction_df_f.csv",
        "AP_direction_df_f.csv",
        "DVT_direction_df_f.csv",
        "TIA_direction_df_f.csv",
        "ICVL_direction_df_f.csv",
    ]

    for file in files:
        input_csv = os.path.join(base_dir, file)
        run_lasso_feature_selection(input_csv, output_dir)

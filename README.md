# Spatial Divides of Wealth Inequality and Politics in Switzerland
## A Kernel and Optimal Transport Approach

This repository reproduces the results of the paper:
**"Spatial Divides of Wealth Inequality and Politics in Switzerland: A Kernel and Optimal Transport Approach".**

We reconstruct municipal income distributions from Swiss federal tax classes,
build economic/political/linguistic/spatial kernels, compare them at the
**kernel level** (RV, z-scores) and at the **factor level** (correlations between first MDS axes),
and evaluate explanatory power via weighted regressions and residual maps.

### In detail, what is done here
This code investigates how economic, political, spatial, and linguistic factors jointly structure the political geography of Switzerland. Leveraging detailed municipal-level income tax data, we reconstruct wealth distributions and compute different inter-municipality economic distances from Swiss Federal Tax data. These, declined into simple mean and optimal transport distances (Wasserstein), are compared to political distances derived from 381 national referenda (1971–2024). Using a kernel-based approach, multidimensional scaling (MDS), spatial autocorrelation index and multiple regression, we explore the alignment between economic inequality, political behavior, and spatial-linguistic structures. Our results reveal that optimal transport income disparities explain political divides better than do the average wealth, while linguistic and urban–rural
divisions remain largely dominant structuring forces of political explanation. This study highlights the potential of MDS visualization and the kernel associations in analyzing spatially embedded political behavior.

## 🧩 Code organisation
- ```01_download_data```
  - "detailed_language_2024": id of each municipality detail of spoken languages
  - "distance_mat_2024": squared matrix of road distances between municipalities in meters
  - "IFD_tax_classes_2024": 10 classes of swiss taxes in CHF per municipality
  - "IFD_tax_wealth_2024": 10 classes of swiss taxes payers per municipality
  - "time_mat_2024": squared matrix of road time between municipalities in seconds
  - "vote_info_2024": name of each vote object
  - "vote_nb_valid_2024": valid ballot
  - "vote_nb_yes_2024":  number of 'yes'
  - "vote_theme_names": theme names of votes
  - "vote_yes_2024": pourcentage of 'yes' per municipality
  - "f": vector of weight of each municipality
- ```02_build_kernels```
  - Build dissimilarity matrices (D) and f-centered kernels (K)
- ```03_produce_results```
  - Compute kernel-level (RV, z) and factor-level (axis-1) associations, partial correlations, and weighted regressions.
  - Export tables & figures
- ```RV_funs```
  - function that produces
    - "Y_list":         MDS factors
    - "Delta_list":     Delta values
    - "RVi":            RV value between kernels
    - "E_RV":           Expected value between kernels
    - "Var_RV"          Variance between kernels
    - "Z_RV"            Z-score between kernels
    - "eigen_val_list": Eigenvalues et eigenvectors of each kernel
    - "nb_lambda_pos":  Number of positive Lambda values per kernel

⚠️🚧 GitHub repo in construction...

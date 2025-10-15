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

⚠️🚧 GitHub repo in construction...

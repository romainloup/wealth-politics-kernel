# Spatial Divides of Wealth Inequality and Politics in Switzerland
_A Kernel and Optimal Transport Approach_

This repository reproduces the results of the paper:
**“Spatial Divides of Wealth Inequality and Politics in Switzerland: A Kernel and Optimal Transport Approach.”**

We reconstruct municipal income distributions from Swiss federal tax classes,
build economic/political/linguistic/spatial kernels, compare them at the
**kernel level** (RV, z-scores) and at the **factor level** (correlations between first MDS axes),
and evaluate explanatory power via weighted regressions and residual maps.

---

## 🧩 What’s inside
- `R/` – analysis pipeline (data ingest → kernels → MDS → associations → regressions → figures/tables)
- `data/` – raw (or scripts), harmonized, and processed data
- `figs/` – all figures (subfolders) exported by the pipeline
- `tables/` – LaTeX tables used in the paper
- `_targets.R` – pipeline definition (targets)
- `params/` – configuration (income bins, midpoints, options)
- `renv.lock` – frozen R package versions for reproducibility

---

## 🔧 Reproducibility

### 1) Restore the R environment
```r
install.packages("renv")
renv::restore()

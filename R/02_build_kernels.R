# =============================================================================
# 02_build_kernels.R
# Build dissimilarity matrices (D) and f-centered kernels (K)
# =============================================================================

source(here::here("R", "RV_funs.R"))  # si tu y ranges des helpers utiles
# (tu peux aussi avoir un "kernels_funs.R" dédié au centrage f, MDS, etc.)

# ---- Helpers noyaux ----
f_center = function(K, f) {
  n = length(f)
  H = diag(n) - matrix(1, n, 1) %*% t(matrix(f, n, 1))
  diag_sqrt_f = diag(sqrt(f), n)
  diag_sqrt_f %*% H %*% K %*% t(H) %*% diag_sqrt_f
}

K_from_D = function(D, f) {
  # K = -0.5 * Π^{1/2} H D H Π^{1/2}
  -0.5 * f_center(D, f)
}

# ---- Fonction principale ----
build_kernels = function(inputs) {
  dl = inputs$data_list
  
  # 1) poids f
  f = if ("f" %in% names(dl)) {
    # soit une colonne unique, soit une table avec colonne 'f'
    vec = dl$f[[1]]
    if (is.null(vec)) stop("f.csv must contain a column with weights.")
    as.numeric(vec) / sum(vec)
  } else stop("Missing f in data_list.")
  
  n = length(f)
  
  # 2) Political distances (squared Euclidean on %yes)
  votes_mat = as.matrix(dl$vote_yes_2024[, -(1:2), drop = FALSE])
  D_pol     = as.matrix(stats::dist(votes_mat))^2
  K_pol     = K_from_D(D_pol, f)
  
  # 3) Road time (symmetrized, squared)
  Traw  = as.matrix(dl$time_mat_2024)
  D_time = ((Traw + t(Traw)) / 2)^2
  diag(D_time) = 0
  K_time = K_from_D(D_time, f)
  
  # 4) Size (distance on log weights)
  D_size = as.matrix(stats::dist(log(f)))^2
  K_size = K_from_D(D_size, f)
  
  # 5) Language (one-hot, categorical kernel)
  Z_pre = tibble::tibble(
    municipality = seq_len(nrow(dl$detailed_language_2024)),
    language     = as.factor(dl$detailed_language_2024$language_region)
  ) |>
    tidyr::pivot_wider(
      names_from = language, values_from = municipality,
      values_fill = 0, values_fn = length, names_prefix = "lang"
    )
  Z = as.matrix(Z_pre[, -1, drop = FALSE])
  rho = as.vector(t(Z) %*% f)  # masses par langue
  K_ling = diag(sqrt(f), n) %*% (Z %*% diag(1 / rho) %*% t(Z) - 1) %*% diag(sqrt(f), n)
  
  # 6) Wealth (mean-based L1/L2)
  # NOTE: adapte selon la signification de tes CSV (montants vs comptes)
  wealth_mat = as.matrix(dl$IFD_tax_classes_2024[, -c(1, 2, 13)])
  payer_mat    = as.matrix(dl$IFD_tax_wealth_2024[, -c(1, 2, 13)])
  
  # moyenne par contribuable (si wealth_mat = montants totaux ; payer_mat = nb contribuables)
  x_income = rowSums(wealth_mat) / pmax(rowSums(payer_mat), 1)
  
  D_wealth1 = abs(outer(x_income, x_income, `-`))
  D_wealth2 = (outer(x_income, x_income, `-`))^2
  diag(D_wealth1) = 0; diag(D_wealth2) = 0
  
  K_wealth1 = K_from_D(D_wealth1, f)
  K_wealth2 = K_from_D(D_wealth2, f)
  
  # 7) OT distances (on full distributions)
  # g : proportions par classe (utilise les comptes, pas les montants)
  G_counts = as.matrix(dl$IFD_tax_wealth_2024[, -c(1, 2, 13)])
  row_sums = pmax(rowSums(G_counts), 1)
  g = G_counts / row_sums
  
  # midpoints (kCHF) — adapte si tu stockes ça dans params/
  mid_y = c(0, 15000, 35000, 45000, 62500, 87500, 150000, 350000, 750000, 2000000)
  # CDF par commune
  Fg = t(apply(g, 1, cumsum))
  
  # longueurs d’intervalle (Δy)
  dy = diff(mid_y)
  # mat [n x (m-1)] en recyclant dy
  W1 = matrix(0, n, n)
  W2 = matrix(0, n, n)
  
  # version vectorisée efficace
  for (i in seq_len(n)) {
    # diff CDF vs toutes les autres communes
    diff_cdf = sweep(Fg, 2, Fg[i, ], `-`)
    # retirer la dernière classe (OT se calcule sur les intervalles)
    diff_cdf_i = diff_cdf[, -ncol(diff_cdf), drop = FALSE]
    W1[i, ] = colSums(t(abs(diff_cdf_i)) * dy)
    W2[i, ] = colSums(t(diff_cdf_i^2) * dy)
  }
  # symétriser
  D_OT1 = (W1 + t(W1)) / 2; diag(D_OT1) = 0
  D_OT2 = (W2 + t(W2)) / 2; diag(D_OT2) = 0
  
  K_OT1 = K_from_D(D_OT1, f)
  K_OT2 = K_from_D(D_OT2, f)
  
  # 8) Retour compact
  list(
    f = f,
    D = list(pol = D_pol, time = D_time, size = D_size, ling = NULL, wealth1 = D_wealth1,
             wealth2 = D_wealth2, OT1 = D_OT1, OT2 = D_OT2),
    K = list(pol = K_pol, time = K_time, size = K_size, ling = K_ling,
             wealth1 = K_wealth1, wealth2 = K_wealth2, OT1 = K_OT1, OT2 = K_OT2)
  )
}

# ---- Exécution (si appelé en script) ----
if (sys.nframe() == 0L) {
  inputs  = readRDS(here::here("data", "interim", "inputs.rds"))
  kernels = build_kernels(inputs)
  saveRDS(kernels, here::here("data", "processed", "kernels.rds"))
  message("Kernels saved to data/processed/kernels.rds")
}



###############################################################################
# =============================================================================
# 02_build_kernels.R
# Author : Romain Loup
# Object : Build dissimilarity matrices (D) and f-centered kernels (K)
# Date   : 2025-09-18
# =============================================================================

# Number of municipality
n = length(f)
# --- Centering matrix
H = diag(n) - rep(1, n) %*% t(f) # centering matrix

# --- Political
D_pol = as.matrix(dist(data_list$vote_yes_2024[,-c(1:2)])^2) # political distances between municipalities
K_pol = -0.5 * diag(sqrt(f)) %*% H %*% D_pol %*% t(H) %*% diag(sqrt(f)) # political kernel

# --- Road time
D_time = as.matrix((data_list$time_mat_2024 + t(data_list$time_mat_2024)) / 2) ^ 2 # road time
K_time = -0.5 * diag(sqrt(f)) %*% H %*% D_time %*% t(H) %*% diag(sqrt(f)) # spatial kernel, time

# --- Weight
# kernel with weights themselves (^2)
D_size = as.matrix(dist(log(f)) ^ 2)
K_size = -0.5 * diag(sqrt(f)) %*% H %*% D_size %*% t(H) %*% diag(sqrt(f))

# --- Linguistic
Z_pre = data.frame(1:n, as.factor(data_list$detailed_language_2024$language_region))
names(Z_pre) = c("municipality", "language")
#pivot wider 
library(tidyr)
Z_pre = Z_pre %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = language, values_from = value, values_fill = list(value = 0), names_prefix = "language")
Z = as.matrix(Z_pre[,-1])

# --- Linguistic kernel
rho = t(Z)%*%f
K_ling = diag(sqrt(f))%*%(Z%*%diag(1/as.vector(rho))%*%t(Z)-1)%*%diag(sqrt(f))

# --- Wealth kernel
wealth_mat = as.matrix(data_list$IFD_tax_classes_2024[, -c(1,2,13)])
payer_mat = as.matrix(data_list$IFD_tax_wealth_2024[, -c(1,2,13)])

richesse_vect = rowSums(wealth_mat)/rowSums(payer_mat)

# L1 Distance : d_{ij}^{wealth_1} = |x_i - x_j|
D_wealth1 = abs(outer(richesse_vect, richesse_vect, "-"))

# L2 Distance : d_{ij}^{wealth_2} = (x_i - x_j)^2
D_wealth2 = (outer(richesse_vect, richesse_vect, "-"))^2

# Convert to proportions (each line = distribution)
row_sums = rowSums(wealth_mat)
g = wealth_mat / row_sums

# Define the tax payment class centers
mid_y = c(0, 15000, 35000, 45000, 62500, 87500, 150000, 350000, 750000, 2000000)

# Cumulative function
cumsum_mat = t(apply(g, 1, cumsum))

# Initialize the distance matrices
D_OT1 = matrix(0, n, n)
D_OT2 = matrix(0, n, n)

# Main loop
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    # Difference in distribution functions
    diff_cdf = cumsum_mat[i, ] - cumsum_mat[j, ]
    
    # Wasserstein-1
    d1 = sum(abs(diff_cdf) * c(diff(mid_y), 0))  # Note: diff(mid_y) = lengths between classes
    # Wasserstein-2
    d2 = sum((diff_cdf^2) * c(diff(mid_y), 0))
    
    D_OT1[i, j] = D_OT1[j, i] = d1
    D_OT2[i, j] = D_OT2[j, i] = d2
  }
  print(paste("i",i))
}

# --- Matrix calculation version
# Assume that `cumsum_mat` is already defined as follows:
# cumsum_mat = t(apply(g, 1, cumsum))

# Define the lengths between classes (delta y)
dy = c(diff(y), 0)  # On ajoute un zéro pour que dy ait longueur 10

# Center the weights according to dy
weighted_diff = function(cdf_mat, dy, power = 1) {
  # (n x 10 x n)  implicit table of all pairs of lines
  n = nrow(cdf_mat)
  # Paired matrix calculations (via broadcasting)
  # => We extend the matrices (via replication) o make paired differences
  
  # Extend cdf (3D implicit)
  A = array(rep(t(cdf_mat), times = n), dim = c(10, n, n))     # A[ , i, j] = cdf[i, ]
  B = array(rep(cdf_mat, each = n), dim = c(n, 10, n))         # B[ i, , j] = cdf[j, ]
  
  diff_array = A - aperm(B, c(2,1,3))  # Dimensions: [10, n, n]
  diff_array = aperm(diff_array, c(2,1,3))  # Dimensions: [n, 10, n]
  
  # Apply the dy weighting and the exponent
  dy_mat = matrix(dy, nrow = 1)
  if (power == 1) {
    w_dist = apply(abs(diff_array), c(1,3), function(x) sum(x * dy))
  } else if (power == 2) {
    w_dist = apply(diff_array^2, c(1,3), function(x) sum(x * dy))
  }
  return(w_dist)
}

# Compute distance matrices
D_OT1 = weighted_diff(cumsum_mat, dy, power = 1)
D_OT2 = weighted_diff(cumsum_mat, dy, power = 2)

# Ensure that the diagonals are zero
diag(D_OT1) = 0
diag(D_OT2) = 0

# Compute kernels
# Simple wealth kernels
K_wealth1 = -0.5 * diag(sqrt(f)) %*% H %*% D_wealth1 %*% t(H) %*% diag(sqrt(f)) # wealth kernel simple 1
K_wealth2 = -0.5 * diag(sqrt(f)) %*% H %*% D_wealth2 %*% t(H) %*% diag(sqrt(f)) # wealth kernel simple 2

#Optimal transport kernels
K_OT1 = -0.5 * diag(sqrt(f)) %*% H %*% D_OT1 %*% t(H) %*% diag(sqrt(f)) # wealth kernel OT1
K_OT2 = -0.5 * diag(sqrt(f)) %*% H %*% D_OT2 %*% t(H) %*% diag(sqrt(f)) # wealth kernel OT2

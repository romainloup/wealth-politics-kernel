## ============================================================
## --- Ressorts politiques ---
## forces politiques -> directions géo normalisées
## + ancrage vers positions initiales
## version 2025/10/10
## ============================================================

## --- Données
# coord_muni : 2126 communes x 2 coordonnées "longitude" et "latitude" 
# dist_pol   : matrice de dissimilarité symétrique, non négative,
#              de diagonale nulle de dimension 2126 x 2126

coord_muni = read.csv("coord_muni.csv", row.names = F)
dist_pol_1 = read.csv("dist_pol_1.csv", row.names = F) # fichier séparé en deux, trop volumineux pour GitHub
dist_pol_2 = read.csv("dist_pol_2.csv", row.names = F)
dist_pol   = rbind(dist_pol_1, dist_pol_2)

## --- Distances euclidiennes en WGS84
euclidean_distance <- function(p1, p2) {
  sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
}

euclidean_dist_matrix <- function(coords) {
  as.matrix(dist(coords))   # en degres (WGS84)
}
## --- Ajustement des positions
adjust_positions_basic_stable <- function(coord_muni, dist_pol,
                                          iterations   = 10,
                                          # force des ressorts politiques
                                          k_spring     = 0.001,
                                          # ancrage vers origine 
                                          # (plus doux par defaut)
                                          k_init       = 0.05,   
                                          # enleve la translation globale
                                          center_each_iter = TRUE,
                                          # en "degrès" (opt.; ex 0.02 ~ qq km)
                                          max_step_per_iter = NA,   
                                          verbose      = TRUE) {
  
  X0 <- as.matrix(coord_muni)              # positions d'origine
  X  <- X0                                 # positions courantes
  n  <- nrow(X)
  eps <- 1e-12
  
  ## 1) Symétriser la matrice politique
  dist_pol <- as.matrix(dist_pol)
  if (!isSymmetric(dist_pol, tol = 1e-12)) {
    dist_pol <- (dist_pol + t(dist_pol)) / 2
  }
  diag(dist_pol) <- 0
  # Normalisation 0-1 sur l'OFF-diagonale
  off <- dist_pol[upper.tri(dist_pol)]
  rng <- range(off, na.rm = TRUE)
  if (diff(rng) > 0) {
    dist_pol_norm <- (dist_pol - rng[1]) / (rng[2] - rng[1])
  } else {
    dist_pol_norm <- dist_pol * 0  # cas dégénéré
  }
  diag(dist_pol_norm) <- 0
  dist_mean <- mean(dist_pol_norm[upper.tri(dist_pol_norm)], na.rm = TRUE)
  
  for (iter in 1:iterations) {
    ## 2) Distances géo courantes (en degrès)
    D_geo <- euclidean_dist_matrix(X)
    diag(D_geo) <- 0
    # éviter divisions par 0
    D_geo[D_geo < eps] <- eps
    
    ## 3) Boucle sur i
    # On accumule les nouvelles positions dans X_new pour
    # éviter d'utiliser des i déjà mis a jour
    X_new <- X
    
    for (i in 1:n) {
      # lignes i
      dist_geo_i    <- D_geo[i, ]                 # n
      dist_target_i <- dist_pol_norm[i, ]         # n
      
      # force politique (centrée par dist_mean)
      force_i <- k_spring * (dist_target_i - dist_mean)  # n
      
      # direction: (x_j - x_i) / d_ij 
      direction_ij <- sweep(X, 2, X[i, ], FUN = "-")     # n x 2
      # division par vecteur (recyclage par lignes)
      direction_ij <- direction_ij / dist_geo_i          
      direction_ij[!is.finite(direction_ij)] <- 0
      
      # ajustement "politique"
      # n x 2 (force_i recycle par colonnes)
      adjustment_mat <- (force_i * direction_ij) 
      adjustment_mat[!is.finite(adjustment_mat)] <- 0
      
      # deplacement total dû aux autres j
      delta_i <- -colSums(adjustment_mat, na.rm = TRUE)  # 1 x 2
      
      # ancrage vers la position d'origine
      d0      <- euclidean_distance(X[i, ], X0[i, ])
      if (d0 > eps) {
        dir0   <- (X0[i, ] - X[i, ]) / d0
      } else {
        dir0   <- c(0, 0)
      }
      adjustment_init <- k_init * dir0 * (d0^2)
      
      # déplacement propose
      move_i <- delta_i + adjustment_init
      
      # cap optionnel (en degrès)
      if (is.finite(max_step_per_iter) && !is.na(max_step_per_iter)) {
        mnorm <- sqrt(sum(move_i^2))
        if (mnorm > max_step_per_iter) {
          move_i <- move_i * (max_step_per_iter / mnorm)
        }
      }
      # appliquer
      X_new[i, ] <- X[i, ] + move_i
    }
    
    # 4) Centrer (supprime la translation globale parasite)
    if (center_each_iter) {
      X_new <- sweep(X_new, 2, colMeans(X_new), "-")
      # recentrer autour du centroide d'origine
      # pour rester proche de la Suisse
      X_new <- X_new + matrix(colMeans(X0), n, 2, byrow = TRUE)
    }
    
    X <- X_new
    
    if (verbose && (iter %% 10 == 0)) {
      cat(sprintf("It %d / %d\n", iter, iterations))
    }
  }
  
  out <- as.data.frame(X)
  names(out) <- names(coord_muni)
  out
}

## ============================================================
## Appel
## ============================================================

# Réglages
# - ancrage doux
# - cap de déplacement
coord_muni_new <- adjust_positions_basic_stable(
  coord_muni  = coord_muni,
  dist_pol   = dist_pol,
  iterations = 50,
  k_spring   = 0.001,
  k_init     = 0.5,
  center_each_iter = TRUE,
  max_step_per_iter = NA,  # 0.02 pour limiter a ~ qq km en lat/long
  verbose    = TRUE
)

## --- Visualisation simple (flèches)
plot(coord_muni$longitude, coord_muni$latitude,
     main = "Carte initiale (bleu) vs ajustee (rouge)",
     col = "blue", pch = 16,
     xlim = range(c(coord_muni$longitude, coord_muni_new$longitude)),
     ylim = range(c(coord_muni$latitude,  coord_muni_new$latitude)))
points(coord_muni_new$longitude, coord_muni_new$latitude, 
       col = "red", pch = 16)
segments(coord_muni$longitude, coord_muni$latitude,
         coord_muni_new$longitude, coord_muni_new$latitude,
         col = adjustcolor("gray30", alpha.f = 0.4))
legend("topright", legend = c("Initial", "Ajuste"),
       col = c("blue","red"), pch = 16)

# --- Visualisation carte
library(ggplot2)
coord_muni_new$language = ch_aggregated_geolevels$language
ggplot() +
  geom_sf(ch,mapping=aes(geometry=geometry),
          fill=alpha("#dfdfdf",0.75),color=alpha("white",0.4)) +
  geom_sf(lakes,mapping=aes(geometry=geometry),
          fill=alpha("#c1c1ce",0.75),color=alpha("white",0.2)) +
  geom_segment(aes(x = coord_muni$longitude, y = coord_muni$latitude, 
                   xend = coord_muni_new$longitude, 
                   yend = coord_muni_new$latitude, 
                   alpha = f), 
               colour = "red") +
  geom_point(data = coord_muni_new, 
             aes(x = longitude, y = latitude, 
                 color = as.character(language), size = f)) +
  scale_color_manual(values = c(
    "1" = "#66C2A5",
    "2" = "#FC8D62",
    "3" = "#8DA0CB",
    "4" = "#E78AC3"),
    name = "Language", 
    labels = c("German", "French", "Italian", "Romansh")
  ) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

# Sauver carte
ggsave("deformation_carte.pdf", width = 10, height = 7)

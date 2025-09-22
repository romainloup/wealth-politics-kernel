# =============================================================================
# 03_produce_results.R
# Compute kernel-level (RV, z) and factor-level (axis-1) associations,
# partial correlations, and weighted regressions. Export tables & figs.
# =============================================================================

# Load functions
source(here::here("R", "RV_funs.R"))

# --- Canonical names for display
dist_names <- c("Political","Road","Size","Language","Wealth1","Wealth2","OT1","OT2")

names(data_list)
# "detailed_language_2024" "distance_mat_2024"      "IFD_tax_classes_2024"
# "IFD_tax_wealth_2024"    "time_mat_2024"          "vote_info_2024"
# "vote_nb_valid_2024"     "vote_nb_yes_2024"       "vote_theme_names"
# "vote_yes_2024"          "f"

# Distance to compute
# Each distance
dist_keys = c("_pol", "_time","_size", "_ling", "_wealth1", "_wealth2", "_OT1", "_OT2")

# --- Run RV function
list_results = RV2(dist_keys, f)
# "Y_list"         "Delta_list"     "E_RV"           "Var_RV"         "Z_RV"           "eigen_val_list" "nb_lambda_pos" 

# --- Language settings
german = which(data_list$detailed_language_2024$language_region == 1)
french = which(data_list$detailed_language_2024$language_region == 2)
italian = which(data_list$detailed_language_2024$language_region == 3)
romansh = which(data_list$detailed_language_2024$language_region == 4)

# add weight
data_list$detailed_language_2024$f = data_list$f

# --- Parameters
# Type to compare (choose the number corresponding to the placement of the kernel dist_keys)
val_1 = 1 # x
val_2 = 3 # y
val_3 = 5 # z
factor_1 = 1
factor_2 = 1
nb_muni = 50

# Select language
select_lang = 1:n # Switzerland (four languages)
select_lang = german
select_lang = french
select_lang = italian
select_lang = romansch


# Run graph
{
  x_axis = list_results$Y_list[[val_1]][select_lang,factor_1]
  x_axis = x_axis/sum(x_axis)
  y_axis = list_results$Y_list[[val_2]][select_lang,factor_2]
  y_axis = y_axis/sum(y_axis)
  
  filtered = as.data.frame(data_list$detailed_language_2024[select_lang,]$f)
  names(filtered) = "f"
  
  filtered$x = x_axis
  filtered$y = y_axis
  
  
  filtered$municipality = data_list$detailed_language_2024[select_lang,]$municipality
  filtered = filtered[filtered$f > sort(filtered$f, decreasing = TRUE)[nb_muni+1], ]
  filtered = filtered[order(filtered$f, decreasing = TRUE),]
  
  lambda_from_RV_1 = list_results$eigen_val_list[[val_1]]$values
  lambda_from_RV_2 = list_results$eigen_val_list[[val_2]]$values
  
  prop_expl_1 = round(100*lambda_from_RV_1 / sum(lambda_from_RV_1), digits = 1 )[factor_1]
  prop_expl_2 = round(100*lambda_from_RV_2 / sum(lambda_from_RV_2), digits = 1 )[factor_2]
  
  x_lab = dist_keys_names[val_1]
  y_lab = dist_keys_names[val_2]
  
  slope_plot <- coef(lm(x_axis ~ y_axis, weights = f[select_lang]*length(select_lang)))
  
  # Graph
  ggplot() +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_point(aes(x = x_axis, y = y_axis, size=f[select_lang], color=mds_I$language[select_lang]),
               alpha = magnif[select_lang]) +
    ggrepel::geom_text_repel(aes(x = filtered$x, y = filtered$y, label = filtered$municipality),
                             box.padding = 0.5,   # Espace autour des étiquettes
                             point.padding = 0.3, # Espace autour des points
                             max.overlaps = Inf ) +
    scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB",  "#E78AC3"),
                       labels = c("German", "French", "Italian","Romansh")) +
    labs(x = paste0(x_lab, ": fact 1: ",prop_expl_1,"%"), y = paste0(y_lab, ": fact 1: ",prop_expl_2,"%")) +
    labs(size = "Reg. weight *f*", color = "Language") +
    scale_size_continuous(range = c(1, 8)) +
    theme_minimal() +
    theme(legend.title = ggtext::element_markdown(lineheight = 1.2)) +
    geom_abline(intercept = slope_plot[1], slope = slope_plot[2], color = "red", linetype="dashed")
}

# ggsave("wasserstein/art_v2/kernel_comparison/mds_Wealth1_OT2.png", width = 9, height = 8)


{
  # --- Regression
  # Weighted regression of the top two MDS factors on all votes
  reg_model_w = lm(x_axis ~ y_axis, weights = f[select_lang]) # faire aussi avec résidu à la place de resultat.jaStimmenInProzent 
  summary(reg_model_w)
  
  # Residuals
  residuals_w = x_axis - predict(reg_model_w)
  
  # Map residuals
  maxData_w = residuals_w[which.max( abs(residuals_w) )]
  
  scale_range_w <- c(-maxData_w, maxData_w)
  palDiv_w <- colorNumeric("PRGn", domain = scale_range_w)
  rev_palDiv_w <- colorNumeric("PRGn", domain = scale_range_w, reverse = TRUE)
  
  leaflet(ch_2024[select_lang,]) %>%
    addProviderTiles(providers$CartoDB.PositronNoLabels) %>%
    setView(lat=46.637785, lng=8.2 , zoom=7) %>%
    # municipality polygons
    addPolygons(
      fillColor = ~palDiv_w(residuals_w),
      fillOpacity = 0.9,
      color = "white",
      weight = 0.1,
      opacity = 1,
      highlight = highlightOptions(
        weight = 2,
        color = "#666",
        fillOpacity = 0.7,
        bringToFront = TRUE,
        sendToBack = TRUE),
      label = paste0(ch_2024$NAME[select_lang], ": ", round(residuals_w, 2)),
      smoothFactor = 0.2,
    ) %>%
    addPolygons(data = lakes,
                weight = 0.7,
                fillColor = "#dddddd",
                fillOpacity = 1,
                color = "white",
                highlight = highlightOptions(
                  weight = 1,
                  color = "#666",
                  fillOpacity = 1,
                  bringToFront = TRUE),
                label = lakes$NAME,
                labelOptions = labelOptions(
                  style = list(
                    "color" = "#666",
                    "font-style" = "italic"
                  )),
                smoothFactor = 0.2,
                group = c("lakes")) %>%
    addLayersControl(
      overlayGroups = c("lakes"),
      options = layersControlOptions(collapsed = TRUE)
    ) %>%
    addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = F)) %>%
    # Ajout de la légende
    addLegend(
      pal = rev_palDiv_w,
      values = scale_range_w,
      position = "bottomright",
      title = "Residuals",
      labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE), suffix = "%"),
      
      opacity = 0.9
    )
}

# Weighted correlation between residuals
weighted.cor <- function(u, v, w) {
  wm <- sum(w)
  mu_u <- sum(w * u) / wm
  mu_v <- sum(w * v) / wm
  cov_uv <- sum(w * (u - mu_u) * (v - mu_v)) / wm
  sd_u <- sqrt(sum(w * (u - mu_u)^2) / wm)
  sd_v <- sqrt(sum(w * (v - mu_v)^2) / wm)
  cov_uv / (sd_u * sd_v)
}


# Weighted regression of the top two MDS factors on all votes
z = list_results$Y_list[[val_3]][select_lang,factor_1]
x = list_results$Y_list[[val_1]][select_lang,factor_1]
reg_model_x = lm(x ~ z, weights = f[select_lang])
# Residuals
residuals_x = x - predict(reg_model_x)

# Reg 2
y = list_results$Y_list[[val_2]][select_lang,factor_1]
reg_model_y = lm(y ~ z, weights = f[select_lang])
# Residuals
residuals_y = y - predict(reg_model_y)

partial_cor <- weighted.cor(residuals_x, residuals_y, f[select_lang])
partial_cor

weighted.cor(x, y, f[select_lang])

dim_mat = length(dist_keys_names)
corr_mat_w = matrix(data = NA, dim_mat,dim_mat)
for (i in 1:dim_mat) {
  x = list_results$Y_list[[i]][select_lang,factor_1]
  for (j in 1:dim_mat) {
    y = list_results$Y_list[[j]][select_lang,factor_1]
    corr_mat_w[i,j] = weighted.cor(x, y, f[select_lang])
  }
}
rownames(corr_mat_w) = dist_keys_names
colnames(corr_mat_w) = dist_keys_names

plot(x,y)

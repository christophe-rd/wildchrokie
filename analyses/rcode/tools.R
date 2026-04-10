# CRD 2 March 2026
library(wesanderson)
# saving my r function because I'm tired of writing it in every code
extract_params <- function(df_fit, pattern, col_prefix, id_col, rename_pattern = NULL) {
  cols <- colnames(df_fit)[grepl(pattern, colnames(df_fit))]
  cols <- cols[!grepl("prior", cols)]  # exclude prior predictive columns
  df <- df_fit[, cols, drop = FALSE]
  
  if (!is.null(rename_pattern)) {
    colnames(df) <- sub(rename_pattern, "\\1", colnames(df))
  }
  
  result <- do.call(rbind, lapply(colnames(df), function(col) {
    x <- df[[col]]
    data.frame(
      id    = col,
      mean  = round(mean(x), 3),
      per5  = round(quantile(x, 0.05), 3),
      per25 = round(quantile(x, 0.25), 3),
      per75 = round(quantile(x, 0.75), 3),
      per95 = round(quantile(x, 0.95), 3),
      row.names = NULL
    )
  }))
  
  colnames(result)[1] <- id_col
  colnames(result)[-1] <- paste0(col_prefix, c("", "_per5", "_per25", "_per75", "_per95"))
  result
}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Figures and stuff ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### Wildchrokie #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# col
wccolslatbi <- c(
  "Alnus incana" = wes_palette("AsteroidCity1")[1],
  "Betula alleghaniensis" = wes_palette("AsteroidCity1")[2],
  "Betula papyrifera" = wes_palette("AsteroidCity1")[3],
  "Betula populifolia" = wes_palette("AsteroidCity1")[4]
)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### CoringTreespotters #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# specify colors
renoir <- c("#17154f", "#2f357c", "#6c5d9e", "#9d9cd5", "#b0799a", "#e48171", 
            "#bf3729", "#e69b00", "#f5bb50", "#ada43b", "#355828")

tscolscommon <- c(
  "Red maple"           = renoir[1],
  "Sugar maple"         = renoir[2],
  "Yellow buckeye"      = renoir[3],
  "Yellow birch"        = renoir[4],
  "River birch"         = renoir[5],
  "Pignut hickory"      = renoir[6],
  "Shagbark hickory"    = renoir[7],
  "Eastern cottonwood"  = renoir[8],
  "White oak"           = renoir[9],
  "Northern red oak"    = renoir[10],
  "American basswood"   = renoir[11]
)

tscolslatbi <- c(
  "Acer rubrum"           = renoir[1],
  "Acer saccharum"        = renoir[2],
  "Aesculus flava"        = renoir[3],
  "Betula alleghaniensis" = renoir[4],
  "Betula nigra"          = renoir[5],
  "Carya glabra"          = renoir[6],
  "Carya ovata"           = renoir[7],
  "Populus deltoides"     = renoir[8],
  "Quercus alba"          = renoir[9],
  "Quercus rubra"         = renoir[10],
  "Tilia americana"       = renoir[11]
)
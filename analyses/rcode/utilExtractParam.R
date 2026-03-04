# CRD 2 March 2026
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
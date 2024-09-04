## functions to format model output for plotting and analysis

## dependencies: tidyverse; plyr

format_output <- function(m,res){

  res_rn <- rename_cols_rows(m,res)
  res_df <- make_df_from_array(m,res_rn)
  res_df_long <- df_long_form(res_df)

  return(res_df_long)

}

rename_cols_rows <- function(m,res){

  rownames(res) <- names(unlist(m$info()$index))
  colnames(res) <- paste0("P",1:ncol(res))

  return(res)
}

make_df_from_array <- function(m,res){

  # collapse from array into df format
  # add a column for variable as this can no longer be a rowname
  res_df <- plyr::adply(res,3,.id = "t") %>%
    mutate(t = as.numeric(t),
           "var" = rep(names(unlist(m$info()$index)),max(t))) %>%
    filter(var!="time")

  return(res_df)
}

df_long_form <- function(res_df){

  # collapse particle realisations
  res_df_long <- pivot_longer(res_df,
                            names_to = "particle",values_to="value",
                            cols=starts_with("P"))

  return(res_df_long)

}



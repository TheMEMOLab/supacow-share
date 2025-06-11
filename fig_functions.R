# General use functions and objects for various scripts

# RCT colors + shapes for general use
RCT_cols <- c("#bfbfbf", "#49ad8b")
RCT_shapes <- c(24, 21)

# General purpose variable name vector
# (corresponding to main metadata file)
sc_vars <- c("animal", "hoko_id", "uk_id", "bcms", "treatment",
             "breed", "sire", "dam", "dob", "updated_pen",
             "performance_test_start_date", "performance_test_end_date",
             "age_start_performance_test", "performance_test_lwg", "performance_test_dmi",
             "fcr", "liveweight", "ch4_g_day", "ch4_g_kg_dmi", "chosen_24",
             "chosen_6", "RCT", "ch4_dmi_binary")
sc_vars <- setNames(sub("Performance test end date", "Performance test group",
                        sub("dmi", "DMI", ignore.case = TRUE,
                            sub("fcr", "FCR", ignore.case = TRUE,
                                sub("Ch4", "CH4",
                                    sub("Updated pen", "Pen",
                                        capitalize(
                                          gsub("_", " ",
                                               sub("g_", "g\\/",
                                                   sub("lwg", "ADG",
                                                       sc_vars))))))))),
                    sc_vars)

# Simple Wilcoxon rank sum tests for comparing rumen community types
basic_rct_wilcox <- function(df, vars){
  
  # Calculate test and get statistics
  wilcoxres <- do.call("rbind", lapply(vars, function(x){
    wc <- wilcox.test(as.formula(paste(x, "~RCT")), df)
    data.frame(feature = x,
               W_statistic = wc$statistic,
               p = wc$p.value)
    }))
  
  # Multiple comparison correction & additional metrics
  res <- data.frame(
    wilcoxres,
    padj = p.adjust(wilcoxres$p, "fdr"),
    gr1_median = colMedians(
      as.matrix(
        subset(df, RCT == "A")[,wilcoxres$feature])),
    gr2_median = colMedians(
      as.matrix(
        subset(df, RCT == "B")[,wilcoxres$feature]))
  )
  
  # RCT direction
  res$up_in <- factor(
    ifelse(res$gr1_median > res$gr2_median,
           "Enriched in\nRCT-A",
           "Enriched in\nRCT-B"))
  
  return(res)
}

# RCT summary for medians + IQR
library("reshape2")
rct_med_iqr <- function(df){
  do.call(
  "rbind", lapply(c("A", "B"), function(x){
    df_gr <- subset(df, RCT == x)
    df_m <- reshape2::melt(df_gr)
    res <- aggregate(value ~ variable, df_m, summary)
    res <- data.frame(res$variable,
                      unlist(res$value))
    res$RCT <- x
    return(res)
  }))
}

# Simple ordination function with robust Aitchison distance
roba_ord <- function(phy){
  # Calculate distance matrix
  ra <- t(prop.table(otu_table(phy), 2))
  dist <- vegdist(ra, method = "robust.aitchison")
  
  # Calculate ordination
  list(dist = dist,
       pcoa = pcoa(dist))}

# Generic function to calculate PCA from selected group_index in imputed proteomics file
# (different "group indices" corresponding to different sample types and imputation
#  with breeds either separated or together)
prot_to_pca <- function(df, gr_i, ntop=1000){
  
  df_sub <- subset(df, group_index == gr_i & database != "contaminant")
  df_final <- reshape2::dcast(df_sub, sample ~ protein, value.var = "intensity")
  
  # sample name matching to generic metadata
  rownames(df_final) <- sub("^[A-Z]", "D", sub("[A-Z]$", "", df_final$sample))
  df_final$sample <- NULL
  
  df_final_vars <- colVars(as.matrix(df_final))
  
  pca <- prcomp(
    df_final[,order(df_final_vars, decreasing = TRUE)[1:ntop]],
    scale = FALSE)
}

# Generic box and dot plot
box_dot_plot <- function(df, varx, vary, cols){
  ggplot(df, aes(y = !!sym(vary), x = !!sym(varx),
                 fill = !!sym(varx))) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center") +
    theme_classic(base_size = 8) +
    xlab(sc_vars[varx]) +
    scale_fill_manual(values = cols) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 12))
}
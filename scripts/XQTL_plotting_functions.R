XQTL_Manhattan_5panel <- function(df, cM = FALSE) {
  # Define UC colors
  uc_colors <- c("#003262", "#FDB515")
  
  # Order chromosomes
  chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
  df$chr <- factor(df$chr, levels = chr_order)
  
  x_var <- if (cM) "cM" else "Mb"
  df$Mb <- df$pos / 1e6
  
  # Create label data and calculate y-axis limits
  label_data <- df %>% 
    group_by(chr) %>% 
    summarise(
      x = max(!!sym(x_var)), 
      y = max(Wald_log10p),
      y_max = max(10, ceiling(max(Wald_log10p)))  # Ensure y_max is at least 10
    )
  
  p <- ggplot(df, aes(x = !!sym(x_var), y = Wald_log10p)) +
    geom_point(color = uc_colors[1], size = 0.25) +
    facet_wrap(~ chr, ncol = 1, scales = "free") +  # Changed to free scales for both x and y
    geom_text(
      data = label_data,
      aes(x = Inf, y = Inf, label = chr),
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +
    labs(x = x_var, y = "-log10(p-value)") +
    theme_bw() +
    theme(
      panel.spacing = unit(0.1, "lines"),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  # Set y-axis limits for each facet
  p <- p + scale_y_continuous(limits = function(y) c(0, max(10, ceiling(max(y)))))
  
  return(p)
}

XQTL_Manhattan <- function(df, cM = FALSE, color_scheme = "KU") {
  # Define color schemes
  color_schemes <- list(
    KU = c("#0051BA", "#E8000D"),
    UCI = c("#003262", "#FDB515"),
    Stanford = c("#8C1515", "#4D4F53"),
    Harvard = c("#A51C30", "#C4D600"),
    MIT = c("#A31F34", "#8A8B8C"),
    Berkeley = c("#003262", "#FDB515"),
    McMaster = c("#7A003C", "#FDBF57"),
    McGill = c("#ED1B2F", "#FFD794"),
    Oxford = c("#002147", "#1C4E91"),
    Cambridge = c("#A3C1AD", "#D6083B"),
    NineInchNails = c("#000000", "#FF0000")
  )
  
  # Select color scheme
  if (tolower(color_scheme) %in% tolower(names(color_schemes))) {
    chosen_colors <- color_schemes[[which(tolower(names(color_schemes)) == tolower(color_scheme))]]
  } else {
    stop("Invalid color scheme. Choose from: ", paste(names(color_schemes), collapse = ", "))
  }

  # Rest of the function remains the same
  # ORDER the chromosome arms
  chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
  df$chr <- factor(df$chr, levels = chr_order)

  # GENETIC DISTANCE plot
  if (cM) {
    # CALCULATE cumulative cM distance
    df <- df %>%
      group_by(chr) %>%
      mutate(chr_max_cM = max(cM)) %>%
      ungroup() %>%
      # THE "* 1.10" ADDS whitespace between chromosomes
      mutate(cumulative_cM = case_when(
        chr == "chrX" ~ cM,
        chr %in% c("chr2L", "chr2R") ~ cM +
          first(chr_max_cM[chr == "chrX"]) * 1.10,
        chr %in% c("chr3L", "chr3R") ~ cM +
          first(chr_max_cM[chr == "chrX"]) * 1.10 + 
          max(chr_max_cM[chr %in% c("chr2L", "chr2R")]) * 1.10))
    
    x_var <- "cumulative_cM"
    x_lab <- "Genetic Distance (cM)"
  }

  # PHYSICAL DISTANCE plot
  else {
    # Calculate Mb and cumulative_Mb
    df <- df %>%
      mutate(Mb = pos / 1e6) %>%
      group_by(chr) %>%
      mutate(chr_max_Mb = max(Mb)) %>%
      ungroup() %>%
      # THE "* 1.05" ADDS whitespace between chromosome arms
      mutate(cumulative_Mb = case_when(
        chr == "chrX" ~ Mb,
        chr == "chr2L" ~ Mb +
          first(chr_max_Mb[chr == "chrX"]) * 1.05,
        chr == "chr2R" ~ Mb + (first(chr_max_Mb[chr == "chrX"]) +
          first(chr_max_Mb[chr == "chr2L"])) * 1.05,
        chr == "chr3L" ~ Mb + (first(chr_max_Mb[chr == "chrX"]) +
          first(chr_max_Mb[chr == "chr2L"]) +
          first(chr_max_Mb[chr == "chr2R"])) * 1.05,
        chr == "chr3R" ~ Mb + (first(chr_max_Mb[chr == "chrX"]) +
          first(chr_max_Mb[chr == "chr2L"]) +
          first(chr_max_Mb[chr == "chr2R"]) +
          first(chr_max_Mb[chr == "chr3L"])) * 1.05))
    
    x_var <- "cumulative_Mb"
    x_lab <- "Physical Distance (Mb)"
  }

  # GET chromosome midpoints for chromosome ID x-axis labels
  chr_midpoints <- df %>%
    group_by(chr) %>%
    summarize(mid = min(!!sym(x_var)) + (max(!!sym(x_var)) - min(!!sym(x_var)))/2)

  # REMOVE all rows containing an NA value
  # (prevents "warning" on plotting that rows will be removed)
  df <- na.omit(df)

  # PLOT
  p <- ggplot(df, aes(x = !!sym(x_var), y = Wald_log10p, color = chr)) +
    geom_point(size = 0.25) +
    scale_color_manual(values = rep(chosen_colors, 3)) +
    scale_x_continuous(breaks = chr_midpoints$mid,
                       labels = gsub("chr","",chr_midpoints$chr)) +
    labs(x = x_lab, y = expression(paste(-log[10],italic(P)))) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(linewidth=0.4),
          panel.grid.minor.y = element_line(linewidth=0.4),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1,
                                     size = 11, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"))

  return(p)
}

	
get_palette <- function(founders, reference_strain = NULL) {
  base_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#000000", "#990099")
  
  if (!is.null(reference_strain) && reference_strain %in% founders) {
    # Add grey for the reference strain
    palette <- c(base_palette[1:(length(founders)-1)], "#999999")
    names(palette) <- c(setdiff(founders, reference_strain), reference_strain)
  } else {
    palette <- base_palette[1:length(founders)]
    names(palette) <- founders
  }
  
  return(palette)
}

XQTL_change_average <- function(df, chr, start, stop, reference_strain = NULL) {
  # Subset the dataframe
  subset_df <- df %>% filter(chr == !!chr, pos > start, pos < stop)
  
  # Pivot the dataframe to wide format and calculate Dfreq
  wide_df <- subset_df %>% 
    pivot_wider(names_from = TRT, values_from = freq, names_prefix = "freq_") %>% 
    mutate(Dfreq = freq_Z - freq_C)
  
  # Calculate average Dfreq over REP
  avg_df <- wide_df %>% 
    group_by(chr, pos, founder) %>% 
    summarize(Dfreq = mean(Dfreq, na.rm = TRUE), .groups = "drop")
  
  # Get the color palette
  color_palette <- get_palette(unique(avg_df$founder), reference_strain)
  
  # Create the plot
  p <- ggplot(avg_df, aes(x = pos, y = Dfreq, color = founder)) + 
    geom_line() +  # Keep original line thickness for the plot
    labs(title = paste("Average Frequency Change by Position (", chr, ")"),
         x = "Position", 
         y = "Δ Frequency (Z - C)") + 
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_color_manual(values = color_palette) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.box = "horizontal",
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Reduced top margin
      legend.spacing.x = unit(0.2, 'cm'),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Added this line
      plot.margin = margin(t = 10, r = 10, b = 20, l = 10, unit = "pt"),  # Adjusted bottom margin
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))  # Added space above x-axis label
      ) +
      guides(color = guide_legend(
      override.aes = list(linewidth = 3),  # Use linewidth instead of size
      nrow = 1
    ))
  
  return(p)
}

XQTL_change_byRep <- function(df, chr, start, stop) {
  # Subset the dataframe
  subset_df <- df %>%
    filter(chr == !!chr, pos > start, pos < stop)
  
  # Pivot the dataframe to wide format and calculate Dfreq
  wide_df <- subset_df %>%
    pivot_wider(names_from = TRT, values_from = freq, names_prefix = "freq_") %>%
    mutate(Dfreq = freq_Z - freq_C)
  
  # Calculate average freq_C by founder and select top 4
  top_founders <- wide_df %>%
    group_by(founder) %>%
    summarize(mean_freq_C = mean(freq_C, na.rm = TRUE)) %>%
    top_n(4, mean_freq_C) %>%
    pull(founder)
  
  # Filter the dataframe for top 4 founders
  plot_df <- wide_df %>%
    filter(founder %in% top_founders)
  
  # Create the 4-panel plot
  p <- ggplot(plot_df, aes(x = pos, y = Dfreq, color = REP, group = REP)) +
    geom_line() +
    facet_wrap(~ founder, ncol = 2) +
    labs(title = paste("Frequency Change by Position and Rep (", chr, ")"),
         x = "Position",
         y = "Δ Frequency (Z - C)",
         color = "Rep") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

XQTL_beforeAfter_selectReps <- function(df, chr, start, stop, reps) {
  # Subset the dataframe
  subset_df <- df %>%
    filter(chr == !!chr, pos > start, pos < stop, REP %in% reps)
  
  # Calculate average freq_C by founder and select top 4
  top_founders <- subset_df %>%
    filter(TRT == "C") %>%
    group_by(founder) %>%
    summarize(mean_freq_C = mean(freq, na.rm = TRUE)) %>%
    top_n(4, mean_freq_C) %>%
    pull(founder)
  
  # Filter the dataframe for top 4 founders
  plot_df <- subset_df %>%
    filter(founder %in% top_founders)
  
  # Create a distinct color palette for reps
  rep_colors <- brewer.pal(length(reps), "Set1")
  names(rep_colors) <- reps
  
  # Create the 4-panel plot
  p <- ggplot(plot_df, aes(x = pos, y = freq, color = as.factor(REP))) +
    geom_line(data = . %>% filter(TRT == "C"), linetype = "dashed") +
    geom_line(data = . %>% filter(TRT == "Z"), linetype = "solid") +
    facet_wrap(~ founder, ncol = 2) +
    scale_color_manual(values = rep_colors) +
    labs(title = paste("Frequency Before and After Treatment (", chr, ")"),
         x = "Position",
         y = "Frequency",
         color = "Rep") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}


XQTL_region <- function(df, chr, start, stop, trait) {
  # UC Berkeley blue
  uc_blue <- "#003262"
  
  # Subset the dataframe
  subset_df <- df %>%
    filter(chr == !!chr, pos >= start, pos <= stop)
  
  # Create the plot
  p <- ggplot(subset_df, aes(x = pos, y = !!sym(trait))) +
    geom_line(color = uc_blue) +
    labs(title = paste("Region plot for", trait, "on", chr),
         x = "Position",
         y = trait) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  return(p)
}

XQTL_combined_plot <- function(df1, df3, chr, start, stop, reference_strain = NULL) {
  
  # Create sqrt_avg_var if it doesn't exist
  if(!"sqrt_avg_var" %in% names(df1)) {
    df1 <- df1 %>% mutate(sqrt_avg_var = sqrt(avg.var))
  }

  # Create individual plots
  p1 <- XQTL_region(df1, chr, start, stop, "Wald_log10p")
  p2 <- XQTL_region(df1, chr, start, stop, "Cutl_H2")
  p3 <- XQTL_region(df1, chr, start, stop, "sqrt_avg_var")
  p4 <- XQTL_change_average(df3, chr, start, stop, reference_strain = reference_strain)  # Pass reference_strain here
  
  # Function to adjust plots
  adjust_plot <- function(p, remove_x = TRUE) {
    p + 
      theme(
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
  }
  
  # Adjust individual plots
  p1 <- adjust_plot(p1)
  p2 <- adjust_plot(p2)
  p3 <- adjust_plot(p3)
  p4 <- p4 + 
    theme(
      plot.title = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
  
  # Combine plots
  combined_plot <- p1 / p2 / p3 / p4 +
    plot_layout(heights = c(1, 1, 1, 1.5))
  
  # Final adjustments
  final_plot <- combined_plot +
    plot_annotation(
      theme = theme(
        plot.margin = margin(5, 5, 5, 5)
      )
    )
  
  return(final_plot)
}


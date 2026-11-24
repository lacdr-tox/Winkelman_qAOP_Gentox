#Documentation on the functions used for response-response analysis as described in 'Advancing AOPs for genotoxicity through quantifying KERs linking DNA alkylation to increase in mutations and chromosome aberrations' by Winkelman et al., 2026.

match_study_data = function(df1, df2) {
  
  # ---- Merge individual studies ----
  df1 = merge_study_data(df1)
  df2 = merge_study_data(df2)
  
  # ---- Determine response columns ----
  response1 = if("FC" %in% names(df1)) "FC" else "outcome_treatment"
  response2 = if("FC" %in% names(df2)) "FC" else "outcome_treatment"
  
  
  # ---- Choose reference dataset ----
  if(length(unique(df1$Concentration)) >=
     length(unique(df2$Concentration))) {
    
    ref_concs = sort(unique(df1$Concentration))
    target_df = df2
    target_name = "df2"
    
  } else {
    
    ref_concs = sort(unique(df2$Concentration))
    target_df = df1
    target_name = "df1"
  }
  
  
  # ---- Map target concentrations onto reference ----
  target_df$mapped_conc = NA_real_
  
  # All possible target-reference combinations
  possible_matches = expand.grid(
    target = seq_len(nrow(target_df)),
    reference = seq_along(ref_concs)
  )
  
  possible_matches$target_conc =
    target_df$Concentration[possible_matches$target]
  
  possible_matches$ref_conc =
    ref_concs[possible_matches$reference]
  
  # Relative distance
  possible_matches$distance =
    ifelse(
      possible_matches$ref_conc == 0,
      abs(possible_matches$target_conc),
      abs(possible_matches$target_conc -
            possible_matches$ref_conc) /
        abs(possible_matches$ref_conc)
    )
  
  # Keep only valid matches within 10% tolerance
  possible_matches =
    possible_matches[
      !is.na(possible_matches$distance) &
        possible_matches$distance <= 0.10,
    ]
  
  # Closest matches first
  possible_matches =
    possible_matches[
      order(possible_matches$distance),
    ]
  
  # Track already-used target and reference concentrations
  matched_targets = integer(0)
  matched_references = integer(0)
  
  # Assign matches
  if(nrow(possible_matches) > 0) {
    
    for(i in seq_len(nrow(possible_matches))) {
      
      target = possible_matches$target[i]
      reference = possible_matches$reference[i]
      
      # Don't allow either concentration to be used twice
      if(target %in% matched_targets ||
         reference %in% matched_references) {
        next
      }
      
      target_df$mapped_conc[target] =
        possible_matches$ref_conc[i]
      
      matched_targets =
        c(matched_targets, target)
      
      matched_references =
        c(matched_references, reference)
    }
  }
  
  
  # ---- Put mapped dataset back ----
  target_df$Concentration = target_df$mapped_conc
  target_df$mapped_conc = NULL
  
  if(target_name == "df1") {
    df1 = target_df
  } else {
    df2 = target_df
  }
  
  
  # ---- Keep only common concentrations ----
  common_concs = intersect(
    df1$Concentration,
    df2$Concentration
  )
  
  df1 = df1[df1$Concentration %in% common_concs, ]
  df2 = df2[df2$Concentration %in% common_concs, ]
  
  
  # ---- Keep only required columns ----
  df1 = df1[, c(
    "Concentration",
    "total_replicates",
    response1,
    "SD",
    "SEM"
  )]
  
  df2 = df2[, c(
    "Concentration",
    "total_replicates",
    response2,
    "SD",
    "SEM"
  )]
  
  
  # ---- Rename columns ----
  names(df1) = c(
    "Concentration",
    "total_replicates_1",
    "Response_1",
    "SD_1",
    "SEM_1"
  )
  
  names(df2) = c(
    "Concentration",
    "total_replicates_2",
    "Response_2",
    "SD_2",
    "SEM_2"
  )
  
  
  # ---- Merge ----
  matched_df = merge(
    df1,
    df2,
    by = "Concentration",
    sort = TRUE
  )
  
  return(matched_df)
}

orthogonal_residuals = function(par,
                                x,
                                response,
                                SD_x,
                                SD_y,
                                model_fun,
                                positive_x = FALSE,
                                search_width = 2,
                                eps = 1e-8) {
  
  # Prevent division by zero
  SD_x = pmax(SD_x, eps)
  SD_y = pmax(SD_y, eps)
  
  res = numeric(length(response))
  x_true = numeric(length(response))
  
  for(i in seq_along(response)) {
    
    x_obs = x[i]
    y_obs = response[i]
    
    sx = SD_x[i]
    sy = SD_y[i]
    
    # Search interval for latent x*
    lower = x_obs - search_width * sx
    upper = x_obs + search_width * sx
    
    if(positive_x) {
      lower = max(lower, eps)
      upper = max(upper, lower + eps)
    }
    
    # Objective function
    obj = function(x_true) {
      
      y_hat = tryCatch(
        model_fun(par, x_true),
        error = function(e) NA_real_
      )
      
      if(!is.finite(y_hat))
        return(1e12)
      
      ((y_obs - y_hat)^2 / sy^2) +
        ((x_true - x_obs)^2 / sx^2)
    }
    
    opt = try(
      optimize(obj, lower = lower, upper = upper),
      silent = TRUE
    )
    
    if(inherits(opt, "try-error") || !is.finite(opt$objective)) {
      res[i] = 1e6
    } else {
      res[i] = sqrt(opt$objective)
    }
    
    x_true[i] = opt$minimum
    
  }
  
  res[!is.finite(res)] = 1e6
  
  return(res)
}

fit_models_ODR = function(x,
                          y,
                          SD_x,
                          SD_y,
                          models,
                          start_list,
                          lower_bounds,
                          upper_bounds,
                          fixed1 = FALSE,
                          biasCorrBG = NULL,
                          positive_x = FALSE) {
  
  fits = lapply(names(models), function(name) {
    
    m = models[[name]]
    
    sp    = unlist(start_list[m$pars])
    lower = unlist(lower_bounds[m$pars])
    upper = unlist(upper_bounds[m$pars])
    
    fit = try(
      
      nls.lm(
        par   = sp,
        lower = lower,
        upper = upper,
        fn    = orthogonal_residuals,
        x = x,
        response = y,
        SD_x = SD_x,
        SD_y = SD_y,
        model_fun = m$fun,
        positive_x = positive_x
      ),
      
      silent = TRUE
    )
    
    if(inherits(fit, "try-error") || is.null(fit$par)) {
      fit = NULL
    }
    
    if(!is.null(fit) && fixed1){
      
      fit$par = c(
        fit$par,
        c = biasCorrBG
      )
      
      fit$par = fit$par[order(names(fit$par))]
    }
    
    list(
      name = name,
      fit = fit
    )
    
  })
  
  names(fits) = names(models)
  
  fits
}

plot_RR_fit = function(df,
                       fit,
                       model_fun,
                       xcol = "Response_1",
                       ycol = "Response_2",
                       sxcol = "SD_1",
                       sycol = "SD_2",
                       measured_col = "Measured",
                       main = "",
                       xlab = "",
                       ylab = "",
                       sqrtSDplot = TRUE,
                       log_x = TRUE,
                       zoom = FALSE){
  
  measured = df[[measured_col]]
  
  x  = df[[xcol]]
  y  = df[[ycol]]
  sx = df[[sxcol]]
  sy = df[[sycol]]
  
  if (sqrtSDplot) {
    sx = sqrt(sx)
    sy = sqrt(sy)
  }
  
  xmin = x - sx
  xmax = x + sx
  ymin = y - sy
  ymax = y + sy
  
  ## Prediction
  x_pred = seq(min(x), max(x), length.out = 300)
  y_pred = model_fun(fit$par, x_pred)
  
  df_pred = data.frame(
    x = x_pred,
    y = y_pred
  )
  
  ## Cap lengths proportional to scale
  x_cap = 0.015 * diff(range(y))
  y_cap = 0.015 * diff(range(x))
  
  df_plot = data.frame(
    x,
    y,
    xmin,
    xmax,
    ymin,
    ymax,
    measured,
    x_cap,
    y_cap
  )
  
  df_plot = df_plot[order(df_plot$measured), ]
  
  p =
    ggplot() +
    
    ## x = y reference line
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "lightgray",
      linewidth = ifelse(zoom, 1, 0.5)
    )
  
  ## Add threshold area
  if (zoom) {
    p = p +
      annotate(
        "rect",
        xmin = -Inf, xmax = 1.5,
        ymin = -Inf, ymax = 1.5,
        fill = "red",
        alpha = 0.1
      )
  }
  
  p = p +
    
    ## Prediction
    geom_line(
      data = df_pred,
      aes(x, y),
      colour = "blue",
      linewidth = 1
    ) +
    
    ## Horizontal error bars
    geom_errorbar(
      data = df_plot,
      aes(
        y = y,
        xmin = xmin,
        xmax = xmax,
        colour = measured
      ),
      orientation = "y",
      width = 0
    ) +
    
    ## Vertical error bars
    geom_errorbar(
      data = df_plot,
      aes(
        x = x,
        ymin = ymin,
        ymax = ymax,
        colour = measured
      ),
      width = 0
    ) +
    
    ## Points
    geom_point(
      data = df_plot,
      aes(x, y, colour = measured),
      size = 2.5
    ) +
    
    scale_colour_manual(
      values = c(
        "TRUE" = "red",
        "FALSE" = "black"
      ),
      guide = "none"
    ) +
    
    labs(
      x = xlab,
      y = ylab,
      title = main
    ) +
    
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  if (log_x)
    p = p + scale_x_log10()
  
  p
}

fit_response_response = function(
    df1,
    df2,
    wpinput1,
    wpinput2,
    main = "",
    xlab = "",
    ylab = "",
    sqrtSDplot = T,
    log_x = FALSE,
    nboot = 200,
    fixed1 = FALSE,
    conc_col  ='Concentration EMS',
    savename = NULL,
    savepath = '../Figures/Response-Response/EMS/',
    log_conc_space = F,
    zoom = NULL){
  
  if (is.null(savename)) {
    savename = paste(
      deparse(substitute(df1)),
      deparse(substitute(df2)),
      sep = "_"
    )
  }
  
  matched = match_study_data(df1, df2)
  
  n_existing = nrow(matched)
  n_needed = ifelse(n_existing < 10, 10 - n_existing, 0)
  
  common_concrange = sort(union(df1[[conc_col]], df2[[conc_col]]))
  
  fitdose_joint = exp(seq(0, log(max(common_concrange) + 1), length.out = 500)) - 1
  n_fitdose = length(fitdose_joint)
  boot_wpBT1 = boot_wpBT2 = numeric(n_fitdose)
  
  nboot = 200
  
  boot_preds1 = boot_preds2 = matrix(NA, nrow = nboot, ncol = length(fitdose_joint))
  
  
  response_col1 = if("FC" %in% names(df1)) "FC" else "Outcome for the treatment"
  response_col2 = if("FC" %in% names(df2)) "FC" else "Outcome for the treatment"
  
  intervalRes1 = c(0, 50 * max(df1[[response_col1]]))
  intervalRes2 = c(0, 50 * max(df2[[response_col2]]))
  
  
  for(i in 1:nboot){
    boot_wp1 = try(
      weighted_prediction(fitdose_joint,
                          wpinput1$fits[[i]],
                          wpinput1$models,
                          wpinput1$best_models[[i]],
                          wpinput1$weights[[i]]),
      silent = TRUE
    )
    
    boot_wp2 = try(
      weighted_prediction(fitdose_joint,
                          wpinput2$fits[[i]],
                          wpinput2$models,
                          wpinput2$best_models[[i]],
                          wpinput2$weights[[i]]),
      silent = TRUE
    )
    
    proxySDFullrange1 = as.numeric(wpinput1$boot_global[[i]]$icept +
                                     wpinput1$boot_global[[i]]$boot_slope *
                                     fitdose_joint)
    
    proxySDFullrange2 = as.numeric(wpinput2$boot_global[[i]]$icept +
                                     wpinput2$boot_global[[i]]$boot_slope *
                                     fitdose_joint)
    
    for(k in seq_len(n_fitdose)) {
      
      boot_wpBT1[k] = backtransform_mu(
        boot_wp1[k],
        SD = proxySDFullrange1[k],
        interval = intervalRes1
      )
      
      boot_wpBT2[k] = backtransform_mu(
        boot_wp2[k],
        SD = proxySDFullrange2[k],
        interval = intervalRes2
      )
    }
    
    boot_preds1[i, ] = boot_wpBT1
    boot_preds2[i, ] = boot_wpBT2
    
  }
  
  valid_rows1 = !is.na(wpinput1$BMD) & apply(!is.na(boot_preds1), 1, any)
  valid_rows2 = !is.na(wpinput2$BMD) & apply(!is.na(boot_preds2), 1, any)
  
  boot_preds1 = boot_preds1[valid_rows1,]
  boot_preds2 = boot_preds2[valid_rows2,]
  
  nb_valid1 = nrow(boot_preds1)
  nb_valid2 = nrow(boot_preds2)
  
  # Build long data.frame for plotting
  df_boot_long1 = melt(boot_preds1)
  df_boot_long2 = melt(boot_preds2)
  
  colnames(df_boot_long1) = c("sample", "fitdose_idx", "response")
  colnames(df_boot_long2) = c("sample", "fitdose_idx", "response")
  
  df_boot_long1$fitdose_idx = as.integer(df_boot_long1$fitdose_idx)
  df_boot_long2$fitdose_idx = as.integer(df_boot_long2$fitdose_idx)
  
  df_boot_long1$Concentration = fitdose_joint[df_boot_long1$fitdose_idx]
  df_boot_long2$Concentration = fitdose_joint[df_boot_long2$fitdose_idx]
  
  df_boot_long1$sample = factor(df_boot_long1$sample)
  df_boot_long2$sample = factor(df_boot_long2$sample)
  
  #Remove rows for which calculation failed (perhaps consider taking out the same rows for both dataframes at a later stage?)
  df_boot_long1 = df_boot_long1[complete.cases(df_boot_long1),]
  df_boot_long2 = df_boot_long2[complete.cases(df_boot_long2),]
  
  ##########Section 2##########
  
  df1 = matched
  df2 = df_boot_long1
  df3 = df_boot_long2
  
  combined = df1 %>%
    select(Concentration, Response_1, Response_2, SD_1, SD_2)
  
  summary_df2 = df2 %>%
    group_by(Concentration) %>%
    summarise(
      Response_1 = mean(response),
      SD_1 = sd(response),
      .groups = "drop"
    )
  
  summary_df3 = df3 %>%
    group_by(Concentration) %>%
    summarise(
      Response_2 = mean(response),
      SD_2 = sd(response),
      .groups = "drop"
    )
  
  candidate_df = summary_df2 %>%
    inner_join(summary_df3, by = "Concentration")
  
  selected_concs = combined$Concentration
  
  candidate_df = candidate_df %>%
    filter(!Concentration %in% selected_concs)
  
  for(i in seq_len(n_needed)) {
    
    if(nrow(candidate_df) == 0)
      break
    
    # Compute distances on the requested scale
    if(log_conc_space) {
      selected = log(selected_concs)
      candidates = log(candidate_df$Concentration)
    } else {
      selected = selected_concs
      candidates = candidate_df$Concentration
    }
    
    # Minimum distance to already selected concentrations
    candidate_df$min_dist = sapply(candidates, function(x) {
      min(abs(x - selected))
    })
    
    # Candidate with largest minimum distance
    best_idx = which.max(candidate_df$min_dist)
    
    best_row = candidate_df[best_idx,
                            c("Concentration",
                              "Response_1",
                              "Response_2",
                              "SD_1",
                              "SD_2")]
    
    combined = bind_rows(combined, best_row)
    
    selected_concs = c(selected_concs, best_row$Concentration)
    
    candidate_df = candidate_df[-best_idx, ]
  }
  
  combined$Measured = combined$Concentration %in% df1$Concentration
  
  ##########Section 3##########
  
  models = build_models(combined$Concentration, combined$Response_2, combined$SD_2, fixed1 = F)
  
  lowresp = min(combined$Response_2)
  hiresp = max(combined$Response_2)
  
  lowSD = min(combined$SD_2)
  hiSD = max(combined$SD_2)
  
  
  start_list = list(
    c = lowresp,
    d = hiresp,
    e = mean(combined$Response_1),
    b = 1
  )
  
  lower_bounds = list(
    c = 0,
    d = 0,
    e = 1e-6,
    b = 1e-6
  )
  
  upper_bounds = list(
    c = lowresp + 2 * lowSD,
    d = 2 * (hiresp + hiSD),
    e = Inf,
    b = Inf
  )
  
  fits = fit_models_ODR(
    x = combined$Response_1,
    y = combined$Response_2,
    SD_x = combined$SD_1,
    SD_y = combined$SD_2,
    models = models,
    start_list = start_list,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    fixed1 = F,
    biasCorrBG = NULL,
    positive_x = F)
  
  AICs = calc_AICs(fits, combined$Response_2)
  print(AICs)
  
  bestmodel = names(AICs[which(AICs == min(AICs))])
  
  cat("\nBest model:", bestmodel, "\n")
  cat("Best-fit parameters:\n")
  print(coef(fits[[bestmodel]]$fit))
  
  ##########Section 4##########
  
  plot_RR_fit(df = combined,
              fit = fits[[bestmodel]]$fit,
              model_fun = get(bestmodel),
              main = main,
              xlab = xlab,
              ylab = ylab,
              sqrtSDplot = sqrtSDplot,
              log_x = log_x)
  
  ggsave(paste0(savepath, savename, '.png'), width = 4, height = 3)
  
  # Optional zoomed plot
  if (!is.null(zoom)) {
    zoom_conc = sort(combined$Concentration)[1:zoom]
    
    combined_zoom = combined[
      combined$Concentration %in% zoom_conc,]
    
    plot_RR_fit(df = combined_zoom,
                fit = fits[[bestmodel]]$fit,
                model_fun = get(bestmodel),
                main = main,
                xlab = xlab,
                ylab = ylab,
                sqrtSDplot = sqrtSDplot,
                log_x = log_x,
                zoom = T)
    
    ggsave(paste0(savepath, savename, '_zoom.png'), width = 3, height = 2.25)
  }
}
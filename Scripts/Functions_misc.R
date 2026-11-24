#Documentation on miscellaneous functions used for the work described in 'Advancing
#AOPs for genotoxicity through quantifying KERs linking DNA alkylation to increase
#in mutations and chromosome aberrations' by Winkelman et al., 2026.

#make_pie:
#  Goal:                  Make a pie chart of extracted versus not extracted studies for one compound.
#Arguments:
#- compound_name:         name of the compound to plot
#- data:                  dataframe with columns 'compound', 'category' ("extracted" or
#                         "not extracted") and 'count'
#Note:                    The legend shows the
#                         counts per category. Returns a ggplot object.


#make_pie2:
#  Goal:                  Make a pie chart of PubMed versus PMC records for one compound.
#Arguments:
#- compound_name:         name of the compound to plot
#- data:                  dataframe with columns 'compound', 'category' ("PubMed" or
#                         "PMC") and 'count'
#Note:                    The legend shows the counts per
#                         category. Returns a ggplot object.


#summarize_all_factors:
#  Goal:                  Summarize the number of rows, studies and DOIs per level of the
#                         main experimental factors.
#Arguments:
#- df:                    input dataframe
#Note:                    Columns 'Study index', 'DOI', 'S9', 'Scoring technique',
#                         'Time of exposure [h]' and 'expression_period' are added as NA if
#                         missing. Expression period (days) is grouped into I (<= 1/24),
#                         E (between 1/24 and 7) and L (>= 7); exposure time is grouped
#                         into S (<= 4 h) and L (> 4 h). Returns a list of summary
#                         dataframes (S9, Scoring_technique, expression_group and
#                         exposure_group), each with n_rows, n_studies and n_dois per
#                         level (missing values shown as "NA").


#summarize_all_facs_specific:
#  Goal:                  Summarize the number of rows, studies and DOIs per level of
#                         user-specified factors.
#Arguments:
#- df:                    input dataframe
#- required_cols:         character vector of the column names to summarize
#Note:                    Adapted version of summarize_all_factors() for assay-specific
#                         variables. Requested columns missing from df are added as NA.
#                         Returns a named list with, for each column, a dataframe with
#                         n_rows, n_studies and n_dois per level (missing values shown
#                         as "NA"). Requires the columns 'Study index' and 'DOI'.


#extract_variance_components:
#  Goal:                  Extract the unique, joint and residual explained variance of a fitted model.
#Arguments:
#- var_obj:               list containing the fitted mixed model ('model') and a dataframe
#                         of R2 estimates ('R2', columns 'term' and 'estimate') for the full
#                         model and for the individual and combined predictors
#                         (combinations denoted with '+')
#- dataset_name:          name of the dataset/assay, stored in the output
#- compound_name:         name of the compound, stored in the output
#Note:                    The unique variance of each predictor is taken from its
#                         single-term R2 entry. Joint explained variance is the full-model
#                         R2 minus the sum of the unique variances (at least 0) and
#                         residual variance is 1 minus the explained variance. The number
#                         of studies is taken from the grouping factor containing 'Study.index'.
#                         Returns a dataframe with columns compound, dataset, component,
#                         variance and n_studies.


#plot_assay:
#  Goal:                  Plot stacked bars of the variance components per compound for one assay.
#Arguments:
#- df:                    dataframe with columns 'dataset', 'compound', 'n_studies',
#                         'component_grouped' and 'percent'
#- assay_name:            name of the assay (value in 'dataset') to plot
#Note:                    Requires a named color vector called 'component_colors' in the
#                         environment (names matching 'component_grouped'). Compound labels
#                         include the number of studies. Only components present in the
#                         assay are shown in the legend. Returns a ggplot object.


#compute_shapley:
#  Goal:                  Calculate Shapley values for the importance of each predictor.
#Arguments:
#- var_obj:               list containing a dataframe of R2 estimates ('R2', columns 'term'
#                         and 'estimate') for the full model and for all subsets of
#                         predictors (subsets denoted with '+')
#Note:                    Shapley values are the average marginal contribution of each
#                         predictor to the R2 over all possible predictor orderings, and
#                         therefore sum to the full-model R2. All subset models must be
#                         present in var_obj$R2, otherwise an error is raised. Negative
#                         values are set to 0 with a warning. Returns a tibble with
#                         columns 'predictor', 'shapley' and 'unique' (R2 of the
#                         predictor alone).


#compute_ir2:
#  Goal:                  Calculate the inclusive R2 (IR2) of each predictor of a mixed model.
#Arguments:
#- var_obj:               list containing the fitted mixed model ('model') and a dataframe
#                         of R2 estimates ('R2') of which the first value is the total R2
#Note:                    For every model term the contribution to the fixed-effects linear
#                         predictor is calculated; the IR2 is the squared correlation
#                         (structure coefficient) between this contribution and the full
#                         fixed-effects prediction, multiplied by the total R2.
#                         Random effects are excluded. Returns a tibble with
#                         columns 'predictor' and 'IR2'.


#extract_importance_metrics:
#  Goal:                  Combine predictor importance metrics (Shapley and/or IR2) for one model.
#Arguments:
#- var_obj:               list containing the fitted model and R2 estimates (see compute_shapley() and compute_ir2())
#- dataset_name:          name of the dataset/assay, stored in the output
#- compound_name:         name of the compound, stored in the output
#- use_shapley:           whether to calculate Shapley values (default: TRUE)
#- use_ir2:               whether to calculate IR2 (default: TRUE)
#Note:                    Returns a dataframe with one row per predictor, the requested
#                         metrics and the columns 'dataset' and 'compound'. At least one
#                         of use_shapley and use_ir2 should be TRUE.


#make_compound_plot:
#  Goal:                  Plot predictor importance (variance explained) for one compound, faceted by assay.
#Arguments:
#- compound_name:         name of the compound to plot
#- plot_imp               dataframe with columns compound, dataset, predictor_label, ordering_value, score and metric
#Note:                    Predictors are ordered within each panel by
#                         'ordering_value'; a trailing ' assay' is removed from the panel
#                         titles. Requires tidytext. Returns a ggplot object.


#combine_var_lists:
#  Goal:                  Combine variance-component results of multiple compounds and assays into one dataframe.
#Arguments:
#- allVARs:               nested named list: compound -> assay -> list containing a
#                         dataframe 'vc' of variance components (row names are the components)
#Note:                    Assays without 'vc' are skipped. A 'prop' column is added with NA
#                         if absent. The compound, assay and component (taken from the row
#                         names) are added as columns. Returns one dataframe (via map_dfr).


#plot_dose_response:
#  Goal:                  Plot dose-response data (means with error bars) on a log10 concentration axis.
#Arguments:
#- df:                    input dataframe
#- conc_col:              name of the concentration column
#- res_col:               name of the response column
#- sem_col:               name of the SEM column
#- unit_col:              name of the column with the concentration unit (first value used in the x-axis label)
#- color_col:             optional name of the variable used to color the data; levels are
#                         labelled by number in order of appearance (default: NULL; black)
#- custom_colors:         optional vector of colors for the study labels (default: NULL)
#- title:                 plot title (default: "")
#- ylim:                  optional y-axis limits, applied by zooming (default: NULL)
#- logtransform:          whether to plot log-transformed responses (default: FALSE)
#- sd_col:                name of the SD column (required if logtransform = TRUE)
#- rep_col:               name of the replicate-count column (required if logtransform = TRUE)
#Note:                    Zero concentrations are displayed at an offset of a quarter of
#                         the lowest non-zero concentration. If logtransform = TRUE,
#                         the response is log-transformed with a bias correction based on
#                         a linear proxy SD model (SD as function of concentration, fitted
#                         with nls; the intercept is the minimum of the first SD and 1.5
#                         times the lowest of the first three SDs), and the SEM is
#                         calculated on the log scale from the coefficient of variation and
#                         replicate count; the y-axis is then labelled 'log(FC)'. The
#                         data are assumed to be ordered by increasing concentration.
#                         Returns a ggplot object.





make_pie = function(compound_name, data) {
  
  plot_data = data %>%
    filter(compound == compound_name) %>%
    mutate(
      legend_label = paste0(
        ifelse(category == "extracted", "Extracted", "Not extracted"),
        " (", count, ")"
      )
    )
  
  colors = c(
    extracted = "steelblue",
    not_extracted = "grey80"
  )
  
  ggplot(plot_data, aes(x = "", y = count, fill = category)) +
    
    geom_col(width = 1, color = "black") +
    
    coord_polar(theta = "y") +
    
    labs(
      title = compound_name,
      fill = NULL
    ) +
    
    scale_fill_manual(
      values = colors,
      labels = plot_data$legend_label
    ) +
    
    theme_void() +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      legend.text = element_text(size = 12)
    )
}



make_pie2 = function(compound_name, data) {
  
  plot_data = data %>%
    filter(compound == compound_name) %>%
    mutate(
      category = factor(category, levels = c("PubMed", "PMC")),
      legend_label = paste0(
        ifelse(category == "PubMed", "PubMed", "PMC"),
        " (", count, ")"
      )
    )
  
  colors = c(
    PubMed = "steelblue",
    PMC = "grey80"
  )
  
  label_vec = setNames(
    plot_data$legend_label,
    as.character(plot_data$category)
  )
  
  ggplot(plot_data, aes(x = "", y = count, fill = category)) +
    
    geom_col(width = 1, color = "black") +
    
    coord_polar(theta = "y") +
    
    labs(
      title = compound_name,
      fill = NULL
    ) +
    
    
    scale_fill_manual(
      values = colors,
      breaks = c("PubMed", "PMC"),
      labels = label_vec[c("PubMed", "PMC")]
    ) +
    
    theme_void() +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      legend.text = element_text(size = 12)
    )
}

#Create function that extracts all required information
summarize_all_factors = function(df) {
  
  # --- Ensure required columns exist ---
  required_cols = c("Study index", "DOI", "S9", "Scoring technique",
                    "Time of exposure [h]", "expression_period")
  
  for (col in required_cols) {
    if (!col %in% names(df)) {
      df[[col]] = NA
    }
  }
  
  # --- Create grouped variables ---
  df2 = df %>%
    mutate(
      expr_group = case_when(
        is.na(expression_period) ~ NA_character_,
        expression_period <= 1/24 ~ "I",
        expression_period > 1/24 & expression_period < 7 ~ "E",
        expression_period >= 7 ~ "L"
      ),
      exposure_group = case_when(
        is.na(`Time of exposure [h]`) ~ NA_character_,
        `Time of exposure [h]` <= 4 ~ "S",
        `Time of exposure [h]` > 4 ~ "L"
      )
    )
  
  # --- Helper summariser ---
  summarise_factor = function(data, var) {
    data %>%
      mutate(var = as.character({{ var }})) %>%
      group_by(var, .drop = FALSE) %>%
      summarise(
        n_rows = n(),
        n_studies = n_distinct(`Study index`, na.rm = TRUE),
        n_dois = n_distinct(DOI, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(var = ifelse(is.na(var), "NA", var)) %>%
      arrange(var)
  }
  
  # --- Compute all summaries ---
  list(
    S9 = summarise_factor(df2, S9),
    Scoring_technique = summarise_factor(df2, `Scoring technique`),
    expression_group = summarise_factor(df2, expr_group),
    exposure_group = summarise_factor(df2, exposure_group)
  )
}

#Adapted function for assay specific variables
summarize_all_facs_specific = function(df, required_cols) {
  
  # --- Ensure required columns exist ---
  for (col in required_cols) {
    if (!col %in% names(df)) {
      df[[col]] = NA
    }
  }
  
  # --- Helper summariser ---
  summarise_factor = function(data, var) {
    data %>%
      mutate(var = as.character({{ var }})) %>%
      group_by(var, .drop = FALSE) %>%
      summarise(
        n_rows = n(),
        n_studies = n_distinct(`Study index`, na.rm = TRUE),
        n_dois = n_distinct(DOI, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(var = ifelse(is.na(var), "NA", var)) %>%
      arrange(var)
  }
  
  # --- Apply to all requested columns ---
  results = lapply(required_cols, function(col) {
    summarise_factor(df, .data[[col]])
  })
  
  names(results) = required_cols
  
  return(results)
}

extract_variance_components = function(var_obj, dataset_name, compound_name) {
  
  model = var_obj$model
  r2_df = var_obj$R2
  
  # ---- Study count (robust) ----
  flist = getME(model, "flist")
  study_term = names(flist)[grepl("Study.index", names(flist))]
  
  n_studies = if (length(study_term) == 0) NA else {
    length(unique(flist[[study_term[1]]]))
  }
  
  # ---- Variance extraction ----
  full_r2 = r2_df %>%
    filter(term == "Full") %>%
    pull(estimate)
  
  unique_df = r2_df %>%
    filter(term != "Full", !str_detect(term, "\\+")) %>%
    transmute(
      compound = compound_name,
      dataset = dataset_name,
      component = term,
      variance = estimate,
      n_studies = n_studies
    )
  
  unique_sum = sum(unique_df$variance, na.rm = TRUE)
  joint_variance = max(0, full_r2 - unique_sum)
  residual_variance = 1 - unique_sum - joint_variance
  
  bind_rows(
    unique_df,
    
    tibble(
      compound = compound_name,
      dataset = dataset_name,
      component = "Joint explained variance",
      variance = joint_variance,
      n_studies = n_studies
    ),
    
    tibble(
      compound = compound_name,
      dataset = dataset_name,
      component = "Residual variance",
      variance = residual_variance,
      n_studies = n_studies
    )
  )
}

plot_assay = function(df, assay_name, component_colors){
  
  assay_df = df %>%
    filter(dataset == assay_name) %>%
    mutate(compound_label = paste0(compound, " (n=", n_studies, ")"))
  
  # Keep only colors actually present
  present = unique(as.character(assay_df$component_grouped))
  cols = component_colors[names(component_colors) %in% present]
  
  assay_df$component_grouped = factor(
    as.character(assay_df$component_grouped),
    levels = names(cols)
  )
  
  assay_df$compound_label = factor(
    assay_df$compound_label,
    levels = unique(assay_df$compound_label)
  )
  
  ggplot(assay_df, aes(x = compound_label, y = percent, fill = component_grouped)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols, drop = TRUE,
                      guide = guide_legend(reverse = TRUE)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = assay_name, x = NULL,
         y = "Percent of total variance", fill = "Variance component") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

compute_shapley = function(var_obj) {
  
  r2_df = var_obj$R2 %>%
    dplyr::select(term, estimate)
  
  predictors = r2_df %>%
    dplyr::filter(term != "Full", !stringr::str_detect(term, fixed("+"))) %>%
    dplyr::pull(term)
  
  n = length(predictors)
  
  canonicalize_term = function(term_vec) {
    term_vec = stringr::str_trim(term_vec)
    term_vec = term_vec[term_vec != ""]
    if (length(term_vec) == 0) return("")
    paste(sort(term_vec), collapse = "+")
  }
  
  r2_lookup = r2_df %>%
    dplyr::mutate(
      term_clean = dplyr::case_when(
        term == "Full" ~ canonicalize_term(predictors),
        TRUE ~ purrr::map_chr(stringr::str_split(term, fixed("+")), canonicalize_term)
      )
    ) %>%
    dplyr::distinct(term_clean, .keep_all = TRUE)
  
  value_map = setNames(r2_lookup$estimate, r2_lookup$term_clean)
  
  get_value = function(term_name) {
    if (term_name == "") return(0)
    if (!(term_name %in% names(value_map))) {
      stop(paste0("Missing subset model: ", term_name))
    }
    unname(value_map[term_name])
  }
  
  shapley_for_predictor = function(pred) {
    others = setdiff(predictors, pred)
    phi = 0
    
    for (k in 0:length(others)) {
      subset_list = combn(others, k, simplify = FALSE)
      
      for (S in subset_list) {
        weight = factorial(length(S)) * factorial(n - length(S) - 1) / factorial(n)
        marginal = get_value(canonicalize_term(c(S, pred))) -
          get_value(canonicalize_term(S))
        phi = phi + weight * marginal
      }
    }
    
    phi
  }
  
  shapley_vals = purrr::map_dbl(predictors, shapley_for_predictor)
  
  # handle negatives centrally
  if (any(shapley_vals < 0)) {
    warning("Negative Shapley values detected (set to 0)")
  }
  
  tibble::tibble(
    predictor = predictors,
    shapley = pmax(0, shapley_vals),
    unique = purrr::map_dbl(predictors, get_value)
  )
}

compute_ir2 = function(var_obj) {
  
  model = var_obj$model
  R2_total = as.numeric(var_obj$R2$estimate)[1]
  
  eta = as.numeric(predict(model, re.form = NA))
  X = model.matrix(model)
  beta = lme4::fixef(model)
  
  assign_vec = attr(X, "assign")
  
  term_labels_full = attr(terms(model), "term.labels")
  term_map = c("(Intercept)", term_labels_full)
  
  term_index = split(seq_along(assign_vec), assign_vec)
  term_index = term_index[names(term_index) != "0"]
  
  groups = term_index
  names(groups) = term_map[as.numeric(names(term_index)) + 1]
  
  IR2 = sapply(names(groups), function(term) {
    cols = groups[[term]]
    if (length(cols) == 0) return(NA_real_)
    
    eta_term = as.numeric(X[, cols, drop = FALSE] %*% beta[cols])
    SC = cor(eta_term, eta)
    
    (SC^2) * R2_total
  })
  
  tibble::tibble(
    predictor = names(groups),
    IR2 = IR2
  )
}

extract_importance_metrics = function(var_obj,
                                      dataset_name,
                                      compound_name,
                                      use_shapley = TRUE,
                                      use_ir2 = TRUE) {
  
  out = NULL
  
  if (use_shapley) {
    shap = compute_shapley(var_obj)
    out = shap
  }
  
  if (use_ir2) {
    ir2 = compute_ir2(var_obj)
    
    out = if (is.null(out)) {
      ir2
    } else {
      dplyr::left_join(out, ir2, by = "predictor")
    }
  }
  
  out %>%
    dplyr::mutate(
      dataset = dataset_name,
      compound = compound_name
    )
}

make_compound_plot = function(compound_name, plot_imp) {
  
  df_sub = plot_imp %>%
    filter(compound == compound_name) %>%
    mutate(
      predictor_in_panel = tidytext::reorder_within(
        predictor_label,
        -ordering_value,
        dataset
      )
    )
  
  component_colors = c(
    "Unique" = "#4C78A8",
    "IR2" = "darkred",
    "Shapley" = "#F58518"
  )
  
  y_lim = c(0, y_max * 1.05)
  
  ggplot(
    df_sub,
    aes(x = predictor_in_panel, y = score, fill = metric)
  ) +
    geom_col(
      width = 0.7,
      position = position_dodge(width = 0.75)
    ) +
    facet_grid(. ~ dataset, scales = "free_x", space = "free_x",
               labeller = labeller(dataset = function(x) sub(' assay$', '', x))) +
    scale_x_reordered() +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    scale_fill_manual(values = component_colors) +
    labs(
      title = compound_name,
      x = NULL,
      y = "Variance explained (%)",
      fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      strip.background = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

combine_var_lists = function(allVARs) {
  
  map_dfr(names(allVARs), function(compound_name) {
    
    comp_list = allVARs[[compound_name]]
    
    map_dfr(names(comp_list), function(assay_name) {
      
      assay_obj = comp_list[[assay_name]]
      
      if (!"vc" %in% names(assay_obj)) {
        return(NULL)
      }
      
      vc_df = assay_obj$vc
      
      # Ensure required column exists
      if (!"prop" %in% colnames(vc_df)) {
        vc_df$prop = NA_real_
      }
      
      vc_df %>%
        mutate(
          compound = compound_name,
          assay = assay_name,
          component = rownames(vc_df)
        )
      
    })
    
  })
}

plot_dose_response = function(
    df,
    conc_col,
    res_col,
    sem_col,
    unit_col,
    color_col = NULL,
    custom_colors = NULL,
    title = "",
    ylim = NULL,
    logtransform = F,
    sd_col = NULL,
    rep_col = NULL
){
  
  df = df  # avoid modifying parent env
  
  dose = df[[conc_col]]
  
  # ============================================================
  # LOG TRANSFORMATION (optional, self-contained)
  # ============================================================
  
  if(logtransform){
    
    if(is.null(sd_col) || is.null(rep_col)){
      stop("logtransform = TRUE requires sd_col and rep_col")
    }
    
    responseUT = df[[res_col]]
    SDUT       = df[[sd_col]]
    reps       = df[[rep_col]]
    
    # --- proxy SD model (simple linear like original) ---
    icept = min(SDUT[1], min(SDUT[1:3] * 1.5))
    
    varmodel = nls(
      SDUT ~ icept + slope * dose,
      start = list(slope = 0),
      control = nls.control(warnOnly = TRUE)
    )
    
    slope = coef(varmodel)[["slope"]]
    proxySD = icept + slope * dose
    
    # --- CV-based log SD ---
    CV2 = (SDUT / responseUT)^2
    logSD = sqrt(log(1 + CV2))
    
    SEM = logSD / sqrt(reps)
    
    biasCorr = -0.5 * log(1 + (proxySD / responseUT)^2)
    
    response = log(responseUT) + biasCorr
    
    df[[res_col]] = response
    df[[sem_col]]  = SEM
    
  } else {
    
    response = df[[res_col]]
    SEM      = df[[sem_col]]
  }
  
  # ============================================================
  # Sequential study labels
  # ============================================================
  
  if(!is.null(color_col)){
    
    unique_groups = unique(df[[color_col]])
    
    df$study_label = factor(
      df[[color_col]],
      levels = unique_groups,
      labels = seq_along(unique_groups)
    )
  }
  
  # ============================================================
  # Handle zero dose on log scale
  # ============================================================
  
  min_nonzero = min(dose[dose > 0], na.rm = TRUE)
  offset_zero = min_nonzero / 4
  
  df$dose_plot = ifelse(dose == 0, offset_zero, dose)
  
  breaks_vec = c(
    offset_zero,
    10^seq(
      floor(log10(min_nonzero)),
      ceiling(log10(max(dose))),
      by = 1
    )
  )
  
  # ============================================================
  # Base plot
  # ============================================================
  
  p = ggplot(df)
  
  # ============================================================
  # Colored version
  # ============================================================
  
  if(!is.null(color_col)){
    
    p = p +
      
      geom_errorbar(
        aes(
          x = dose_plot,
          ymin = .data[[res_col]] - .data[[sem_col]],
          ymax = .data[[res_col]] + .data[[sem_col]],
          color = study_label
        ),
        width = 0.1,
        alpha = 0.7
      ) +
      
      geom_point(
        aes(
          x = dose_plot,
          y = .data[[res_col]],
          color = study_label
        ),
        size = 2
      ) +
      
      labs(color = "Study")
    
    if(!is.null(custom_colors)){
      p = p + scale_color_manual(values = custom_colors)
    }
    
  } else {
    
    # ==========================================================
    # Black version
    # ==========================================================
    
    p = p +
      
      geom_errorbar(
        aes(
          x = dose_plot,
          ymin = .data[[res_col]] - .data[[sem_col]],
          ymax = .data[[res_col]] + .data[[sem_col]]
        ),
        width = 0.1,
        color = "black"
      ) +
      
      geom_point(
        aes(
          x = dose_plot,
          y = .data[[res_col]]
        ),
        size = 2,
        color = "black"
      )
  }
  
  # ============================================================
  # Formatting
  # ============================================================
  
  p = p +
    
    scale_x_continuous(
      trans = "log10",
      breaks = breaks_vec,
      labels = c("0", formatC(breaks_vec[-1], format = "fg"))
    ) +
    
    labs(
      x = paste0(
        gsub("_", " ", conc_col),
        " (", df[[unit_col]][1], ")"
      ),
      y = ifelse(logtransform,'log(FC)', 'FC'),
      title = title
    ) +
    
    theme_classic() +
    
    theme(plot.title = element_text(hjust = 0.5))
  
  # ============================================================
  # Optional y-limits
  # ============================================================
  
  if(!is.null(ylim)){
    p = p + coord_cartesian(ylim = ylim)
  }
  
  return(p)
}
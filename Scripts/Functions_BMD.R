# Documentation on the functions used for benchmark dose (BMD) analysis as
# described in 'Advancing AOPs for genotoxicity through quantifying KERs
# linking DNA alkylation to increase in mutations and chromosome aberrations'
# by Winkelman et al., 2026.


#plot_residual_diagnostics:
#  Goal:                  Assess the distribution of model residuals using diagnostic plots and a normality test.
#Arguments:
#- residuals:             numeric vector of model residuals
#Note:                    Returns a list containing a histogram with an overlaid normal
#                         distribution, a normal QQ plot and the p-value from a
#                         Shapiro-Wilk test.


#Dose-response functions (HM4, HM3, LN4, LN3, EXP4, EXP3, GM4, GM3):
#  Goal:                  Define the dose-response functions used for BMD modelling.
#Arguments:
#- params:                named vector/list of model parameters (b, c, d and/or e, depending on the model)
#- dose:                  numeric vector of dose/concentration values
#Note:                    HM denotes Hill, LN denotes log-normal, EXP denotes exponential
#                         and GM denotes Gompertz. Models ending in 4 include
#                         four parameters (b, c, d, e); models ending in 3 omit the
#                         upper-response parameter d. Each function returns the predicted
#                         response at the supplied doses.


#build_models:
#  Goal:                  Assemble the dose-response models, parameter starting values and fitting bounds.
#Arguments:
#- dose:                  numeric vector of dose/concentration values
#- response:              numeric vector of responses used for model fitting
#- SEM:                   numeric vector of standard errors of the mean
#- model_set:             model set to use (currently only 'default' is supported but may be expanded upon in the future)
#- fixed1:                whether to fix the background parameter c (default: TRUE)
#- biasCorrBG:            optional background bias correction used when c is fixed (default: NULL)
#Note:                    Returns a named list containing the model function, parameter names,
#                         starting values and lower/upper bounds for each model. When
#                         fixed1 = TRUE, c is removed from the fitted parameter set and
#                         fixed to the supplied bias correction.


#weighted_residuals:
#  Goal:                  Calculate residuals weighted by the supplied SEM values for model fitting.
#Arguments:
#- params:                parameter values supplied to the model function
#- dose:                  numeric vector of dose/concentration values
#- response:              numeric vector of observed responses
#- SEM:                   numeric vector of standard errors of the mean
#- model_fun:             dose-response model function
#Note:                    Applies inverse-variance weights, with a floor on the weights
#                         based on the median SEM to limit the influence of observations
#                         with exceptionally small SEM values. Returns the square-root-
#                         weighted residuals.


#get_residuals:
#  Goal:                  Calculate standardized residuals for a fitted model.
#Arguments:
#- fit_obj:               fitted model object
#- dose:                  numeric vector of dose/concentration values
#- response:              numeric vector of observed responses
#- SEM:                   numeric vector of standard errors of the mean
#- model_fun:             dose-response model function
#Note:                    Returns (observed - predicted) / SEM. Returns NA if the fit
#                         object inherits from "try-error".


#calc_AICs:
#  Goal:                  Calculate an AIC score for a collection of fitted models.
#Arguments:
#- fits:                  named list of fitted model objects (or objects containing a fit element)
#- response:              numeric vector of responses used for fitting
#Note:                    Calculates n * log(RSS / n) + 2 * k, where RSS is the residual
#                         sum of squares, n is the number of responses and k is the number
#                         of fitted coefficients. Failed or missing fits receive NA.


#calc_weights:
#  Goal:                  Convert model AIC values into normalized model weights.
#Arguments:
#- AICs:                  numeric vector of AIC values, preferably named by model
#Note:                    Calculates weights proportional to exp(-0.5 * delta AIC) and
#                         returns weights rounded to four decimal places. If one or more
#                         AIC values are -Inf, weight is shared equally among those models.


#plot_fit:
#  Goal:                  Plot observed dose-response data and the fitted dose-response curve.
#Arguments:
#- fit:                   fitted model object containing parameter estimates in $par
#- model_fun:             dose-response model function
#- dose:                  numeric vector of observed dose/concentration values
#- response:              numeric vector of observed responses
#- SEM:                   numeric vector of standard errors of the mean
#- main:                  plot title (default: "")
#- xlab:                  x-axis label (default: "")
#- ylab:                  y-axis label (default: "")
#Note:                    Displays observations with SEM error bars and the fitted curve
#                         on a logarithmic dose axis. Zero-dose observations are shown
#                         at an offset below the lowest non-zero dose; a dashed segment
#                         connects the model predictions at zero and at the curve start.
#                         Returns a ggplot object.


#weighted_prediction:
#  Goal:                  Calculate model-averaged predictions over a dose grid.
#Arguments:
#- dose_grid:             numeric vector of doses at which predictions are required
#- fits:                  named list of fitted model objects
#- models:                named list of model definitions returned by build_models()
#- selected:              names of the models to include in model averaging
#- weights:               vector of weights for selected models
#Note:                    Predicts with each selected model and calculates the weighted
#                         average prediction at every dose in dose_grid.


#fit_models:
#  Goal:                  Fit all supplied dose-response models to the data.
#Arguments:
#- dose:                  numeric vector of dose/concentration values
#- response:              numeric vector of responses used for fitting
#- SEM:                   numeric vector of standard errors of the mean
#- models:                named list of model definitions returned by build_models()
#- fixed1:                whether the background parameter c is fixed (default supplied by caller)
#- biasCorrBG:            optional background bias correction added to c when fixed (default: NULL)
#Note:                    Fits each model by nonlinear least squares using the
#                         Levenberg-Marquardt algorithm. Returns a named list of fits;
#                         unsuccessful fits are stored with fit = NULL.


#get_lognorm_params:
#  Goal:                  Convert an arithmetic mean and standard deviation to lognormal parameters.
#Arguments:
#- m:                     arithmetic mean (must be non-negative)
#- s:                     arithmetic standard deviation
#Note:                    Returns a list containing meanlog and sdlog for the corresponding
#                         lognormal distribution.


#merge_study_data:
#  Goal:                  Merge study-level observations at matching concentrations and summarize them.
#Arguments:
#- TK6_df:                input dataframe containing study, concentration, replicate,
#                         response and variability information
#- conc_col:              name of the concentration column (default: NULL; detected from
#                         "Concentration EMS", "Concentration MMS" or "Concentration ENU")
#Note:                    Maps concentrations across studies using a 1% relative tolerance,
#                         using the study with the largest number of distinct concentrations
#                         as the initial reference grid. Summarizes observations by mapped
#                         concentration and concentration unit, calculating replicate-
#                         weighted responses and relative survival, pooled within-study
#                         and between-study variability, SD and SEM. Uses FC/FC_SD/FC_SEM
#                         when an FC column is present; otherwise uses the treatment outcome
#                         and its SD/SEM. Returns a concentration-level dataframe.


#generate_subsets:
#  Goal:                  Divide an assay dataset into subsets defined by study-design variables.
#Arguments:
#- df:                    input dataframe for one assay
#- assay_name:            assay name used to derive assay-specific grouping variables
#Note:                    Calls derive_group_vars() and splits the data by available,
#                         non-missing exposure, expression, scoring and extra grouping
#                         variables. Returns a tibble containing a subset label, the
#                         corresponding dataframe and its row count for each subset.


#process_subset:
#  Goal:                  Merge and analyse one assay subset using the BMD bootstrap workflow.
#Arguments:
#- df:                    input dataframe for one assay subset
#- assay_name:            assay name; "alkyl" uses the treatment outcome, while other assays use FC
#- conc_col:              name of the concentration column passed to run_bmd_boot() (default: "Concentration_EMS")
#- BMR:                   benchmark response; numeric percentage or "SD" (default: 50)
#- subset_name:           optional subset label used in the plot title (default: NULL)
#- compound:              compound name used in the plot title (default: "EMS")
#Note:                    First combines study data with merge_study_data(), then runs
#                         run_bmd_boot(). Returns the BMD estimate and confidence limits,
#                         together with bootstrap predictions, weights, model selections
#                         and fit objects. Returns NULL if merging yields no usable rows.


#make_global_compare_plot:
#  Goal:                  Compare a single global bootstrap quantity with its original estimate.
#Arguments:
#- df:                    dataframe containing the bootstrap and original estimates
#- boot_col:              name of the column containing bootstrap values
#- true_col:              name of the column containing the original estimate
#- title:                 plot title
#- ylim:                  optional y-axis limits (default: NULL)
#Note:                    Displays the bootstrap distribution as a boxplot, the original
#                         estimate as a red point and the bootstrap mean as a black diamond.
#                         Optional limits are applied with coord_cartesian(). Returns a
#                         ggplot object.


#make_compare_plot:
#  Goal:                  Compare bootstrap predictions and original estimates across doses.
#Arguments:
#- df:                    dataframe containing dose, bootstrap predictions and original estimates
#- boot_col:              name of the column containing bootstrap values
#- true_col:              name of the column containing the original estimate
#- title:                 plot title
#Note:                    Displays a boxplot of bootstrap predictions, the original estimate
#                         as a red point and the bootstrap mean as a black diamond at each
#                         dose. Limits the number of displayed dose labels to ten and
#                         adjusts the y-axis display when values extend beyond -10 to 10.
#                         Returns a ggplot object.


#plot_BMDresults:
#  Goal:                  Plot BMD estimates and their confidence intervals for multiple subsets.
#Arguments:
#- data:                  dataframe containing BMD, BMD_low, BMD_high, label and subset columns
#- y_mode:                y-axis label mode, either "subset" or "label" (default: "subset")
#Note:                    Displays BMD estimates and confidence intervals on a logarithmic
#                         concentration axis. In "subset" mode, the y-axis labels are
#                         replaced by the corresponding subset names. Returns a ggplot object.


#backtransform_mu:
#  Goal:                  Back-transform a bias-corrected log-scale prediction to the response scale.
#Arguments:
#- yhat:                  predicted response on the log scale
#- SD:                    standard deviation parameter used for the bias correction
#- interval:              numeric interval in which to search for the response-scale root
#Note:                    Uses uniroot() with frootResp() to solve for the response-scale
#                         value. Returns NA if no root can be bracketed or root finding fails.


#frootResp:
#  Goal:                  Define the equation used to back-transform a log-scale prediction.
#Arguments:
#- mu:                    candidate response-scale mean
#- yhat:                  predicted response on the log scale
#- SD:                    standard deviation parameter used for the bias correction
#Note:                    Returns log(mu) - 0.5 * (SD / mu)^2 - yhat; a zero of this
#                         expression is the back-transformed response.


#run_bmd_boot:
#  Goal:                  Fit dose-response models, calculate a model-averaged BMD and quantify its uncertainty by bootstrap.
#Arguments:
#- doseres_df:            input dataframe containing dose, response, SEM, SD, concentration unit and replicate count
#- conc_col:              name of the concentration column (default: "Concentration")
#- res_col:               name of the response column (default: "FC")
#- sem_col:               name of the SEM column (default: "SEM")
#- sd_col:                name of the SD column (default: "SD")
#- unit_col:              name of the concentration-unit column (default: "Concentration_unit")
#- bias_summ:             optional filename stem for the bias-summary plot (default: "")
#- rep_col:               name of the replicate-count column (default: "total_replicates")
#- fixed1:                whether to fix the background parameter c (default: TRUE)
#- nboot:                 number of bootstrap iterations; use 0 to skip bootstrap (default: 200)
#- BMR:                   benchmark response, either a numeric percentage or "SD" (default: 50)
#- seed:                  random-number seed (default: 20251208)
#- sampling:              bootstrap sampling distribution setting (default: "lognorm")
#- pboottitle:            title used for plots and their filenames
#- savepath:              base directory for saving diagnostic and bootstrap plots
#Note:                    Estimates a proxy SD trend, transforms the response and SEM to
#                         the log scale with bias correction, fits four families of models
#                         selects the best three- or four-parameter model within each family 
#                         using AIC, and calculates model-averaged predictions. The BMD is 
#                         the dose at which the back-transformed prediction reaches
#                         the specified benchmark response. When nboot > 0, bootstrap 
#                         resampling is used to obtain BMD confidence limits and prediction
#                         distributions. The returned list includes AICs, selected-model
#                         weights, model definitions, BMD and confidence limits, and 
#                         (when available) bootstrap predictions, fits and related intermediate
#                         results. Diagnostics and bootstrap plots are saved under savepath.


#derive_group_vars:
#  Goal:                  Apply preprocessing filters and derive assay-specific subset variables.
#Arguments:
#- df:                    input dataframe for one assay
#- assay_name:            assay name used to determine assay-specific grouping
#Note:                    Filters exposure duration, gene locus, S9 condition and selected
#                         DNA-adduct types when the corresponding columns are present.
#                         Creates exposure, expression and scoring groups, and an additional
#                         group for MN assay type, mutation gene locus or ECDS target region
#                         where applicable. Returns the modified dataframe.


#convert_results:
#  Goal:                  Combine BMD results from multiple processed subsets into one dataframe.
#Arguments:
#- results_list:          list of result objects containing assay, subset, BMD, BMD_low and BMD_high
#- compound_name:         compound name added to each output row
#Note:                    Extracts the compound, assay, subset, BMD estimate and confidence
#                         limits from each result object and combines them into a tibble.


#Documentation on the functions used for BMD analysis as described in 'Advancing AOPs for genotoxicity through quantifying KERs linking DNA alkylation to increase in mutations and chromosome aberrations' by Winkelman et al., 2026.

plot_residual_diagnostics = function(residuals) {
  
  # compute stats
  res_mean = mean(residuals)
  res_sd   = sd(residuals)
  
  # turn into data frame
  df = data.frame(residuals = residuals)
  
  # ---------- Histogram with Normal Fit ---------- -> All works
  breaks = pretty(residuals)
  
  curve_df = data.frame(x = seq(min(breaks) - 1, max(breaks) + 1, length.out = 500))
  
  p_hist = ggplot(df, aes(x = residuals)) +
    geom_histogram(
      aes(y = ..density..),
      breaks = breaks,
      fill = "lightblue",
      color = "white"
    ) +
    stat_function(
      data = curve_df,
      aes(x = x),
      fun = dnorm,
      args = list(mean = res_mean, sd = res_sd),
      color = "red",
      linewidth = 1
    ) +
    labs(
      title = "Histogram of Residuals with Normal Fit",
      x = "Residuals",
      y = "Density"
    ) +
    theme_classic()
  
  # ---------- QQ Plot ----------
  # generate theoretical quantiles and sample quantiles
  qq_df = data.frame(
    theoretical = qqnorm(residuals, plot.it = FALSE)$x,
    sample      = qqnorm(residuals, plot.it = FALSE)$y
  )
  
  p_qq = ggplot(qq_df, aes(x = theoretical, y = sample)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
    labs(
      title = "QQ Plot of Residuals",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_classic()
  
  # ---------- Shapiro-Wilk Test ----------
  shapiro_p = shapiro.test(residuals)$p.value
  
  # return all plots in a list (with p-value as attribute)
  out = list(
    histogram = p_hist,
    qqplot    = p_qq, 
    shapiro_pvalue = shapiro_p
  )
  
  return(out)
}

#Model definitions
HM4 = function(params, dose) {
  with(as.list(params), {
    c + (d - c) * (dose^b) / (e^b + dose^b)
  })
}

HM3 = function(params, dose) {
  with(as.list(params), {
    c + (dose / e)^b
  })
}

LN4 = function(params, dose) {
  with(as.list(params), {
    c + (d - c) * pnorm(b * (log(dose) - log(e)))
  })
}

LN3 = function(params, dose) {
  with(as.list(params), {
    c + -log(1 - pnorm(b * (log(dose) - log(e))))
  })
}

EXP4 = function(params, dose) {
  with(as.list(params), {
    c + (d - c) * (1 - exp(-(dose / e)^b))
  })
}

EXP3 = function(params, dose) {
  with(as.list(params), {
    c + exp((dose / e)^b) - 1
  })
}

GM4 = function(params, dose) {
  with(as.list(params), {
    c + (d - c) * exp(-exp(-b * (log(dose) - log(e))))
  })
}

GM3 = function(params, dose) {
  with(as.list(params), {
    c + -log(1 - exp(-exp(-b * (log(dose) - log(e)))))
  })
}

build_models = function(dose, response, SEM, model_set = 'default', fixed1 = TRUE, biasCorrBG = NULL) {
  # Define starting parameters, and lower and upper bounds for fitting.
  lowresp = response[which.min(response)]
  hiresp = response[which.max(response)]
  
  lowSEM = SEM[which.min(SEM)]
  hiSEM = SEM[which.max(SEM)]
  
  # Model definitions default (all parameters and bounds included) for the chosen set of models.
  if(model_set == "default"){
    # Define bounds and starting position for the default models.
    start_list_default = list(
      c = lowresp - lowSEM,
      d = hiresp,
      e = mean(dose),
      b = 1
    )
    
    lower_bounds_default = list(
      c = -Inf,
      d = 0,
      e = 0,
      b = 1e-6
    )
    
    upper_bounds_default = list(
      c = lowresp + 2 * lowSEM,
      d = 2 * (hiresp + hiSEM),
      e = Inf,
      b = Inf
    )
    
    # Set models as the default models with the corresponding bounds and starting positions.
    base_models = list(
      HM4 = list(fun = HM4, pars = c("b","c","d","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default),
      HM3 = list(fun = HM3, pars = c("b","c","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default),
      
      LN4 = list(fun = LN4, pars = c("b","c","d","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default),
      LN3 = list(fun = LN3, pars = c("b","c","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default),
      
      EXP4 = list(fun = EXP4, pars = c("b","c","d","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default),
      EXP3 = list(fun = EXP3, pars = c("b","c","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default),
      
      GM4 = list(fun = GM4, pars = c("b","c","d","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default),
      GM3 = list(fun = GM3, pars = c("b","c","e"), start_list = start_list_default, upper_bounds = upper_bounds_default, lower_bounds = lower_bounds_default)
    )
  } 
  
  
  # If c is fixed, transform each model.
  if (fixed1) {
    models = lapply(base_models, function(m) {
      
      # Remove "c" from the parameter vector.
      new_pars = setdiff(m$pars, "c")
      
      # Wrap the model function to force a fixed background.
      new_fun = function(params, dose) {
        m$fun(c(params, c = 0 + biasCorrBG), dose)
      }
      
      list(fun = new_fun, pars = new_pars, start_list = m$start_list, lower_bounds = m$lower_bounds, upper_bounds = m$upper_bounds)
    })
  } else {
    models = base_models
  }
  return(models)
}



weighted_residuals = function(params, dose, response, SEM, model_fun) {
  
  SEM_floor = median(SEM) * .2 #Avoid huge effect of single obs with tiny SEM
  w_floor = 1/SEM_floor
  
  w = pmin(1 / SEM^2, w_floor)
  resid = response - model_fun(params, dose)
  sqrt(w) * resid  # weighted residual
}


get_residuals = function(fit_obj, dose, response, SEM, model_fun) {
  if (inherits(fit_obj, "try-error")) return(NA)
  par_est = coef(fit_obj)
  pred = model_fun(par_est, dose)
  res = (response - pred) / SEM  # same form as weighted residuals, if applicable
  return(res)
}

calc_AICs = function(fits, response) {
  sapply(fits, function(fit_obj) {
    # If fit_obj failed, NA is returned.
    if (inherits(fit_obj, "try-error") || is.null(fit_obj)) return(NA_real_)
    
    # Handle nested structure like fit$fit if needed.
    fit = if (!is.null(fit_obj$fit)) fit_obj$fit else fit_obj
    
    if (inherits(fit, "try-error")) return(NA_real_)
    
    # Calculate the AIC score.
    RSS = sum(stats::resid(fit)^2)
    k = length(stats::coef(fit))   # number of estimated parameters
    n = length(response)
    
    AIC = n * log(RSS / n) + 2 * k
    return(AIC)
  })
}



calc_weights = function(AICs){
  
  is_inf = is.infinite(AICs) & AICs == -Inf
  
  if (any(is_inf)) {
    n_inf = sum(is_inf)
    w = rep(0, length(AICs))
    w[is_inf] = 1 / n_inf
    return(round(w, 4))
  }
  delts = AICs - min(AICs)
  rawweights = exp(-0.5 * delts)
  weights = round(rawweights / sum(rawweights), 4)
  
  return(weights)
}


plot_fit = function(fit, model_fun, dose, response, SEM, main = '', xlab = '', ylab = '') {
  
  # Determine minimal positive dose and offset positions
  min_nonzero = min(dose[dose > 0])
  offset_zero = min_nonzero / 4    # visual position for 0-dose
  offset_start = min_nonzero / 2   # where the fitted curve begins
  
  # --- Compute model predictions ---
  fitdose_log = exp(seq(log(offset_start), log(max(dose)), length.out = 200))
  fitresp_log = model_fun(fit$par, fitdose_log)
  
  # Also compute model predictions at 0 and offset_start
  fitresp_zero = model_fun(params = fit$par, dose = 0)       
  fitresp_start = model_fun(fit$par, offset_start)
  breaks_vec = c(offset_zero,
                 10^seq(floor(log10(min_nonzero)), ceiling(log10(max(dose))),
                        by = 1))
  
  # --- Data preparation ---
  df_pred = data.frame(dose = fitdose_log, response = fitresp_log)
  df_data = data.frame(dose = dose, response = response, SEM = SEM)
  df_data$dose_plot = ifelse(df_data$dose == 0, offset_zero, df_data$dose)
  
  # --- Plot ---
  p = ggplot() +
    # 1. Error bars
    geom_errorbar(
      data = df_data,
      aes(x = dose_plot, ymin = response - SEM, ymax = response + SEM),
      color = "black",
      width = .1
    ) +
    # 2. Points
    geom_point(
      data = df_data,
      aes(x = dose_plot, y = response),
      size = 2
    ) +
    # 3. Dashed linear segment (from model prediction at 0 to curve start)
    geom_segment(
      aes(
        x = offset_zero,
        y = fitresp_zero,              # =- use model prediction instead of data
        xend = offset_start,
        yend = fitresp_start
      ),
      color = "blue",
      linetype = "dashed"
    ) +
    # 4. Blue model curve
    geom_line(
      data = df_pred,
      aes(x = dose, y = response),
      color = "blue",
      linewidth = 1
    ) +
    # 5. Axes and theme
    scale_x_continuous(
      trans = "log10",
      breaks = breaks_vec,
      labels = c("0", formatC(breaks_vec[-1], format = "fg"))
    ) +
    labs(x = xlab, y = ylab, title = main) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}


weighted_prediction = function(dose_grid, fits, models, selected, weights) {
  
  # For each selected model, compute predictions (DON'T simplify to vector)
  preds = sapply(seq_along(selected), function(i) {
    model_name = selected[i]
    model_fun = models[[model_name]]$fun
    params = fits[[model_name]]$fit$par
    model_fun(params = params, dose = dose_grid)
  }, simplify = FALSE)
  
  # Convert list of vectors into a matrix with each model as a column
  preds = do.call(cbind, preds)
  
  # Weighted average across models
  rowSums(preds * weights[col(preds)])
}


fit_models = function(dose, response, SEM, models, fixed1, biasCorrBG = NULL) {
  fits = lapply(names(models), function(name) {
    m = models[[name]]
    
    start_list = m$start_list
    lower_bounds = m$lower_bounds
    upper_bounds = m$upper_bounds
    
    # Get the starting position, and lower and upper limit for fitting.
    sp = unlist(start_list[m$pars])
    lower = unlist(lower_bounds[m$pars])
    upper = unlist(upper_bounds[m$pars])
    
    # Fit the model using non-linear least squares with the Levenberg-Marquardt algorithm.
    fit = try(
      minpack.lm::nls.lm(
        par = sp,
        fn = function(p){
          weighted_residuals(
            p,
            dose = dose,
            response = response,
            SEM = SEM,
            model_fun = m$fun
          )
        },
        lower = lower,
        upper = upper,
        control = minpack.lm::nls.lm.control(maxiter = 500)
      ),
      silent = TRUE
    )
    
    # If fit failed, keep fit as NULL (so downstream code can detect).
    if(inherits(fit, "try-error") || is.null(fit$par)) {
      fit = NULL
    }
    
    if ((!(is.null(fit))) && fixed1){
      fit$par = c(fit$par, c = 0 + biasCorrBG)
      
      # Keep parameter order consistent.
      fit$par = fit$par[order(names(fit$par))]
    }
    list(name = name, fit = fit)
  })
  
  # Return models as named list of fits (so fits$HM4$fit etc). Returns NULL if fitting was unsuccessful.
  names(fits) = names(models)
  fits_named = stats::setNames(fits, names(models))
  return(fits_named)
}



get_lognorm_params = function(m, s){
  if(any(m < 0)) stop("mean must be > 0 for lognormal")
  sigma2 = log(1 + (s^2) / (m^2))
  sdlog  = sqrt(sigma2)
  meanlog = log(m) - 0.5 * sigma2
  list(meanlog = meanlog, sdlog = sdlog)
}


#Incorporates both within and between group variance (Newest version also takes between study variance of ctrl into account)
merge_study_data = function(TK6_df, conc_col = NULL) {
  
  # Ensure 'replicates' is numeric
  TK6_df$replicates = as.numeric(TK6_df$replicates)
  
  # ---- Detect concentration column ----
  if (is.null(conc_col)) {
    possible_cols = c("Concentration EMS", "Concentration MMS", "Concentration ENU")
    found_cols = possible_cols[possible_cols %in% colnames(TK6_df)]
    
    if (length(found_cols) == 0) {
      stop("No valid concentration column found.")
    } else if (length(found_cols) > 1) {
      stop(paste("Multiple concentration columns found:", paste(found_cols, collapse = ", "),
                 "\nPlease specify conc_col manually."))
    }
    
    conc_col = found_cols[1]
  }
  
  # ---- Decide measurement column ----
  use_FC = "FC" %in% colnames(TK6_df)
  
  # ---- Select SD/SEM columns dynamically ----
  if (use_FC) {
    sd_col  = "FC_SD"
    sem_col = "FC_SEM"
  } else {
    sd_col  = "SD"
    sem_col = "SEM"
  }
  
  # ---- Base columns ----
  keep_cols = c(
    "Study index",
    conc_col,
    "Concentration unit (eg μg/ml or μM)",
    "replicates",
    "Relative survival (%)",
    sd_col,
    sem_col
  )
  
  if (use_FC) {
    keep_cols = c(keep_cols, "FC", "Outcome for the treatment")
  } else if ("Outcome for the treatment" %in% colnames(TK6_df)) {
    keep_cols = c(keep_cols, "Outcome for the treatment")
    message("'FC' column not found. Using 'Outcome for the treatment' instead.")
  } else {
    stop("Neither 'FC' nor 'Outcome for the treatment' column found.")
  }
  
  # Keep only necessary columns
  TK6_df = TK6_df[, keep_cols]
  
  # ---- Controlled concentration mapping ----
  
  # helper tolerance function
  is_match = function(x, ref,
                      rel_tol = 0.01,
                      abs_tol = NULL) {
    
    # near-zero safeguard
    if (is.null(abs_tol)) {
      abs_tol = rel_tol * max(ref, na.rm = TRUE)
    }
    
    if (ref == 0) {
      return(abs(x - ref) <= abs_tol)
    }
    
    abs(x - ref) / abs(ref) <= rel_tol
  }
  
  # ---- determine reference study ----
  
  ref_study = TK6_df %>%
    group_by(`Study index`) %>%
    summarise(
      n_conc = n_distinct(.data[[conc_col]]),
      .groups = "drop"
    ) %>%
    arrange(desc(n_conc)) %>%
    slice(1) %>%
    pull(`Study index`)
  
  # reference concentration grid
  ref_concs = TK6_df %>%
    filter(`Study index` == ref_study) %>%
    pull(.data[[conc_col]]) %>%
    unique() %>%
    sort()
  
  # ---- map concentrations iteratively ----
  
  TK6_df$mapped_conc = TK6_df[[conc_col]]
  
  # Start with reference concentrations
  dynamic_concs = sort(ref_concs)
  
  all_studies = unique(TK6_df$`Study index`)
  
  for (study in all_studies) {
    
    idx = which(TK6_df$`Study index` == study)
    
    for (i in idx) {
      
      current_conc = TK6_df[[conc_col]][i]
      
      # Find candidate matches in CURRENT dynamic grid
      matches = dynamic_concs[
        sapply(dynamic_concs, function(ref) {
          is_match(current_conc, ref,
                   rel_tol = 0.01)
        })
      ]
      
      # Case 1: exactly one match -> map to existing bin
      if (length(matches) == 1) {
        
        TK6_df$mapped_conc[i] = matches
        
        # Case 2: ambiguous matches -> preserve original
      } else if (length(matches) > 1) {
        
        TK6_df$mapped_conc[i] = current_conc
        
        # optionally add to grid
        if (!(current_conc %in% dynamic_concs)) {
          dynamic_concs = sort(c(dynamic_concs, current_conc))
        }
        
        # Case 3: no match -> create NEW bin
      } else {
        
        TK6_df$mapped_conc[i] = current_conc
        
        # Add this concentration to future reference grid
        dynamic_concs = sort(c(dynamic_concs, current_conc))
      }
    }
  }
  
  # replace concentration column
  TK6_df[[conc_col]] = TK6_df$mapped_conc
  TK6_df$mapped_conc = NULL
  
  # ---- Summarise ----
  merged_df = TK6_df %>%
    group_by(.data[[conc_col]], `Concentration unit (eg μg/ml or μM)`) %>%
    summarise(
      total_replicates = sum(replicates, na.rm = TRUE),
      
      weighted_fc = if (use_FC) {
        weighted.mean(FC, replicates, na.rm = TRUE)
      } else {
        weighted.mean(`Outcome for the treatment`, replicates, na.rm = TRUE)
      },
      
      Rel_surv = weighted.mean(`Relative survival (%)`, replicates, na.rm = TRUE),
      
      SD = {
        
        means = if (use_FC) FC else `Outcome for the treatment`
        
        # --- Within-study variance ---
        within_var_numer = sum((replicates - 1) * (.data[[sd_col]]^2), na.rm = TRUE)
        within_var_denom = sum((replicates - 1)[!is.na(.data[[sd_col]])], na.rm = TRUE)
        
        within_var = ifelse(within_var_denom > 0,
                            within_var_numer / within_var_denom,
                            0)
        
        # --- Between-study variance ---
        mean_fc = weighted.mean(means, replicates, na.rm = TRUE)
        
        if (all(.data[[conc_col]] == 0) && use_FC) {
          
          # Special handling for control
          raw_means = TK6_df %>%
            filter(.data[[conc_col]] == 0) %>%
            pull(`Outcome for the treatment`)
          
          if (sum(!is.na(raw_means)) < 2) {
            between_var = 0
          } else {
            between_var = var(raw_means, na.rm = TRUE) /
              (mean(raw_means, na.rm = TRUE)^2)
          }
          
        } else {
          
          between_var = sum(replicates * (means - mean_fc)^2, na.rm = TRUE) /
            sum(replicates, na.rm = TRUE)
        }
        
        sqrt(within_var + between_var)
      },
      
      SEM = SD / sqrt(total_replicates),
      
      sourceStudyIndex = paste(unique(`Study index`[!is.na(`Study index`)]), collapse = ", "),
      
      .groups = "drop"
    ) %>%
    filter(!is.na(SEM) & SEM != 0) %>%
    arrange(.data[[conc_col]]) %>%
    rename(
      Concentration = .data[[conc_col]],
      Concentration_unit = `Concentration unit (eg μg/ml or μM)`
    )
  
  # ---- Rename output ----
  if (use_FC) {
    merged_df = merged_df %>% rename(FC = weighted_fc)
  } else {
    merged_df = merged_df %>% rename(outcome_treatment = weighted_fc)
  }
  
  return(merged_df)
}

generate_subsets = function(df, assay_name) {
  
  df = derive_group_vars(df, assay_name)
  
  # columns that define subsets
  group_cols = c("exposure_group", "expression_group", "scoring_group", "extra_group")
  
  # keep only those that exist and are not all NA
  group_cols = group_cols[group_cols %in% colnames(df)]
  
  group_cols = group_cols[
    sapply(df[group_cols], function(x) !all(is.na(x)))
  ]
  
  
  subsets_tbl = df %>%
    group_by(across(all_of(group_cols))) %>%
    group_split() %>%
    purrr::map_dfr(function(sub_df) {
      
      vals = sub_df[1, group_cols, drop = FALSE]
      
      label = paste0(
        vals$exposure_group,
        vals$expression_group,
        vals$scoring_group,
        if ("extra_group" %in% names(vals)) vals$extra_group else ""
      )
      
      tibble(
        subset_name = label,
        data = list(sub_df),
        n = nrow(sub_df)
      )
    })
  
  subsets_tbl
}

process_subset = function(df, assay_name, conc_col = 'Concentration_EMS', BMR = 50, subset_name = NULL, compound = 'EMS') {
  
  merged = merge_study_data(df)
  
  # Skip if no usable data after merging
  if (is.null(merged) || nrow(merged) == 0) {
    return(NULL)
  }
  
  bmd = run_bmd_boot(
    merged,
    res_col = if (assay_name == "alkyl") "outcome_treatment" else "FC",
    conc_col = conc_col,
    fixed1 = ifelse(assay_name == 'alkyl', F, T),
    BMR = BMR,
    sampling = "lognorm",
    pboottitle = paste(compound, assay_name, subset_name)
  )
  
  bmdres = data.frame(
    BMD_low  = bmd$BMDCI[1],
    BMD_high = bmd$BMDCI[2],
    BMD = bmd$BMD
  )
  
  return(list(bmdres = bmdres, boot_long = bmd$boot_long, wpinput = list(
    weights = bmd$boot_weights, models = bmd$models, best_models = bmd$boot_best_models, fits = bmd$boot_fits, BMD = bmd$boot_BMD, boot_dfs = bmd$boot_dfs, boot_global = bmd$boot_global)))
}

make_global_compare_plot = function(
    df,
    boot_col,
    true_col,
    title,
    ylim = NULL
){
  # The true response is plotted as a red point and the bootstrap as a boxplot and black point per concentration.
  p = ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = "",
      y = .data[[boot_col]]
    )
  ) +
    
    ggplot2::geom_boxplot(
      fill = "lightblue",
      alpha = 0.7,
      outlier.alpha = 0.2
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(y = .data[[true_col]]),
      color = "red",
      size = 3
    ) +
    
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 3,
      color = "black"
    ) +
    
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = title,
    ) +
    
    ggplot2::theme_bw()
  # Apply limits only if necessary.
  if(!is.null(ylim)){
    p = p + ggplot2::coord_cartesian(ylim = ylim)
  }
  
  return(p)
}

make_compare_plot = function(
    df,
    boot_col,
    true_col,
    title
){
  
  # Determine whether limits are needed by checking the bootstrap and true values.
  all_vals = c(
    df[[boot_col]],
    df[[true_col]]
  )
  
  use_limits = any(all_vals < -10 | all_vals > 10)
  
  # Determine x-axis breaks.
  n_labels = 10
  
  dose_levels = unique(df$dose)
  
  if (length(dose_levels) > n_labels) {
    keep = round(seq(1, length(dose_levels), length.out = n_labels))
    breaks = dose_levels[keep]
  } else {
    breaks = dose_levels
  }
  
  # The true response is plotted as a point and the bootstrap as a boxplot per concentration.
  p = ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = factor(.data$dose),
      y = .data[[boot_col]]
    )
  ) +
    
    ggplot2::geom_boxplot(
      fill = "lightblue",
      alpha = 0.7,
      outlier.alpha = 0.2
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(y = .data[[true_col]]),
      color = "red",
      size = 2
    ) +
    
    ggplot2::stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 2.5,
      color = "black"
    ) +
    
    ggplot2::scale_x_discrete(
      breaks = as.character(breaks),
      labels = function(x) round(as.numeric(x))
    ) +
    
    ggplot2::labs(
      x = "Concentration (μM)",
      y = NULL,
      title = title
    ) +
    
    ggplot2::theme_bw()
  
  # Apply limits only if necessary.
  if(use_limits){
    p = p + ggplot2::coord_cartesian(ylim = c(-10, 160))
  }
  
  return(p)
}

plot_BMDresults = function(data, y_mode = c("subset", "label")) {
  
  y_mode = match.arg(y_mode)
  
  ggplot(data, aes(y = label, x = BMD)) +
    
    # whiskers / confidence interval
    geom_errorbar(
      aes(xmin = BMD_low, xmax = BMD_high),
      height = 0.4,
      color = "lightblue",
      linewidth = 1.2,
      orientation = "y"
    ) +
    
    # central estimate
    geom_point(
      color = "darkblue",
      size = 3
    ) +
    
    scale_x_log10() +
    
    labs(
      x = "Concentration (μM)",
      y = ""
    ) +
    
    theme_classic() +
    
    scale_y_discrete(
      labels = if (y_mode == "subset") {
        setNames(data$subset, data$label)
      } else {
        waiver()
      }
    )
}

backtransform_mu = function(yhat, SD, interval) {
  
  f = function(mu) frootResp(mu, yhat = yhat, SD = SD)
  
  vals = try(sapply(interval, f), silent = TRUE)
  
  if (inherits(vals, "try-error") || any(is.na(vals))) {
    return(NA_real_)
  }
  
  # No root in interval
  if (vals[1] * vals[2] > 0) {
    return(NA_real_)
  }
  
  # Safe to solve
  tryCatch({
    uniroot(f, interval = interval)$root
  }, error = function(e) NA_real_)
}

frootResp = function(mu, yhat, SD){
  log(mu) - .5 * (SD / mu)^2 - yhat
}

run_bmd_boot = function(doseres_df,
                        conc_col = 'Concentration',
                        res_col = 'FC',
                        sem_col = 'SEM',
                        sd_col = 'SD',
                        unit_col = 'Concentration_unit',
                        bias_summ = '',
                        rep_col = 'total_replicates',
                        fixed1 = T,
                        nboot = 200,
                        BMR = 50,
                        seed = 20251208,
                        sampling = 'lognorm',
                        pboottitle = 'Bootstrap weighted predictions',
                        savepath = '../Figures/BMD/'){
  
  set.seed(seed)
  
  dose = doseres_df[[conc_col]]
  response = responseUT = doseres_df[[res_col]]
  SEM = SEMUT = doseres_df[[sem_col]]
  SD = SDUT = doseres_df[[sd_col]]
  concunit = doseres_df[[unit_col]]
  reps = doseres_df[[rep_col]]
  
  # Define a name for saving if not specified.
  basesave = paste(gsub(' ', '_', pboottitle), ifelse(fixed1, 'fixed1', ''), sep = '_')
  
  # Set near-zero limit to response to avoid infinities.
  response[response < 10^-8] = 10^-8
  responseUT[responseUT < 10^-8] = 10^-8
  
  # Get proxy standard deviation and bias correction.
  icept = min(SDUT[1], min(SDUT[1:3] * 1.5))
  
  varmodel = stats::nls(SDUT ~ icept + slope * dose,
                        start = list(slope = 0),
                        control = stats::nls.control(maxiter = 500, warnOnly = TRUE))
  
  slope = stats::coef(varmodel)[["slope"]]
  proxySD = icept + slope * dose
  
  CV2 = (SDUT / responseUT)^2
  SD = sqrt(log(1 + CV2))

  SEM = SD / sqrt(reps)
  
  biasCorr = -0.5 * log(1 + (proxySD / responseUT)^2)
  
  #Remove extreme outliers
  keep = abs(biasCorr) <= 10 * mean(abs(biasCorr))
  
  dose       = dose[keep]
  response   = response[keep]
  responseUT = responseUT[keep]
  SEM        = SEM[keep]
  SEMUT      = SEMUT[keep]
  SD         = SD[keep]
  proxySD    = proxySD[keep]
  SDUT       = SDUT[keep]
  reps       = reps[keep]
  biasCorr   = biasCorr[keep]
  doseres_df = doseres_df[keep,]
  
  response = log(responseUT) + biasCorr 
  
  #lowresp = response[which.min(response)]
  #hiresp = response[which.max(response)]
  
  #lowSEM = SEM[which.min(SEM)]
  #hiSEM = SEM[which.max(SEM)]
  
  
  models = build_models(fixed1, biasCorr[1], model_set = 'default', dose = dose, response = response, SEM = SEM)
  
  fits = fit_models(dose, response, SEM, models, fixed1, biasCorr[1])
  
  if (all(sapply(fits, function(x) is.null(x$fit)))){
    stop("Could not fit any models.")
  }
  
  AICs = calc_AICs(fits, response)
  
  cat("AIC per model: \n")
  print(AICs)
  
  xlab = paste0(gsub('_', ' ', conc_col), ' (', doseres_df[[unit_col]][1], ')')
  ylab = paste0('log(', res_col, ')')
  
  p_list = list(
    plot_fit(fits$HM4$fit, HM4, dose, response, SEM, main = paste0('HM4', ' (AIC = ', round(AICs['HM4'],2), ')'), xlab = xlab, ylab = ylab),
    plot_fit(fits$HM3$fit, HM3, dose, response, SEM, main = paste0('HM3', ' (AIC = ', round(AICs['HM3'],2), ')'), xlab = xlab, ylab = ylab),
    plot_fit(fits$LN4$fit, LN4, dose, response, SEM, main = paste0('LN4', ' (AIC = ', round(AICs['LN4'],2), ')'), xlab = xlab, ylab = ylab),
    plot_fit(fits$LN3$fit, LN3, dose, response, SEM, main = paste0('LN3', ' (AIC = ', round(AICs['LN3'],2), ')'), xlab = xlab, ylab = ylab),
    plot_fit(fits$EXP4$fit, EXP4, dose, response, SEM, main = paste0('EXP4', ' (AIC = ', round(AICs['EXP4'],2), ')'), xlab = xlab, ylab = ylab),
    plot_fit(fits$EXP3$fit, EXP3, dose, response, SEM, main = paste0('EXP3', ' (AIC = ', round(AICs['EXP3'],2), ')'), xlab = xlab, ylab = ylab),
    plot_fit(fits$GM4$fit, GM4, dose, response, SEM, main = paste0('GM4', ' (AIC = ', round(AICs['GM4'],2), ')'), xlab = xlab, ylab = ylab),
    plot_fit(fits$GM3$fit, GM3, dose, response, SEM, main = paste0('GM3', ' (AIC = ', round(AICs['GM3'],2), ')'), xlab = xlab, ylab = ylab)
  )
  
  grid_fits = 
    (p_list[[1]] | p_list[[3]] | p_list[[5]] | p_list[[7]]) /
    (p_list[[2]] | p_list[[4]] | p_list[[6]] | p_list[[8]])
  
  # Define a name for saving if not specified.
  basesave = paste(gsub(' ', '_', pboottitle), ifelse(fixed1, 'fixed1', ''), sep = '_')
  
  print(grid_fits)
  ggsave(paste0(savepath, 'fitgrids/', basesave, '_fitgrid.png'), width = 12, height = 6, dpi = 500) 
  
  families = unique(gsub("[34]$", "", names(models)))
  best_models = sapply(families, function(fam) {
    family_models = grep(fam, names(AICs), value = TRUE)
    family_AICs = AICs[family_models]
    best = family_models[which.min(family_AICs)]
    best
  }, simplify = TRUE)
  
  AICs_sel = AICs[best_models]
  
  weights = calc_weights(AICs_sel)
  
  cat("\nWeights of the selected models: \n")
  print(weights)
  
  min_nonzero = min(dose[dose > 0])
  offset_zero = min_nonzero / 4
  offset_start = min_nonzero / 2
  
  fitdose = exp(seq(0, log(max(dose) + 1), length.out = 500)) - 1 
  fitdose = sort(unique(c(dose, offset_zero, offset_start, fitdose)))
  
  wp = weighted_prediction(fitdose, fits, models, best_models, weights)
  
  #Store the best fit for later response-response
  #best_model = names(AICs)[which.min(AICs)]
  #best_fit = fits[[best_model]]
  
  intervalRes = c(0, 50 * max(responseUT))
  
  proxySDFullrange = icept + varmodel$m$getPars() * fitdose
  
  wpBT = numeric(length(fitdose))
  for(i in 1:length(fitdose)){
    wpBT[i] = backtransform_mu(wp[i], proxySDFullrange[i], interval = intervalRes)
  }
  
  # Define background as the lowest weighed prediction back-transformed data point.
  background = wpBT[1]
  
  if(BMR == 'SD'){
    target = background + SEMUT[1] * sqrt(reps[1])
  }else if(is.numeric(BMR)){
    target = background + background * (1/100 * BMR) #Make if statement to support SD based CIs
  }else{
    stop('BMR should be either a number between 1 and 100 or "SD"')
  }
  

  frootBT = function(d) {
    yhat_d = weighted_prediction(d, fits, models, best_models, weights)
    sd_d = icept + varmodel$m$getPars() * d
    mu_d = backtransform_mu(yhat_d, SD = sd_d, interval = intervalRes)
    mu_d - target
  }
  
  d_range = range(fitdose)
  vals = sapply(d_range, frootBT)
  
  # ---- EARLY EXIT ----
  if (any(is.na(vals)) || vals[1] * vals[2] > 0) {
    message(sprintf(
      "No BMD: target not reached (f(min)=%.4f, f(max)=%.4f)",
      vals[1], vals[2]
    ))
    
    return(list(
      AICs   = AICs,
      weights= weights,
      BMD    = NA_real_,
      BMDCI  = c(NA_real_, NA_real_)
    ))
  }
  
  # Find BMD result, and print it.
  BMD = stats::uniroot(frootBT, interval = d_range)$root
  
  cat(paste(" \nResult BMD: ", BMD, "\n", sep=''))
  
  # Stop and return the original BMD results if nboot is zero.
  if (nboot == 0){
    return(list(AICs = AICs, models = models, BMD = BMD))
  }
  
  
  # ---------------------- Bootstrap -----------------------
  
  # Find bootstrapping parameters.
  n_fitdose = length(fitdose)
  boot_BMD = numeric(nboot)
  
  # Configure data structures for bootstrapping results.
  boot_preds = matrix(NA, nrow = nboot, ncol = n_fitdose)
  boot_store = vector("list", nboot)
  boot_global_store = vector("list", nboot)
  boot_weights_store     = vector("list", nboot)
  boot_fits_store        = vector("list", nboot)
  boot_dfs_store          = vector("list", nboot)
  boot_best_models_store = vector("list", nboot)
  
  # 1. Draw a new value nboot times.
  for(i in seq_len(nboot)) {
    list_ress = vector("list", nrow(doseres_df))
    
    # Perform bootstrapping per concentration.
    for (j in seq_len(nrow(doseres_df))) {
      nrep = reps[j]
      mu   = responseUT[j]
      sig  = SDUT[j]
      
      if(sampling == 'lognorm'){
        lp = get_lognorm_params(mu, sig)
        orig_ress = rlnorm(nrep, meanlog = lp$meanlog, sdlog = lp$sdlog)
      }else{
        orig_ress = rnorm(nrep, mu, sig)
      }
      list_ress[[j]] = orig_ress
    }
    
    # 2. Determine original data summaries.
    control_idx = which(dose == min(dose))
    control_mean_boot = mean(unlist(list_ress[control_idx]))
    
    boot_MEAN_orig = numeric(length(dose))
    boot_SD_orig   = numeric(length(dose))
    boot_SEM_orig  = numeric(length(dose))
    
    for (j in seq_along(dose)) {
      
      if(fixed1){
        ress_norm = list_ress[[j]] / control_mean_boot
      } else {
        ress_norm = list_ress[[j]]
      }
      
      boot_MEAN_orig[j] = mean(ress_norm)
      boot_SD_orig[j]   = sd(ress_norm)
      boot_SEM_orig[j]  = sd(ress_norm) / sqrt(length(ress_norm))
    }
    
    ## 3. Identical data transform as the original BMD.
    
    boot_varmodel = try(stats::nls(boot_SD_orig ~ icept + slope * dose,
                                   start = list(slope = 0),
                                   control = stats::nls.control(warnOnly = TRUE)), silent = TRUE)
    
    if(inherits(boot_varmodel, "try-error")) next
    
    boot_slope = coef(boot_varmodel)[["slope"]]
    proxySD_boot = icept + boot_slope * dose
    
    proxySDFullrange = as.numeric(icept + boot_slope * fitdose)
      
    boot_CV2 = (boot_SD_orig / boot_MEAN_orig)^2
    boot_SD = sqrt(log(1 + boot_CV2))
    boot_SEM = boot_SD / sqrt(reps)
    biasCorr_boot = -0.5 * log(1 + (proxySD_boot / boot_MEAN_orig)^2)
    
    boot_response = log(boot_MEAN_orig) + biasCorr_boot
      
    
    boot_df = doseres_df
    boot_df[[res_col]] = boot_response
    boot_df[[sem_col]] = boot_SEM
    boot_df$SD         = boot_SD
    
    # Store bootstrapped dataset
    boot_dfs_store[[i]] = boot_df
    
    ## 4. Fit models
    boot_fits = fit_models(
      dose,
      boot_df[[res_col]],
      boot_df[[sem_col]],
      models, fixed1,
      biasCorr[1]
    )
    
    ## 5. Calculate AICs
    boot_AICs = try(calc_AICs(boot_fits,
                              boot_df[[res_col]]),
                    silent = TRUE)
    
    if(inherits(boot_AICs, "try-error") || any(is.na(boot_AICs))) {
      next
    }
    
    ## 6. Select best model per family
    boot_best_models = sapply(
      families,
      function(fam) {
        family_models = grep(fam, names(boot_AICs), value = TRUE)
        family_AICs   = boot_AICs[family_models]
        family_models[which.min(family_AICs)]
      },
      simplify = TRUE
    )
    
    boot_AICs_sel = boot_AICs[boot_best_models]
    
    boot_weights = try(calc_weights(boot_AICs_sel), silent = TRUE)
    if(inherits(boot_weights, "try-error") || any(is.na(boot_weights))) next
    
    
    ## 7. Weighted prediction on model scale
    
    boot_wp = try(
      weighted_prediction(fitdose,
                          boot_fits,
                          models,
                          boot_best_models,
                          boot_weights),
      silent = TRUE
    )
    if(inherits(boot_wp, "try-error") || length(boot_wp) != n_fitdose) next
    
    
    ## 8. Back-transform weighted predictions if needed
    boot_wpBT = numeric(n_fitdose)
    
    for(k in seq_len(n_fitdose)) {
      
      boot_wpBT[k] = backtransform_mu(
        boot_wp[k],
        SD = proxySDFullrange[k],
        interval = intervalRes
      )
    }
    
    if (any(is.na(boot_wpBT))) {
      boot_BMD[i] = NA_real_
      next   # skip to next bootstrap iteration
    }
    
    background_boot = boot_wpBT[1]
      
    
    ## 9. Compute bootstrap target level
    if(BMR == "SD") {
      target_boot = background_boot +
        boot_SEM_orig[1] * sqrt(boot_df[[rep_col]][1])
      
    } else if(is.numeric(BMR)) {
      target_boot = background_boot + background_boot * (BMR/100)
      
    } else {
      stop('BMR should be either numeric or "SD"')
    }
    
    
    ## 10. Solve for bootstrap BMD
    froot_boot = function(d) {
      
      yhat_d = weighted_prediction(
        d, boot_fits, models, boot_best_models, boot_weights
      )
      
      sd_d = icept + boot_varmodel$m$getPars() * d
      
      mu_d = backtransform_mu(yhat_d, SD = sd_d, interval = intervalRes)
      
      mu_d - target_boot
    }
      
    d_range = c(fitdose[1], tail(fitdose, 1))
    
    # Evaluate function at interval boundaries
    vals = try(sapply(d_range, froot_boot), silent = TRUE)
    
    if (inherits(vals, "try-error") || any(is.na(vals))) {
      
      # true failure (model broke, etc.)
      boot_BMD[i] = NA_real_
      
    } else if (vals[1] * vals[2] > 0) {
      
      # No root in interval → BMD beyond range
      boot_BMD[i] = Inf
      
    } else {
      
      # Safe to solve
      boot_BMD[i] = tryCatch({
        uniroot(froot_boot,
                lower = d_range[1],
                upper = d_range[2])$root
      }, error = function(e) NA_real_)
    }
    
    boot_weights_store[[i]]     = boot_weights
    boot_fits_store[[i]]        = boot_fits
    boot_best_models_store[[i]] = boot_best_models
    
    
    ## 11. Store weighted prediction curve (model scale)
    boot_preds[i, ] = boot_wp
    
    boot_store[[i]] = data.frame(
      iter = i,
      dose = dose,
      
      boot_MEAN_orig = boot_MEAN_orig,
      boot_SEM_orig  = boot_SEM_orig,
      
      boot_response  = boot_response,
      boot_SEM       = boot_SEM,
      
      proxySD_boot   = proxySD_boot,
      biasCorr_boot  = biasCorr_boot,
      
      true_MEAN_orig = responseUT,
      true_SEM_orig  = SEMUT,
      
      true_response  = response,
      true_SEM       = SEM,
      
      true_proxySD   = proxySD,
      true_biasCorr  = biasCorr
    )
    
    boot_global_store[[i]] = data.frame(
      iter = i,
      
      boot_slope = boot_slope,
      true_slope = slope,
      icept = icept,
      
      boot_BMD = boot_BMD[i]
    )
    
  }#end of bootstrap
  
  boot_results_df = dplyr::bind_rows(boot_store)
  boot_global_df  = dplyr::bind_rows(boot_global_store)
  
  if (isTRUE(bias_summ) || is.character(bias_summ)){
    # Compare the bootstrap results to the true values.
    bias_summary = boot_results_df %>%
      dplyr::group_by(dose) %>%
      dplyr::summarise(
        
        # Calculate the bias of the bootstrap estimator on the original scale.
        mean_orig_bias =
          mean(boot_MEAN_orig, na.rm = TRUE) -
          dplyr::first(.data$true_MEAN_orig),
        
        # Calculate the SEM bias of the bootstrap estimator on the original scale.
        sem_orig_bias =
          mean(boot_SEM_orig, na.rm = TRUE) -
          dplyr::first(.data$true_SEM_orig),
        
        # Calculate the response bias of the bootstrap estimator.
        response_bias =
          mean(boot_response, na.rm = TRUE) -
          dplyr::first(.data$true_response),
        
        # Calculate the SEM bias of the bootstrap estimator on a logarithmic scale.
        sem_log_bias =
          mean(boot_SEM, na.rm = TRUE) -
          dplyr::first(.data$true_SEM),
        
        # Calculate the proxy standard deviation bias of the bootstrap estimator.
        proxySD_bias =
          mean(proxySD_boot, na.rm = TRUE) -
          dplyr::first(.data$true_proxySD),
        
        # Calculate the bias correction bias of the bootstrap estimator.
        biasCorr_bias =
          mean(biasCorr_boot, na.rm = TRUE) -
          dplyr::first(.data$true_biasCorr),
        
        .groups = "drop"
      )
    
    # Compare the bootstrap results of the fitted model.
    global_bias_summary = boot_global_df %>%
      dplyr::summarise(
        
        # Calculate the slope bias of the bootstrap estimator.
        slope_bias =
          mean(boot_slope, na.rm = TRUE) -
          dplyr::first(.data$true_slope),
      )
    
    # Make plots comparing the bootstrap results to the true results.
    # 1, Compare original mean.
    p1 = make_compare_plot(
      boot_results_df,
      "boot_MEAN_orig",
      "true_MEAN_orig",
      "Original scale means"
    )
    
    # 2, Compare original SEM
    p2 = make_compare_plot(
      boot_results_df,
      "boot_SEM_orig",
      "true_SEM_orig",
      "Original scale SEMs"
    )
    
    # 3. Compare log-transformed response.
    p3 = make_compare_plot(
      boot_results_df,
      "boot_response",
      "true_response",
      "Log-scale means"
    )
    
    # 4. Compare log-transformed SEM.
    p4 = make_compare_plot(
      boot_results_df,
      "boot_SEM",
      "true_SEM",
      "Log-scale SEMs"
    )
    
    # 5. Compare SD (log-transformed).
    p5 = make_compare_plot(
      boot_results_df,
      "proxySD_boot",
      "true_proxySD",
      "proxySD comparison"
    )
    
    # 6. Compare bias correction.
    p6 = make_compare_plot(
      boot_results_df,
      "biasCorr_boot",
      "true_biasCorr",
      "biasCorr comparison"
    )
    
    # 7. Compare slope.
    slope_plot = make_global_compare_plot(
      boot_global_df,
      "boot_slope",
      "true_slope",
      "Slope comparison"
    )
    
    # Create plot legend.
    legend = grid::grobTree(
      # Red circle with label.
      grid::pointsGrob(x = 0.2, y = 0.55, pch = 16, gp = grid::gpar(col = "red", cex = 1.2)),
      grid::textGrob("Original value", x = 0.25, y = 0.55, just = "left", gp = grid::gpar(fontsize = 11)),
      
      # Black diamond with label.
      grid::pointsGrob(x = 0.2, y = 0.45, pch = 18, gp = grid::gpar(col = "black", cex = 1.2)),
      grid::textGrob("Bootstrap Mean", x = 0.25, y = 0.45, just = "left", gp = grid::gpar(fontsize = 11))
    )
    
    # Combine plots.
    combined_plot =
      (p1 | p2 | p3 | p4) /
      (p5 | p6 | slope_plot | legend)
    
    # Print and save the plot.
    print(combined_plot)
    
    if (is.character(bias_summ)){
      if (nzchar(bias_summ)){
        ggplot2::ggsave(
          paste0(savepath, 'bias_summaries/', bias_summ, '.png'), plot = combined_plot,
          width = 12,
          height = 6,
          dpi = 500,
        )
      } else{
        ggplot2::ggsave(
          paste0(savepath, 'bias_summaries/', basesave, '_bias_summary.png'), plot = combined_plot,
          width = 12,
          height = 6,
          dpi = 500,
        )
      }
    }
  }
  
  
  
  
  # Keep rows where prediction exists AND BMD is not NA
  valid_rows = !is.na(boot_BMD) & apply(!is.na(boot_preds), 1, any)
  
  boot_preds = boot_preds[valid_rows,]
  boot_BMD   = boot_BMD[valid_rows]
  
  nb_valid = nrow(boot_preds)
  
  # Case 1: nothing usable
  if (length(boot_BMD) == 0) {
    
    message("No valid bootstrap BMDs — stopping analysis")
    
    return(list(
      AICs   = AICs,
      weights= weights,
      BMD    = BMD,
      BMDCI  = c(NA_real_, NA_real_)
    ))
  }
  
  # Compute CI (this now includes Inf if present)
  BMDCI = quantile(boot_BMD, c(.05, .95 ), na.rm = TRUE)
  
  print(BMDCI)
  
  
  # Case 2: all infinite -> no information
  if (all(is.infinite(boot_BMD))) {
    
    message("All bootstrap BMDs are infinite — target not reached in any sample")
    
    return(list(
      AICs   = AICs,
      weights= weights,
      BMD    = BMD,
      BMDCI  = c(NA_real_, NA_real_)
    ))
  }
  
  # Case 3: upper CI infinite -> valid but informative
  if (is.infinite(BMDCI[2])) {
    message("Upper CI is infinite — BMD exceeds tested dose range")
  }
  
  # Build long data.frame for plotting
  df_boot_long = melt(boot_preds)
  # melt gives: Var1 = bootstrap index, Var2 = column index (fitdose idx), value = prediction
  colnames(df_boot_long) = c("sample", "fitdose_idx", "response")
  df_boot_long$fitdose_idx = as.integer(df_boot_long$fitdose_idx)
  df_boot_long$dose = fitdose[df_boot_long$fitdose_idx]
  df_boot_long$sample = factor(df_boot_long$sample)
  df_boot_longP = df_boot_long[df_boot_long$dose >= offset_start, ] #Store a version for plotting that removes estimates from the concentrations below offset_start
  
  # data.frame for original weighted prediction (red line)
  df_orig = data.frame(dose = fitdose, response = wp)
  df_orig = df_orig[df_orig$dose >= offset_start, ]
  
  # bootstrap predictions at 0
  boot_zero = boot_preds[, fitdose == 0, drop = TRUE]
  # bootstrap predictions at offset_start
  boot_start = boot_preds[, fitdose == offset_start, drop = TRUE]
  
  # segment data frame for bootstrap (blue)
  df_seg_boot = data.frame(
    sample  = factor(seq_along(boot_zero)),
    x       = offset_zero,
    xend    = offset_start,
    y       = boot_zero,
    yend    = boot_start
  )
  
  predresp_zero = wp[fitdose == 0]
  predresp_start = wp[fitdose == offset_start]
  
  # segment data frame for original (red)
  df_seg_orig = data.frame(
    x = offset_zero,
    y = predresp_zero,
    xend = offset_start,
    yend = predresp_start
  )
  
  df_data = data.frame(dose = dose, response = response, SEM = SEM)
  df_data$dose_plot = ifelse(df_data$dose == 0, offset_zero, df_data$dose)
  maxplotres = max(c(max(response) * 1.5, max(response) + 2 * SEM[response == max(response)]))
  minplotres = min(c(min(response) * .75, min(response) - 2 * SEM[response == min(response)]))
  
  breaks_vec = c(offset_zero,
                 10^seq(floor(log10(min_nonzero)), ceiling(log10(max(dose))),
                        by = 1))
  
  fulltitle = paste0(pboottitle, " (n=", nb_valid, ")")
  
  
  
  # Plot: bootstrap lines (blue, semi-transparent), original weighted prediction (red)
  pboot = ggplot() +
    # 1. BOOTSTRAP dashed segments
    geom_segment(
      data = df_seg_boot,
      aes(x = x, y = y, xend = xend, yend = yend, group = sample),
      color = "blue",
      linetype = "dashed",
      alpha = 0.10
    ) +
    # 2. ORIGINAL dashed (red) segment
    geom_segment(
      data = df_seg_orig,
      aes(x = x, y = y, xend = xend, yend = yend),
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    
    # 3. Bootstrap curves (starting at offset_start)
    geom_line(
      data = df_boot_longP,
      aes(x = dose, y = response, group = sample),
      color = "blue",
      alpha = 0.08,
      linewidth = 0.6
    ) +
    
    # 4. Original weighted curve (red)
    geom_line(
      data = df_orig,
      aes(x = dose, y = response),
      color = "red",
      linewidth = 1.2
    ) +
    # 5. Error bars
    geom_errorbar(
      data = df_data,
      aes(x = dose_plot, ymin = response - SEM, ymax = response + SEM),
      color = "black",
      width = .1
    ) +
    # 6. Points
    geom_point(
      data = df_data,
      aes(x = dose_plot, y = response),
      size = 2
    ) +
    
    scale_x_continuous(
      trans = "log10",
      breaks = breaks_vec,
      labels = c("0", formatC(breaks_vec[-1], format = "fg"))
    ) +
    
    #if(df_data$response < 0)
    coord_cartesian(
      ylim = c(minplotres, maxplotres)
    ) +
    
    labs(x = xlab, y = ylab, title = fulltitle) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(pboot)
  ggsave(paste0(savepath, 'bootplots/', basesave, '_bootplot.png'), width = 5, height = 5, dpi = 500)
  
  return(list(AICs = AICs, boot_dfs = boot_dfs_store, boot_weights = boot_weights_store, models = models, boot_best_models = boot_best_models_store, boot_fits = boot_fits_store, BMD = BMD, BMDCI = BMDCI, boot_long = df_boot_long, boot_BMD = boot_BMD, boot_global = boot_global_store))
}

derive_group_vars = function(df, assay_name) {
  
  # -------------------------
  # Preprocessing
  # -------------------------
  #Exposure time filter
  if ("Time of exposure [h]" %in% colnames(df)) {
    df = df %>%
      filter(`Time of exposure [h]` >= 1)
  }
  
  #Remove unwanted loci
  if ("Gene locus" %in% colnames(df)) {
    df = df %>%
      filter(`Gene locus` %in% c("HPRT", "TK"))
  }
  
  #Remove S9 == yes
  if ("S9" %in% colnames(df)) {
    df = df %>%
      filter(S9 != "yes")
  }
  
  #Remove specific DNA adduct types
  if ("distinct" %in% colnames(df)) {
    df = df %>%
      filter(!distinct %in% c("N-7G", "06G", "N-3A"))
  }
  
  # -------------------------
  # Exposure
  # -------------------------
  df = df %>%
    mutate(
      exposure_group = case_when(
        `Time of exposure [h]` <= 10 ~ "S",
        `Time of exposure [h]` > 10  ~ "L",
        TRUE ~ NA_character_
      )
    )
  
  
  # -------------------------
  # Expression
  # -------------------------
  if ("Expression period[days]" %in% colnames(df)) {
    
    df = df %>%
      mutate(
        expression_group = case_when(
          `Expression period[days]` <= (1/24) ~ "I",
          `Expression period[days]` > (1/24) & `Expression period[days]` < 6 ~ "E",
          `Expression period[days]` >= 6 ~ "L",
          TRUE ~ "E"
        )
      )
    
  } else if ("expression_period" %in% colnames(df)) {
    
    df = df %>%
      mutate(
        expression_group = case_when(
          expression_period <= (1/24) ~ "I",
          expression_period > (1/24) & expression_period < 6 ~ "E",
          expression_period >= 6 ~ "L",
          TRUE ~ "E"
        )
      )
    
  } else {
    df$expression_group = "I"
  }
  
  
  # -------------------------
  # Scoring
  # -------------------------
  if ("Scoring technique" %in% colnames(df)) {
    
    df = df %>%
      mutate(
        scoring_group = case_when(
          grepl("flow", `Scoring technique`, ignore.case = TRUE)  ~ "F",
          grepl("micro", `Scoring technique`, ignore.case = TRUE) ~ "M",
          TRUE ~ "O"
        )
      )
    
  } else if (grepl("comet", assay_name, ignore.case = T)) {
    df$scoring_group = "M"
  } else {
    df$scoring_group = "O"
  }
  
  
  # -------------------------
  # Extra variable
  # -------------------------
  df$extra_group = NA_character_
  
  if (assay_name == "MN" && "Assay type" %in% colnames(df)) {
    
    df = df %>%
      mutate(
        extra_group = case_when(
          grepl("CBMN", `Assay type`, ignore.case = TRUE) ~ "B",
          grepl("conventional", `Assay type`, ignore.case = TRUE) ~ "C",
          TRUE ~ NA_character_
        )
      )
  }
  
  if (assay_name == "muta" && "Gene locus" %in% colnames(df)) {
    
    df = df %>%
      mutate(
        extra_group = case_when(
          `Gene locus` == "TK"   ~ "T",
          `Gene locus` == "HPRT" ~ "H",
          TRUE ~ NA_character_
        )
      )
  }
  
  if (assay_name == "ECDS" && "Target region" %in% colnames(df)) {
    
    df = df %>%
      mutate(
        extra_group = case_when(
          grepl("Whole genome", `Target region`, ignore.case = TRUE) ~ "W",
          grepl("20 target", `Target region`, ignore.case = TRUE) ~ "T",
          TRUE ~ NA_character_
        )
      )
  }
  
  df
}

# Helper to convert one results list to a dataframe
convert_results = function(results_list, compound_name) {
  map_dfr(results_list, function(x) {
    tibble(
      compound = compound_name,
      assay = x$assay,
      subset = x$subset,
      BMD = x$BMD,
      BMD_low = x$BMD_low,
      BMD_high = x$BMD_high
    )
  })
}
#Documentation on the functions used for data wrangling as described in 
#'Advancing AOPs for genotoxicity through quantifying KERs linking DNA 
#'alkylation to increase in mutations and chromosome aberrations' by Winkelman
#et al., 2026.


#generalwrangle:
#  Goal:                  Perform some steps for general data wrangling
#Arguments:                 
#- df:                    input dataframe  
#- conc_col:              name of the column specifying the concentration
#- outcome_col:           name of the column specifying the outcome measure
#- rs_col:                name of the column specifying the relative survival
#- ctthresh:              threshold of relative survival at which responses are still deemed relevant (in %)
#Note:                    df requires a last column called 'Comments', a column titled: 
#                         `Relative survival (%)` and a column called 'Reference'
#                         structured like: "<Author>., <Year>, doi:<DOI>"


#convert_outcome_units:
#  Goal:                  Standardize micronucleus (MN) outcome units to percentages.
#Arguments:
#- df:                    input dataframe
#- unit_col:              name of the column specifying the outcome unit
#- outcome_col:           name of the column specifying the treatment outcome
#Note:                    Recognized MN units are mapped to either %MN (mononucleated
#                         cells) or %BNMN (binucleated cells). Counts per 1,000 or
#                         2,000 cells are converted to percentages. Unrecognized
#                         unit labels cause an error.


#convert_muta_to_perm:
#  Goal:                  Convert mutation-frequency estimates to mutations per million cells.
#Arguments:
#- df:                    input dataframe
#- outcome_col:           name of the column specifying the treatment outcome
#- lower_col:             name of the column specifying the lower variability estimate
#- upper_col:             name of the column specifying the upper variability estimate
#- unit_col:              name of the column specifying the mutation-frequency unit
#Note:                    The exponent in the unit label is used to rescale the treatment
#                         outcome and both variability estimates to a 10^-6 (per million)
#                         scale. The unit column is then standardized to
#                         "Mutant frequency x 10^6 (per million cells)".


#ObsToMeans:
#  Goal:                  Convert individual observations into means and standard deviations.
#Arguments:
#- TK6df:                 input dataframe containing a 'replicates' column
#- varOI:                 name of the variable used to define groups (eg assay or study condition)
#- ct:                    whether to calculate and retain mean relative survival (%) (default: FALSE)
#- varextra:              optional name of an additional numeric variable to average (default: NA)
#Note:                    Rows with replicates == 1 are treated as individual observations.
#                         They are summarized by varOI and concentration; all other rows
#                         are retained as already summarized data. The function calculates
#                         the mean, SD and number of observations for each group and adds
#                         lower/upper variability estimates. Metadata must be consistent
#                         within each varOI/concentration group. If ct = TRUE, the input
#                         must contain 'Relative survival (%)'. If varextra is supplied,
#                         that column is averaged and added back as 'Amount of DNA analyzed'.


#remove_few_obs_studies:
#  Goal:                  Remove studies with two observations or fewer.
#Arguments:
#- df:                    input dataframe containing a 'Study index' column
#Note:                    Retains only rows belonging to studies with more than two
#                         observations in the input dataframe.


#convert_varest:
#  Goal:                  Standardize variability estimates by calculating SD and SEM.
#Arguments:
#- df:                    input dataframe containing variability estimates, their unit,
#                         the treatment outcome and the replicate count
#Note:                    Recognizes SD, SEM, IQR, 1.5 IQR and 95% CI estimates. Variability
#                         estimates are first converted to SD and then used to calculate
#                         SEM. The resulting values are added in columns named 'SD' and
#                         'SEM'. Conversion assumes the supplied bounds represent a
#                         symmetric variability estimate around the treatment outcome.


#convert_concs:
#  Goal:                  Convert concentration values to µM.
#Arguments:
#- df:                    input dataframe containing a concentration column and its unit
#- molar_weight:          molecular weight of the tested compound (g/mol)
#Note:                    Detects a concentration column named 'Concentration EMS',
#                         'Concentration MMS' or 'Concentration ENU'. Supports µM, µg/mL,
#                         mg/mL and mM, including corresponding log10(conc) unit labels.
#                         Log10-transformed concentrations are first back-transformed.
#                         All concentration values and unit labels are then standardized
#                         to the linear µM scale.


#correctclass:
#  Goal:                  Convert selected dataframe columns to numeric or factor classes.
#Arguments:
#- df:                    input dataframe
#- extranumcols:          optional additional column names to convert to numeric (default: NULL)
#- extrafaccols:          optional additional column names to convert to factors (default: NULL)
#Note:                    Converts recognized standard columns (eg Year, concentration,
#                         outcome and variability estimates) to numeric and study/reference
#                         metadata columns to factors. Columns not present in df are ignored.


#convert_CTassaytype:
#  Goal:                  Standardize labels describing the cytotoxicity assay type.
#Arguments:
#- df:                    input dataframe
#Note:                    Applies predefined mappings to columns whose names start with
#                         'Type of cytotoxicity assay' (case-insensitive). Recognized
#                         assay descriptions are replaced by standardized labels and
#                         converted to factors. Unmatched labels are left unchanged.


#convert_scoringTechnique:
#  Goal:                  Standardize labels describing the scoring or measurement technique.
#Arguments:
#- df:                    input dataframe
#Note:                    Applies predefined mappings to columns whose names start with
#                         'technique for scoring' or 'type of measurement'
#                         (case-insensitive). Recognized descriptions are replaced by
#                         standardized labels (eg CA4, Komet6, Metafer, flowCyto or manual)
#                         and converted to factors. Unmatched labels are left unchanged.


#convert_locus_naming:
#  Goal:                  Standardize gene-locus names.
#Arguments:
#- df:                    input dataframe
#Note:                    Applies predefined mappings to columns whose names start with
#                         'Gene locus' (case-insensitive). Recognized locus names are
#                         standardized (eg TK, HPRT, NAKATPase, TUB and RNAsyn) and
#                         converted to factors. Unmatched labels are left unchanged.


#convert_to_FC:
#  Goal:                  Calculate fold changes relative to the concurrent control.
#Arguments:
#- df:                    input dataframe containing study, concentration, outcome,
#                         SD and replicate information
#- conc_col:              name of the column specifying the concentration
#- outcome_col:           name of the column specifying the treatment outcome
#Note:                    Calculates fold change (FC) within each study using the
#                         zero-concentration response as the control. Propagates the
#                         treatment and control SDs to calculate FC_SD and derives
#                         FC_SEM using the replicate count. Also adds FC_lower and
#                         FC_upper as FC ± FC_SEM. Each study therefore requires a
#                         zero-concentration control observation.


#doseres_plot:
#  Goal:                  Plot dose-response data for individual studies.
#Arguments:
#- TK6_df:                input dataframe containing concentration, response, SEM,
#                         'Study index', 'Author' and 'Year' columns
#- title:                 plot title suffix
#- assay:                 assay name used in the plot title and output filename
#- conc_col:              name of the concentration column (default: "Concentration EMS")
#- response_col:          name of the response column (default: "FC")
#- sem_col:               name of the SEM column (default: "FC_SEM")
#- width:                 width of the saved plot in inches (default: 7)
#- height:                height of the saved plot in inches (default: 5)
#- save:                  whether to save the plot as a PNG (default: FALSE)
#- path:                  directory in which to save the plot (default: NULL)
#Note:                    Displays a log-log dose-response plot with log-symmetric error
#                         bars. Zero-concentration observations are displayed at an
#                         offset below the lowest non-zero concentration. If a 'distinct'
#                         column is present, it is used to distinguish studies by the
#                         same author in the legend; otherwise, Author is used.


#doseres_contrast:
#  Goal:                  Plot dose-response data while contrasting groups defined by dataframe variables.
#Arguments:
#- TK6_df:                input dataframe containing concentration, response, SEM and 'Study index'
#- title:                 plot title suffix
#- assay:                 assay name used in the plot title and output filename
#- conc_col:              name of the concentration column (default: "Concentration EMS")
#- response_col:          name of the response column (default: "FC")
#- sem_col:               name of the SEM column (default: "FC_SEM")
#- color_col:             name of the variable mapped to point/line color (required)
#- shape_col:             optional name of the variable mapped to point shape (default: NULL)
#- linetype_col:          optional name of the variable mapped to line type (default: NULL)
#- solid_level:           optional level of linetype_col to display as a solid line (default: NULL)
#- factor_colors:         optional named vector specifying colors for color_col levels (default: NULL)
#- width:                 width of the saved plot in inches (default: 7)
#- height:                height of the saved plot in inches (default: 5)
#- save:                  whether to save the plot as a PNG (default: FALSE)
#- path:                  directory in which to save the plot (default: NULL)
#Note:                    Displays a log-log dose-response plot with log-symmetric error
#                         bars. Zero-concentration observations are displayed at an
#                         offset below the lowest non-zero concentration. Numeric color
#                         variables are treated as ordered factors; categorical variables
#                         are treated as factors. Optional shape and line-type mappings
#                         allow additional grouping variables to be displayed.


#plot_lmer_dose_response:
#  Goal:                  Plot observed data and fixed-effect predictions from a linear mixed-effects model.
#Arguments:
#- model:                 fitted model accepted by predict() and model.frame()
#- data:                  dataframe used to construct the prediction grid and display observations
#- conc_var:              name of the concentration variable
#- response_var:          name of the response variable (must be on the log scale)
#- group_var:             optional grouping variable for separate prediction curves (default: NULL)
#- sem_var:               optional name of the SEM variable (must be on the log scale; default: NULL)
#- factor_reference:      named list specifying reference levels for categorical predictors (default: list())
#- grid_length:           number of non-zero concentration values in the prediction grid (default: 200)
#- offset_fraction:       fraction of the minimum non-zero concentration used to display zero doses (default: 0.25)
#- save_plot:             whether to save the plot as a PNG (default: FALSE)
#- save_dir:              directory in which to save the plot (default: ".")
#- width:                 width of the saved plot in inches (default: 6)
#- height:                height of the saved plot in inches (default: 5)
#- dpi:                   resolution of the saved plot (default: 300)
#Note:                    Numeric predictors are fixed at their mean. Categorical predictors
#                         are fixed at the supplied factor_reference level, or otherwise
#                         their most frequent level. If group_var is supplied, predictions
#                         are generated separately for each of its levels. Predictions
#                         exclude random effects (re.form = NA). Returns a list containing
#                         the ggplot object and the prediction dataframe.


#prepare_varcomp:
#  Goal:                  Prepare response and variability data for variance decomposition on the log scale.
#Arguments:
#- df:                    input dataframe
#- moddf:                 optional dataframe supplying additional columns (default: NULL; if NULL, df is used)
#- out_col:               name of the response column to transform (default: "FC")
#- sem_col:               name of the SEM column (default: "FC_SEM")
#- sd_col:                name of the SD column (default: "FC_SD")
#- rep_col:               name of the replicate-count column (default: "replicates")
#Note:                    Adds columns from moddf that are not already present in df.
#                         Applies a log-scale bias correction to the response and calculates
#                         the log-scale SEM from the coefficient of variation and replicate
#                         count. For the default FC/FC_SEM combination, the output columns
#                         are named 'FClog' and 'FC_SEMlog'; otherwise, names are formed
#                         from the supplied column names.


#doseres_contrast_multi:
#  Goal:                  Combine multiple datasets and plot their dose-response data in separate panels.
#Arguments:
#- df_list:               list of input dataframes
#- df_names:              optional names assigned to the datasets (default: NULL; uses list names or Dataset_1, Dataset_2, ...)
#- title:                 plot title suffix
#- assay:                 assay name used in the plot title and output filename
#- response_col:          name of the response column (default: "FC")
#- sem_col:               name of the SEM column (default: "FC_SEM")
#- color_col:             name of the variable mapped to point/line color (required)
#- shape_col:             optional name of the variable mapped to point shape (default: NULL)
#- linetype_col:          optional name of the variable mapped to line type (default: NULL)
#- solid_level:           optional level of linetype_col to display as a solid line (default: NULL)
#- width:                 width of the saved plot in inches (default: 10)
#- height:                height of the saved plot in inches (default: 5)
#- save:                  whether to save the plot as a PNG (default: TRUE)
#Note:                    Each dataframe must contain exactly one concentration column
#                         matching 'Concentration [one to three letters]'. The datasets
#                         are combined, concentrations are displayed on a common log
#                         scale, and panels are faceted by dataset (EMS, MMS and ENU).
#                         Error bars are calculated from the response and SEM columns.
#                         The plot is saved in the working directory when save = TRUE.




generalwrangle = function(df,
                          conc_col = "Concentration EMS",
                          outcome_col = "Outcome for the treatment",
                          rs_col = "Relative survival (%)",
                          ctthresh = 50){
  
  # 1 remove empty columns at the end of the sheet
  df = df[,1:which(colnames(df) == 'Comments')]
  
  # 2 Reformat column names
  colnames(df) = gsub('[\r\n]', '', colnames(df))
  
  # 3 Reformat the Reference column
  df$Reference = gsub('[\r\n]', '', df$Reference)
  
  # 4 remove rows with missing study index
  df = df[!is.na(df$`Study index`),]
  
  # 5 Ensure that NAs are read in correctly
  df[df == 'NA'] = NA
  
  # 6 Properly read in the reference
  df = separate_wider_delim(df, cols = Reference, delim = '.,',
                            names = c('Author', 'Reference'))
  df = separate_wider_delim(df, cols = Reference, delim = ',',
                            names = c('Year', 'DOI'))
  
  df$DOI = trimws(df$DOI)
  df$DOI = ifelse(
    grepl("^PMID", df$DOI, ignore.case = TRUE),
    df$DOI,
    gsub(".*:", "", df$DOI)
  )
  
  # 7 Remove all columns containing NA values only
  df = df[, colSums(is.na(df)) < nrow(df)]
  
  # 8 Remove observations for which no outcome was measured
  df = df[!is.na(df[[outcome_col]]), ]
  
  # 9 Mark overly cytotoxic observations
  if (rs_col %in% names(df)) {
    df = df %>%
      group_by(`Study index`, .data[[conc_col]]) %>%
      mutate(
        ct = mean(as.numeric(.data[[rs_col]]), na.rm = TRUE) >= ctthresh
      ) %>%
      ungroup()
  } else {
    df = df %>%
      mutate(ct = NA)
  }
  
  
  # 10 Remove redundant control columns
  df = df[, !(names(df) %in% 
                c("Outcome for the control",
                  "Number of cells scored control"))]
  
  # 11 Replicates correction
  df$replicates[
    df$`Type of outcome: (mean / median / individual observation)` == 
      'Individual observation'
  ] = 1
  
  # 12 convert to dataframe
  df = as.data.frame(df)
  
  return(df)
}

convert_outcome_units = function(
    df,
    unit_col = "Outcome unit (eg percentage of binucleated cells containing one or more micronuclei)",
    outcome_col = "Outcome for the treatment"
){
  
  assay_map = list(
    MNperc = c("% MN (no of micronuclei / no of nuclei * 100)",
               "% of MN",
               "Percent of cells with micronuclei",
               "Percentage of micronucleated cells",
               "% MN",
               "% MN Frequency",
               "% MN Frequency of mononucleated cells",
               "%MN",
               "%MN mononucleates",
               "Micronucleated cells (%)",
               "Fold increase of %MN"),
    
    MNBNperc = c("(% binucleated cells with MN",
                 "% binucleated cells with MN",
                 "% MNBNC",
                 "Bi-nucleated cell MN frequency (%)",
                 "% BnMN freqeuncy"),
    
    MNBNproK = c("MN-cells per 1000 BN-cells",
                 "MN/1000 BN"),
    
    MNpro2K = c("MN/2000 cells"),
    
    MNBNpro2K = c("Number of BN cells with micronuclei per 2000 binucleated cells")
  )
  
  # Build reverse lookup
  map_df = stack(assay_map)
  colnames(map_df) = c("raw_label", "map_class")
  
  # Match each row's unit to mapping
  match_idx = match(df[[unit_col]], map_df$raw_label)
  
  # Fail if unknown labels appear
  if(any(is.na(match_idx))){
    unknown = unique(df[[unit_col]][is.na(match_idx)])
    stop("Unmapped outcome unit(s) found:\n", paste(unknown, collapse = "\n"))
  }
  
  # Assign mapped class
  mapped_class = map_df$map_class[match_idx]
  
  # --- Convert to percent ---
  
  # per 1000 BN cells → percent
  df[[outcome_col]][mapped_class == "MNBNproK"] =
    df[[outcome_col]][mapped_class == "MNBNproK"] / 10
  
  # per 2000 mononucleated cells → percent
  df[[outcome_col]][mapped_class == "MNpro2K"] =
    df[[outcome_col]][mapped_class == "MNpro2K"] / 20
  
  # per 2000 binucleated cells → percent
  df[[outcome_col]][mapped_class == "MNBNpro2K"] =
    df[[outcome_col]][mapped_class == "MNBNpro2K"] / 20
  
  
  # --- Collapse to final classes ---
  df[[unit_col]] = dplyr::case_when(
    mapped_class %in% c("MNperc", "MNpro2K") ~ "%MN",
    mapped_class %in% c("MNBNperc", "MNBNproK", "MNBNpro2K") ~ "%BNMN"
  )
  
  df[[unit_col]] = as.factor(df[[unit_col]])
  
  return(df)
}


convert_muta_to_perm = function(df,
                                outcome_col = "Outcome for the treatment",
                                lower_col = "Outcome variability lower estimate",
                                upper_col = "Outcome variability upper estimate",
                                unit_col = "Outcome unit (eg Mutation frequency in percentage)") {
  
  extract_exponent = function(unit) {
    
    exp = stringr::str_extract(unit, "\\^\\s*-?\\d+")
    exp = gsub("\\^", "", exp)
    exp[is.na(exp)] = 0
    
    exp = abs(as.numeric(exp))
    
    -exp
  }
  
  exponent = extract_exponent(df[[unit_col]])
  
  # Convert exponent relative to 10^-6 scale
  rel_exp = exponent + 6
  
  scale_factor = 10^(rel_exp)
  
  df[[outcome_col]] = df[[outcome_col]] * scale_factor
  df[[lower_col]]   = df[[lower_col]]   * scale_factor
  df[[upper_col]]   = df[[upper_col]]   * scale_factor
  
  df[[unit_col]] = "Mutant frequency x 10^6 (per million cells)"
  
  return(df)
}

#This function converts a dataframe with individual observations to one with means and SDs
ObsToMeans = function(TK6df, varOI, ct = FALSE, varextra = NA){
  
  # Split raw vs already summarized data
  if(!"replicates" %in% colnames(TK6df)){
    stop("Input dataframe must contain a 'replicates' column.")
  }
  
  #Convert factors to character to aid eventual merging of dataframes
  TK6df = TK6df %>%
    mutate(across(where(is.factor), as.character))
  
  raw_df  = TK6df[!is.na(TK6df$replicates) & TK6df$replicates == 1, ]
  keep_df = TK6df[TK6df$replicates != 1 | is.na(TK6df$replicates), ]
  
  levelsOI = unique(raw_df[[varOI]])
  
  # If no raw data, return unchanged
  if(nrow(raw_df) == 0){
    return(TK6df)
  }
  
  uniq = levelmean = levelsd = ctmean = meanextra = check = replicates = partial_dfs =
    vector(mode = 'list', length = length(levelsOI))
  
  # Optional cytotoxicity column
  if(ct == FALSE){
    ctmean = NULL
    colct = NULL
  } else {
    colct = which(colnames(raw_df) == 'Relative survival (%)')
  }
  
  # Optional extra variable
  if(is.na(varextra)){
    colextra = NULL
  } else {
    colextra = which(colnames(raw_df) == varextra)
  }
  
  # Required columns
  colconc = which(colnames(raw_df) %in% c('Concentration EMS', 'Concentration MMS', 'Concentration ENU'))
  colOI   = which(colnames(raw_df) == varOI)
  colout  = which(colnames(raw_df) == 'Outcome for the treatment')
  
  # Loop through levels of interest
  for(i in seq_along(levelsOI)){
    
    uniq[[i]] = unique(raw_df[raw_df[, colOI] == levelsOI[i], colconc])
    levelmean[[i]] = levelsd[[i]] = ctmean[[i]] = meanextra[[i]] = replicates[[i]] = check[[i]] =
      numeric(length(uniq[[i]]))
    
    # Compute means, SDs, replicate counts
    for(j in seq_along(uniq[[i]])){
      subset_data = raw_df[
        raw_df[, colOI] == levelsOI[i] & raw_df[, colconc] == uniq[[i]][j], 
      ]
      
      levelmean[[i]][j]  = mean(subset_data[, colout], na.rm = TRUE)
      levelsd[[i]][j]    = sd(subset_data[, colout], na.rm = TRUE)
      replicates[[i]][j] = nrow(subset_data)
      
      if(ct == TRUE){
        ctmean[[i]][j] = mean(subset_data[, colct], na.rm = TRUE)
      }
      if(!is.na(varextra)){
        meanextra[[i]][j] = mean(subset_data[, colextra], na.rm = TRUE)
      }
    }
    
    # Variation check (only metadata columns)
    for(j in seq_along(uniq[[i]])){
      check[[i]][j] = nrow(
        distinct(
          raw_df[
            raw_df[, colOI] == levelsOI[i] & raw_df[, colconc] == uniq[[i]][j], 
            -c(colct, colextra, colout)
          ]
        )
      )
    }
    
    if(sum(check[[i]]) != length(check[[i]])){
      stop('The dataframe varies on more variables than the ones provided')
    }
    
    # Build summarized dataframe
    drop_cols = c(colct, colextra, colout)
    drop_cols = drop_cols[!is.na(drop_cols)]
    
    partial_dfs[[i]] = as.data.frame(
      matrix(
        nrow = length(uniq[[i]]), 
        ncol = ncol(raw_df) - length(drop_cols)
      )
    )
    
    colnames(partial_dfs[[i]]) = colnames(raw_df)[-drop_cols]
    
    # Fill metadata rows
    for(j in seq_len(nrow(partial_dfs[[i]]))){
      partial_dfs[[i]][j, ] = distinct(
        raw_df[
          raw_df[, colOI] == levelsOI[i] & raw_df[, colconc] == uniq[[i]][j], 
          -drop_cols
        ]
      )
    }
    
    # Add back computed columns
    if(ct == TRUE){
      partial_dfs[[i]] = add_column(
        partial_dfs[[i]], 
        `Relative survival (%)` = ctmean[[i]], 
        .after = colct - 1
      )
    }
    
    if(!is.na(varextra)){
      partial_dfs[[i]] = add_column(
        partial_dfs[[i]], 
        `Amount of DNA analyzed` = meanextra[[i]], 
        .after = colextra - 1
      )
    }
    
    partial_dfs[[i]] = add_column(
      partial_dfs[[i]], 
      `Outcome for the treatment` = levelmean[[i]], 
      .after = colout - 1
    )
    
    partial_dfs[[i]]$`Outcome variability lower estimate` = levelmean[[i]] - levelsd[[i]]
    partial_dfs[[i]]$`Outcome variability upper estimate` = levelmean[[i]] + levelsd[[i]]
    
    partial_dfs[[i]]$`Outcome variability unit (eg std deviation)` = "std dev"
    partial_dfs[[i]]$`Type of outcome: (mean / median / individual observation)` = "mean"
    partial_dfs[[i]]$replicates = replicates[[i]]
  }
  
  summarized_raw = bind_rows(partial_dfs)
  
  # Combine summarized raw + untouched summarized rows
  final_out = bind_rows(keep_df, summarized_raw)
  #Convert data back to desired classes
  final_out = correctclass(final_out)
  final_out = final_out %>%
    dplyr::mutate(across(where(is.numeric), ~ replace(., is.nan(.), NA)))
  
  return(final_out)
}


#This function removes studies with 2 observations or less in the dataset
remove_few_obs_studies = function(df){
  counts = table(df$`Study index`)
  df[counts[df$`Study index`] > 2, ]
}

#Convert estimates of variability to the same unit
convert_varest = function(df) {
  
  SDnames      = c("Std deviation", "std deviation", "SD", "std dev", "Std dev")
  SEMnames     = c("Std error", "std error", "standard error", "SEM") 
  IQRnames     = c("IQR")
  IQR1.5names  = c("1.5 IQR", "1.5IQR", "1.5 IQR (interquartile range)")
  CI95names    = c("95% CI", "95%CI", "95% CI calculated from Poisson distribution")
  
  unit_col = "Outcome variability unit (eg std deviation)"
  
  SD_index     = which(df[[unit_col]] %in% SDnames)
  SEM_index    = which(df[[unit_col]] %in% SEMnames)
  IQR_index    = which(df[[unit_col]] %in% IQRnames)
  IQR1.5_index = which(df[[unit_col]] %in% IQR1.5names)
  CI95_index   = which(df[[unit_col]] %in% CI95names)
  
  ## symmetric variability estimate (half-width)
  varest = (df$`Outcome variability upper estimate`
            - df$`Outcome for the treatment`
            + df$`Outcome for the treatment`
            - df$`Outcome variability lower estimate`) / 2
  
  ## convert EVERYTHING to SD first
  SD = varest
  
  SD[SEM_index]    = SD[SEM_index] * sqrt(df$replicates[SEM_index])
  SD[IQR_index]    = SD[IQR_index] / 1.35
  SD[IQR1.5_index] = SD[IQR1.5_index] / (1.35 * 1.5)
  SD[CI95_index]   = (varest[CI95_index] / 1.96) * sqrt(df$replicates[CI95_index])
  
  ## derive SEM from SD
  SEM = SD / sqrt(df$replicates)
  
  ## append columns
  df$SD  = SD
  df$SEM = SEM
  
  return(df)
}

convert_concs = function(df, molar_weight){
  
  # supported unit variants
  um_units   = c('µM', "μM", 'uM')
  ugml_units = c('µg/mL', 'µg/ml', 'ug/mL', 'ug/ml', 'μg/ml', 'μg/mL')
  mgml_units = c('mg/mL', 'mg/ml')
  mm_units   = c('mM')
  
  # log10 variants
  log_um_units   = paste0("log10(conc) ", um_units)
  log_ugml_units = paste0("log10(conc) ", ugml_units)
  log_mgml_units = paste0("log10(conc) ", mgml_units)
  log_mm_units   = paste0("log10(conc) ", mm_units)
  
  allunits = c(
    um_units, ugml_units, mgml_units, mm_units,
    log_um_units, log_ugml_units, log_mgml_units, log_mm_units
  )
  
  # detect concentration column dynamically
  conc_col = grep("^Concentration (EMS|MMS|ENU)$", colnames(df), value = TRUE)
  
  if(length(conc_col) != 1){
    stop("Could not find a concentration column named 'Concentration EMS', 'MMS', or 'ENU'")
  }
  
  unit_col = "Concentration unit (eg μg/ml or μM)"
  
  if(!(unit_col %in% colnames(df))){
    stop("Missing unit column: 'Concentration unit (eg μg/ml or μM)'")
  }
  
  # check for invalid units
  offending_units = unique(df[[unit_col]])[!(unique(df[[unit_col]]) %in% allunits)]
  
  if(length(offending_units) > 0){
    stop(paste("Unidentified concentration unit(s):", paste(offending_units, collapse = ", ")))
  }
  
  # convert log10 concentrations to linear scale 
  log_idx = grepl("^log10\\(conc\\)", df[[unit_col]])
  df[[conc_col]][log_idx] = 10^(df[[conc_col]][log_idx])
  df[[unit_col]] = gsub("^log10\\(conc\\) ", "", df[[unit_col]])
  
  # µg/mL -> µM
  ugml_idx = df[[unit_col]] %in% ugml_units
  df[[conc_col]][ugml_idx] = (df[[conc_col]][ugml_idx] / molar_weight) * 1000
  
  # mg/mL -> µM
  mgml_idx = df[[unit_col]] %in% mgml_units
  df[[conc_col]][mgml_idx] = (df[[conc_col]][mgml_idx] / molar_weight) * 1e6
  
  # mM -> µM
  mm_idx = df[[unit_col]] %in% mm_units
  df[[conc_col]][mm_idx] = df[[conc_col]][mm_idx] * 1000
  
  # standardize all units to µM
  df[[unit_col]] = "µM"
  
  return(df)
}


#Convert data classes to desirable types
correctclass = function(df, extranumcols = NULL, extrafaccols = NULL){
  
  numcols = c('Year',
              'Time of exposure [h]',
              'Concentration EMS',
              'Concentration MMS',
              'Concentration ENU',
              'Time of harvest [h]',
              'Percentage DMSO in medium (%)',
              'Relative survival (%)',
              'Time of cytotoxicity assessment [h]',
              'Outcome for the treatment',
              'Outcome for the control',
              'Outcome variability lower estimate',
              'Outcome variability upper estimate',
              'Number of cells scored treatment',
              'Number of cells scored control',
              'cell count treatment',
              'cell count control',
              'replicates',
              'Amount of DNA analyzed',
              'Expression period[days]')
  
  faccols = c('Study index',
              'Author',
              'DOI',
              'Evaluator')
  
  # Add extra requested columns
  numcols = union(numcols, extranumcols)
  faccols = union(faccols, extrafaccols)
  
  # Keep only columns that actually exist in df
  numcols = intersect(numcols, colnames(df))
  faccols = intersect(faccols, colnames(df))
  
  # Convert classes
  df = df %>% mutate(across(all_of(numcols), as.numeric))
  df = df %>% mutate(across(all_of(faccols), as.factor))
  
  return(df)
}

convert_CTassaytype = function(df){
  
  assay_map = list(
    trypan = c("trypan blue dye-exclusion assay", "Trypan Blue exclusion", "Trypan blue exclusion", "trypan blue dye-exclusion", "Trypan Blue", "trypan"),
    RSG = c("RSG (relative suspension growth)", "Relative Growth"),
    celltiter = c("CellTiter-Glo Luminescent Cell Viability Assay"),
    MTS = c("MTS", "MTS assay"),
    RCC = c("Relative Cell Count", "RCC"),
    RICC = c("RICC", "Relative growth compared to the 0 conc group", "RCG",
             "RCG (Sysmex microcell counter F-300)",
             "relative increase in cell counts (RICC)",
             "RICC (% of control)"),
    RNC = c("RNC"),
    BTC = c("beads:2:cell ratio", "Ratio of viable cells to counting beads (Microflow)", "nuclei-to-bead ratio of all treated samples by that of the corresponding vehicle control"),
    CE = c("Cloning efficiency at day 0",
           "Comparing colonoy-forming ability of treated cells to control",
           "Comparing colony-forming ability of treated cells to control", "Cloning efficiency", "Plating efficiency"),
    RPD = c("RPD", 'relative population doubling', '%RPD', 'Relative Population Doubling', 'Relative population doubling'),
    FDA = c("FDA viability staining", "fluorescein diacetate/GelRed Nucleic Acid Gel Stain"),
    cytostasis = c("Cell growth inhibition (%)", "percentage of healthy proliferating cell", "Cytostasis %", "Calculated from binucleate frequencies")
  )
  
  # Find relevant columns (case-insensitive)
  assay_cols = names(df)[
    grepl("^Type of cytotoxicity assay", names(df), ignore.case = TRUE)
  ]
  
  # Apply mapping to each matching column
  for (col in assay_cols) {
    
    x = df[[col]]
    
    for (name in names(assay_map)) {
      for (term in assay_map[[name]]) {
        hit = grepl(term, x, fixed = TRUE)
        x[hit] = name
      }
    }
    
    x = as.factor(x)
    df[[col]] = x
  }
  return(df)
}

convert_scoringTechnique = function(df){
  scoring_map = list(
    
    # Comet assay IV family
    CA4 = c(
      "Automated (Comet assay IV)",
      "Semi-automated image analysis system Comet Assay IV (Instem Perceptive Instruments)",
      "Comet Assay IV",
      "Automated (Comet V v4.1)",
      "Automated (Comet assay 4.1)",
      "Automated (Comet assay IV 4.11/4.2)"
    ),
    
    # Komet software
    Komet6 = c(
      "Automated (Komet6)",
      "Komet6"
    ),
    
    # Metafer
    Metafer = c(
      "Automated (+ MetaFer4)",
      "Metafer",
      "Automated (Metafer)",
      "Automated (MetaFer)"
    ),
    
    # Colourmorph
    clrmorph = c(
      "Colourmorph image analysis system",
      "Automated (Colourmorph image analysis system)"
    ),
    
    # Celigo imaging cytometer
    CeligoS = c(
      "Celigo S imaging cytometer (Nexcelom Bioscience, Lawrence, MA)"
    ),
    
    # Cell Analyzer
    cellAnalyzer6K = c(
      "Automated (Cell Analyzer 6000)",
      "Cell Analyzer 6000"
    ),
    
    # Coulter impedance counting
    Coulter = c(
      "Automated - Coulter counter",
      "Coulter"
    ),
    
    # Classical flow cytometry
    flowCyto = c(
      "Automated, flow cytometric scoring",
      "Flow cytometry",
      "Automated (FACSDiva v8.0.1)",
      "Microflow",
      "MicroFlow",
      "flowCyto",
      "flow cytometry (Becton Dickinson FACScalibur)",
      "flow cytometry (Miltenyi Biotec\r\nMACSQuant Analyzer 10)"
    ),
    
    # Imaging flow cytometry (NEW)
    imageStream = c(
      "Automated (Imagestream)",
      "Imagestream"
    ),
    
    # Fully manual scoring
    manual = c(
      "manual",
      "Manual",
      "manual (fluorescence microscopy)",
      "Fluorescence microscopy"
    ),
    
    # Semi-manual microscope scoring (NEW)
    semiManualMicroscope = c(
      "Semi-Automated (Olympus BX50 Fluorescence Microscope)",
      "Semi-Automated (Olympus Fluorescence Microscope)"
    )
  )

  scoring_cols = names(df)[
    grepl("^(technique for scoring|type of measurement)", names(df), ignore.case = TRUE)
  ]
  
  # Apply mapping to each matching column
  for (col in scoring_cols) {
    
    x = df[[col]]
    
    for (name in names(scoring_map)) {
      for (term in scoring_map[[name]]) {
        hit = grepl(term, x, fixed = TRUE)
        x[hit] = name
      }
    }
    
    x = as.factor(x)
    df[[col]] = x
  }
  return(df)
}

convert_locus_naming = function(df){
  
  locus_map = list(
    TK = c("TK", "Tk", "tk"),
    HPRT = c("HGPRT", "HPRT", "hprt", "hypoxanthine guanine phosphoribosyltransferase"),
    NAKATPase = c("NA/K ATPase"),
    TUB = c("tubulin"),
    RNAsyn = c("RNA synthesis gene")
  )
  
  # Find relevant columns (case-insensitive)
  locus_cols = names(df)[
    grepl("^Gene locus", names(df), ignore.case = TRUE)
  ]
  
  # Apply mapping to each matching column
  for (col in locus_cols) {
    
    x = df[[col]]
    
    for (name in names(locus_map)) {
      for (term in locus_map[[name]]) {
        hit = grepl(term, x, fixed = TRUE)
        x[hit] = name
      }
    }
    
    x = as.factor(x)
    df[[col]] = x
  }
  return(df)
}


#Convert response to FC
convert_to_FC = function(df,
                         conc_col = "Concentration EMS",
                         outcome_col = "Outcome for the treatment") {
  
  df = df %>%
    group_by(`Study index`) %>%
    mutate(
      control_response = .data[[outcome_col]][.data[[conc_col]] == 0],
      control_SD       = SD[.data[[conc_col]] == 0],
      
      FC = .data[[outcome_col]] / control_response,
      
      FC_SD = FC * sqrt(
        (SD / .data[[outcome_col]])^2 +
          (control_SD / control_response)^2
      ),
      
      FC_SEM = FC_SD / sqrt(replicates)
    ) %>%
    ungroup()
  
  df$FC_lower = df$FC - df$FC_SEM
  df$FC_upper = df$FC + df$FC_SEM
  
  return(df)
}

doseres_plot = function(TK6_df,
                        title,
                        assay,
                        conc_col = "Concentration EMS",
                        response_col = "FC",
                        sem_col = "FC_SEM",
                        width = 7,
                        height = 5,
                        save = F,
                        path = NULL) {
  
  # --- Offset zero handling ---
  min_nonzero = min(TK6_df[[conc_col]][TK6_df[[conc_col]] > 0])
  offset_zero = min_nonzero / 4
  
  # --- Legend labels ---
  if ("distinct" %in% names(TK6_df)) {
    
    label_df = TK6_df %>%
      distinct(`Study index`, Author, Year, distinct) %>%
      group_by(Author) %>%
      mutate(n_author = n_distinct(`Study index`)) %>%
      group_by(Author, Year) %>%
      mutate(n_author_year = n_distinct(`Study index`)) %>%
      ungroup() %>%
      mutate(
        legend_label = case_when(
          
          # Only one study by this author total
          n_author == 1 ~ Author,
          
          # Multiple studies by author, but only one in this year
          n_author > 1 & n_author_year == 1 ~ paste0(Author, ", ", Year),
          
          # Multiple studies same author same year
          n_author_year > 1 & !is.na(distinct) ~ paste0(Author, ", ", distinct),
          
          # Fallback if distinct is NA
          TRUE ~ paste0(Author, ", ", Year)
        )
      ) %>%
      select(`Study index`, legend_label)
    
  } else {
    
    # Fallback: use Author only
    label_df = TK6_df %>%
      distinct(`Study index`, Author) %>%
      mutate(legend_label = Author) %>%
      select(`Study index`, legend_label)
    
  }
  
  TK6_df = TK6_df %>%
    left_join(label_df, by = "Study index")
  
  # --- Offset zero ---
  TK6_df = TK6_df %>%
    mutate(
      dose_plot = ifelse(.data[[conc_col]] == 0,
                         offset_zero,
                         .data[[conc_col]])
    )
  
  # --- Log-symmetric whiskers ---
  TK6_df = TK6_df %>%
    mutate(
      response_val = .data[[response_col]],
      sem_val      = .data[[sem_col]],
      logFC  = log10(response_val),
      logSEM = sem_val / (response_val * log(10)),
      ymin   = 10^(logFC - logSEM),
      ymax   = 10^(logFC + logSEM)
    )
  
  log_start = ceiling(log10(min_nonzero))
  log_end   = ceiling(log10(max(TK6_df[[conc_col]])))
  
  breaks_vec = c(offset_zero, 10^seq(log_start, log_end, by = 1))
  whiskerwidth = .2
  
  if(length(breaks_vec) < 4){
    log_start = log_start = ceiling(log2(min_nonzero))
    log_end   = ceiling(log2(max(TK6_df[[conc_col]])))
    
    breaks_vec = c(offset_zero, 2^seq(log_start, log_end, by = 1))
    whiskerwidth = .05
  }
  
  # --- Plot ---
  p = ggplot(
    TK6_df,
    aes(x = dose_plot,
        y = response_val,
        group = `Study index`,
        color = legend_label)
  ) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = whiskerwidth, alpha = 0.7) +
    geom_line() +
    geom_point(size = 2) +
    
    scale_x_continuous(
      trans = "log10",
      breaks = breaks_vec,
      labels = c("0", formatC(breaks_vec[-1], format = "fg"))
    ) +
    scale_y_log10() +
    
    ggtitle(paste(gsub("_", " ", assay), title)) +
    xlab(paste0(conc_col, " (µM)")) +
    ylab("Fold Change") +
    labs(color = "Study") +
    theme_classic()
  
  print(p)
  
  if(save){
    ggsave(
      filename = paste0(path, '/', assay, "_", title, ".png"),
      width = width,
      height = height,
      plot = p,
      dpi = 500
    )
  }
}

#Makes a plot where one can filter on specific variables within the dataframe
doseres_contrast = function(TK6_df,
                            title,
                            assay,
                            conc_col = "Concentration EMS",
                            response_col = "FC",
                            sem_col = "FC_SEM",
                            color_col,
                            shape_col = NULL,
                            linetype_col = NULL,
                            solid_level = NULL,
                            factor_colors = NULL,
                            width = 7,
                            height = 5,
                            save = F,
                            path = NULL) {
  
  # --- Offset zero handling ---
  min_nonzero = min(TK6_df[[conc_col]][TK6_df[[conc_col]] > 0], na.rm = TRUE)
  offset_zero = min_nonzero / 4
  
  TK6_df = TK6_df %>%
    mutate(
      dose_plot = ifelse(.data[[conc_col]] == 0,
                         offset_zero,
                         .data[[conc_col]])
    )
  
  # --- Log-symmetric whiskers ---
  TK6_df = TK6_df %>%
    mutate(
      response_val = .data[[response_col]],
      sem_val      = .data[[sem_col]],
      logFC  = log10(response_val),
      logSEM = sem_val / (response_val * log(10)),
      ymin   = 10^(logFC - logSEM),
      ymax   = 10^(logFC + logSEM)
    )
  
  # --- Axis breaks ---
  log_start = ceiling(log10(min_nonzero))
  log_end   = ceiling(log10(max(TK6_df[[conc_col]], na.rm = TRUE)))
  
  breaks_vec = c(offset_zero, 10^seq(log_start, log_end, by = 1))
  
  if(length(breaks_vec) < 4){
    log_start = log_start = ceiling(log2(min_nonzero))
    log_end   = ceiling(log2(max(TK6_df[[conc_col]])))
    
    breaks_vec = c(offset_zero, 2^seq(log_start, log_end, by = 1))
    whiskerwidth = .05
  }else if(length(breaks_vec) < 6){
    whiskerwidth = .10
  }else{
    whiskerwidth = .2
  }
  
  
  # --- Determine color scale ---
  
  color_is_numeric = is.numeric(TK6_df[[color_col]])
  
  if (color_is_numeric) {
    
    # sort unique numeric values
    unique_vals = sort(unique(TK6_df[[color_col]]))
    
    # convert to ordered factor (preserves numeric order)
    TK6_df[[color_col]] = factor(
      TK6_df[[color_col]],
      levels = unique_vals,
      ordered = TRUE
    )
    
    # generate evenly spaced colors from Spectral palette
    base_cols = RColorBrewer::brewer.pal(9, "YlOrRd")[3:9]
    
    spectral_fun = colorRampPalette(base_cols)
    factor_colors = spectral_fun(length(unique_vals))
    names(factor_colors) = unique_vals
    
  } else {
    
    TK6_df[[color_col]] = as.factor(TK6_df[[color_col]])
    
    if (is.null(factor_colors)) {
      n_levels = length(levels(TK6_df[[color_col]]))
      factor_colors = scales::hue_pal()(n_levels)
      names(factor_colors) = levels(TK6_df[[color_col]])
    }
  }
  
  # --- Shape variable ---
  if (!is.null(shape_col)) {
    TK6_df[[shape_col]] = as.factor(TK6_df[[shape_col]])
  }
  
  # --- Linetype variable ---
  if (!is.null(linetype_col)) {
    
    TK6_df[[linetype_col]] = as.factor(TK6_df[[linetype_col]])
    levels_lt = levels(TK6_df[[linetype_col]])
    
    base_linetypes = c(
      "solid",
      "dashed",
      "dotted",
      "dotdash",
      "longdash",
      "twodash"
    )
    
    if (is.null(solid_level)) {
      
      # assign unique linetype to every level
      linetypes = setNames(
        base_linetypes[seq_along(levels_lt)],
        levels_lt
      )
      
    } else {
      
      non_solid_types = base_linetypes[base_linetypes != "solid"]
      
      linetypes = setNames(
        non_solid_types[seq_along(levels_lt)],
        levels_lt
      )
      
      if (solid_level %in% levels_lt) {
        linetypes[solid_level] = "solid"
      }
      
    }
    
  }
  
  # --- Dynamic aesthetic mapping ---
  
  aes_map = list(
    x = as.name("dose_plot"),
    y = as.name("response_val"),
    group = as.name("Study index"),
    color = as.name(color_col)
  )
  
  if (!is.null(shape_col)) {
    aes_map$shape = as.name(shape_col)
  }
  
  if (!is.null(linetype_col)) {
    aes_map$linetype = as.name(linetype_col)
  }
  
  p = ggplot(
    TK6_df,
    do.call(aes, aes_map)
  )
  
  # --- Errorbars (no linetype inheritance) ---
  
  p = p +
    
    geom_errorbar(
      aes(
        x = dose_plot,
        ymin = ymin,
        ymax = ymax,
        color = .data[[color_col]]
      ),
      width = whiskerwidth,
      alpha = 0.7,
      inherit.aes = FALSE,
      data = TK6_df
    ) +
    
    geom_line(alpha = 0.7) +
    
    geom_point(size = 2) +
    
    scale_x_continuous(
      trans = "log10",
      breaks = breaks_vec,
      labels = c("0", formatC(breaks_vec[-1], format = "fg"))
    ) +
    
    scale_y_log10() +
    
    ggtitle(paste(gsub("_", " ", assay), title)) +
    
    xlab(paste0(conc_col, " (µM)")) +
    
    ylab("Fold Change") +
    
    labs(color = color_col) +
    
    theme_classic()
  
  # --- Color scale ---
  
  p = p +
    scale_color_manual(values = factor_colors)
  
  # --- Shape legend ---
  
  if (!is.null(shape_col)) {
    p = p + labs(shape = shape_col)
  }
  
  # --- Linetype legend ---
  
  if (!is.null(linetype_col)) {
    
    p = p +
      scale_linetype_manual(values = linetypes) +
      labs(linetype = linetype_col)
    
  }
  
  print(p)
  
  # --- Save ---
  
  if(save){
    
    ggsave(
      filename = paste0(path, '/', assay, "_", title, ".png"),
      width = width,
      height = height,
      plot = p,
      dpi = 500
    )
    
  }
}

#Plots dose response curve(s) based on model predictions of the mixed effects model
plot_lmer_dose_response = function(model,
                                   data,
                                   conc_var,
                                   response_var,   # must be log(FC)
                                   group_var = NULL,
                                   sem_var = NULL, # must be SEM on log scale
                                   factor_reference = list(),
                                   grid_length = 200,
                                   offset_fraction = 0.25,
                                   save_plot = FALSE,
                                   save_dir = ".",
                                   width = 6,
                                   height = 5,
                                   dpi = 300) {
  
  # helper: statistical mode
  get_mode = function(x) {
    ux = na.omit(x)
    tab = table(ux)
    names(tab)[which.max(tab)]
  }
  
  # extract model name and clean it
  model_name =
    deparse(substitute(model)) |>
    gsub("^model_", "", x = _) |>
    gsub("^model", "", x = _)
  
  # extract model frame
  mf = model.frame(model)
  
  conc = data[[conc_var]]
  
  min_nonzero = min(conc[conc > 0], na.rm = TRUE)
  offset_zero = min_nonzero * offset_fraction
  
  conc_grid = c(
    0,
    10^seq(log10(min_nonzero),
           log10(max(conc, na.rm = TRUE)),
           length.out = grid_length)
  )
  
  predictors = names(mf)[-1]
  
  base_newdata = data.frame(tmp = conc_grid)
  names(base_newdata)[1] = conc_var
  
  for (var in predictors) {
    
    if (var == conc_var) next
    
    # numeric predictors → mean
    if (is.numeric(data[[var]])) {
      
      base_newdata[[var]] =
        mean(data[[var]], na.rm = TRUE)
      
    }
    
    # factor predictors → reference level OR most frequent level
    else if (is.factor(data[[var]]) || is.character(data[[var]])) {
      
      if (var %in% names(factor_reference)) {
        
        ref_level = factor_reference[[var]]
        
      } else {
        
        ref_level = get_mode(data[[var]])
        
      }
      
      base_newdata[[var]] =
        factor(ref_level,
               levels = levels(as.factor(data[[var]])))
    }
  }
  
  # prediction dataframe
  if (!is.null(group_var)) {
    
    levels_group =
      levels(as.factor(data[[group_var]]))
    
    pred_all =
      dplyr::bind_rows(lapply(levels_group, function(level) {
        
        nd = base_newdata
        
        nd[[group_var]] =
          factor(level,
                 levels = levels(as.factor(data[[group_var]])))
        
        nd$pred =
          predict(model,
                  newdata = nd,
                  re.form = NA)
        
        nd
        
      }))
    
  } else {
    
    pred_all = base_newdata
    
    pred_all$pred =
      predict(model,
              newdata = pred_all,
              re.form = NA)
  }
  
  # prepare plotting data
  plotdata =
    data %>%
    dplyr::mutate(
      dose_plot =
        ifelse(.data[[conc_var]] == 0,
               offset_zero,
               .data[[conc_var]])
    )
  
  pred_all =
    pred_all %>%
    dplyr::mutate(
      dose_plot =
        ifelse(.data[[conc_var]] == 0,
               offset_zero,
               .data[[conc_var]])
    )
  
  # error bars (NOW SIMPLE ± SEM ON LOG SCALE)
  if (!is.null(sem_var)) {
    
    plotdata =
      plotdata %>%
      dplyr::mutate(
        ymin = .data[[response_var]] - .data[[sem_var]],
        ymax = .data[[response_var]] + .data[[sem_var]]
      )
  }
  
  log_start = ceiling(log10(min_nonzero))
  log_end   = ceiling(log10(max(conc)))
  
  breaks_vec =
    c(offset_zero,
      10^seq(log_start, log_end, by = 1))
  
  # initialize plot
  if (!is.null(group_var)) {
    
    p =
      ggplot2::ggplot(plotdata,
                      ggplot2::aes(x = dose_plot,
                                   y = .data[[response_var]],
                                   color = .data[[group_var]]))
    
  } else {
    
    p =
      ggplot2::ggplot(plotdata,
                      ggplot2::aes(x = dose_plot,
                                   y = .data[[response_var]]))
  }
  
  # prediction lines
  if (!is.null(group_var)) {
    
    p =
      p +
      ggplot2::geom_line(data = pred_all,
                         ggplot2::aes(x = dose_plot,
                                      y = pred,
                                      color = .data[[group_var]]),
                         linewidth = 1.3,
                         alpha = 0.9)
    
  } else {
    
    p =
      p +
      ggplot2::geom_line(data = pred_all,
                         ggplot2::aes(x = dose_plot,
                                      y = pred),
                         linewidth = 1.3,
                         alpha = 0.9)
  }
  
  # error bars
  if (!is.null(sem_var)) {
    
    p =
      p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ymin,
                     ymax = ymax),
        width = .15,
        alpha = 0.7,
        linewidth = 0.6)
  }
  
  # datapoints
  p =
    p +
    ggplot2::geom_point(size = 2.2,
                        alpha = 0.9)
  
  # axes
  p =
    p +
    ggplot2::scale_x_continuous(
      trans = "log10",
      breaks = breaks_vec,
      labels =
        c("0",
          formatC(breaks_vec[-1],
                  format = "fg"))
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = conc_var,
      y = paste0("log(", gsub("^log\\(|\\)$", "", response_var), ")"),
      color = group_var
    )
  
  # auto save
  if (save_plot) {
    
    group_part =
      ifelse(is.null(group_var),
             "",
             paste0("_", group_var))
    
    filename =
      paste0(model_name,
             group_part,
             "_predictionplot_logscale.png")
    
    filepath =
      file.path(save_dir, filename)
    
    ggplot2::ggsave(filepath,
                    plot = p,
                    width = width,
                    height = height,
                    dpi = dpi)
  }
  
  return(list(
    plot = p,
    prediction_data = pred_all
  ))
}

#Get data in appropriate format to run a lmm on log transformed responses
prepare_varcomp = function(df, moddf = NULL,
                           out_col = "FC",
                           sem_col = "FC_SEM",
                           sd_col = "FC_SD",
                           rep_col = 'replicates') {
  
  #Check for presence of moddf
  if(is.null(moddf)){
    moddf = df
  }
  
  # Copy columns from moddf into df
  #df[[sem_col]] = moddf[[sem_col]]
  #df[[sd_col]]  = moddf[[sd_col]]
  new_cols = setdiff(names(moddf), names(df))
  df[new_cols] = moddf[new_cols]
  
  # Bias correction
  #biasCorr = -0.5 * (df[[sd_col]] / df[[out_col]])
  biasCorr = -0.5 * log(1 + (df[[sd_col]] / df[[out_col]])^2)
  biasCorr[is.na(biasCorr)] = median(na.omit(biasCorr))
  
  # Determine output column names
  if (out_col == "FC" && sem_col == "FC_SEM") {
    fc_log_col  = "FClog"
    sem_log_col = "FC_SEMlog"
  } else {
    fc_log_col  = paste0(out_col, "log")
    sem_log_col = paste0(sem_col, "log")
  }
  
  # Compute transformed values
  df[[fc_log_col]]  = log(df[[out_col]]) + biasCorr
  
  CV2 <- (df[[sd_col]] / df[[out_col]])^2
  SD <- sqrt(log(1 + CV2))
  
  df[[sem_log_col]] = SD / sqrt(df[[rep_col]])
  
  return(df)
}

doseres_contrast_multi = function(df_list,
                                  df_names = NULL,
                                  title,
                                  assay,
                                  response_col = "FC",
                                  sem_col = "FC_SEM",
                                  color_col,
                                  shape_col = NULL,
                                  linetype_col = NULL,
                                  solid_level = NULL,
                                  width = 10,
                                  height = 5,
                                  save = F) {
  
  
  # =========================
  # 1. NAME HANDLING
  # =========================
  
  if (is.null(df_names)) {
    df_names = names(df_list)
    if (is.null(df_names)) {
      df_names = paste0("Dataset_", seq_along(df_list))
    }
  }
  
  # =========================
  # 2. DETECT CONCENTRATION COLUMN + MERGE
  # =========================
  
  df_list_processed = mapply(function(df, nm) {
    
    # find concentration column via regex
    conc_candidates = grep("^Concentration [A-Za-z]{1,3}$",
                           names(df),
                           value = TRUE)
    
    if (length(conc_candidates) == 0) {
      stop(paste("No valid concentration column found in", nm))
    }
    
    if (length(conc_candidates) > 1) {
      stop(paste("Multiple concentration columns found in", nm))
    }
    
    conc_col = conc_candidates[1]
    
    df %>%
      mutate(
        concentration = .data[[conc_col]],
        dataset = nm
      )
    
  }, df_list, df_names, SIMPLIFY = FALSE)
  
  combined_df = bind_rows(df_list_processed)
  
  # =========================
  # 3. OFFSET ZERO HANDLING
  # =========================
  
  min_nonzero = min(combined_df$concentration[combined_df$concentration > 0], na.rm = TRUE)
  offset_zero = min_nonzero / 4
  
  combined_df = combined_df %>%
    mutate(
      dose_plot = ifelse(concentration == 0,
                         offset_zero,
                         concentration)
    )
  
  # =========================
  # 4. LOG-SYMMETRIC WHISKERS
  # =========================
  
  combined_df = combined_df %>%
    mutate(
      response_val = .data[[response_col]],
      sem_val      = .data[[sem_col]],
      logFC  = log10(response_val),
      logSEM = sem_val / (response_val * log(10)),
      ymin   = 10^(logFC - logSEM),
      ymax   = 10^(logFC + logSEM)
    )
  
  # =========================
  # 5. AXIS BREAKS (GLOBAL)
  # =========================
  
  log_start = ceiling(log10(min_nonzero))
  log_end   = ceiling(log10(max(combined_df$concentration, na.rm = TRUE)))
  
  breaks_vec = c(offset_zero, 10^seq(log_start, log_end, by = 1))
  
  if (length(breaks_vec) < 4) {
    
    log_start = ceiling(log2(min_nonzero))
    log_end   = ceiling(log2(max(combined_df$concentration, na.rm = TRUE)))
    
    breaks_vec = c(offset_zero, 2^seq(log_start, log_end, by = 1))
    whiskerwidth = .05
    
  } else if (length(breaks_vec) < 6) {
    
    whiskerwidth = .10
    
  } else {
    
    whiskerwidth = .2
  }
  
  # =========================
  # 6. GLOBAL LEVELS
  # =========================
  
  # --- COLOR ---
  if (is.numeric(combined_df[[color_col]])) {
    color_levels = sort(unique(combined_df[[color_col]]))
    combined_df[[color_col]] = factor(combined_df[[color_col]],
                                      levels = color_levels,
                                      ordered = TRUE)
  } else {
    color_levels = levels(as.factor(combined_df[[color_col]]))
    combined_df[[color_col]] = factor(combined_df[[color_col]],
                                      levels = color_levels)
  }
  
  base_cols = RColorBrewer::brewer.pal(9, "YlOrRd")[3:9]
  spectral_fun = colorRampPalette(base_cols)
  
  n_levels = length(color_levels)
  
  if (n_levels == 1) {
    factor_colors = spectral_fun(3)[2]
  } else if (n_levels == 2) {
    full_palette = spectral_fun(101)
    factor_colors = full_palette[c(20, 80)]
  } else {
    factor_colors = spectral_fun(n_levels)
  }
  
  names(factor_colors) = color_levels
  
  # --- SHAPE ---
  if (!is.null(shape_col)) {
    shape_levels = levels(as.factor(combined_df[[shape_col]]))
    combined_df[[shape_col]] = factor(combined_df[[shape_col]], levels = shape_levels)
  }
  
  # --- LINETYPE ---
  if (!is.null(linetype_col)) {
    
    linetype_levels = levels(as.factor(combined_df[[linetype_col]]))
    combined_df[[linetype_col]] = factor(combined_df[[linetype_col]], levels = linetype_levels)
    
    base_linetypes = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    
    if (is.null(solid_level)) {
      
      linetypes = setNames(base_linetypes[seq_along(linetype_levels)], linetype_levels)
      
    } else {
      
      non_solid = base_linetypes[base_linetypes != "solid"]
      linetypes = setNames(non_solid[seq_along(linetype_levels)], linetype_levels)
      
      if (solid_level %in% linetype_levels) {
        linetypes[solid_level] = "solid"
      }
    }
  }
  
  # =========================
  # 7. PLOT
  # =========================
  
  aes_map = list(
    x = as.name("dose_plot"),
    y = as.name("response_val"),
    group = as.name("Study index"),
    color = as.name(color_col)
  )
  
  if (!is.null(shape_col)) aes_map$shape = as.name(shape_col)
  if (!is.null(linetype_col)) aes_map$linetype = as.name(linetype_col)
  
  combined_df$dataset = factor(combined_df$dataset, levels = c('EMS', 'MMS', 'ENU'))
  
  p = ggplot(combined_df, do.call(aes, aes_map)) +
    
    geom_errorbar(
      aes(x = dose_plot, ymin = ymin, ymax = ymax, color = .data[[color_col]]),
      width = whiskerwidth,
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
    
    geom_line(alpha = 0.7) +
    geom_point(size = 2) +
    
    scale_x_continuous(
      trans = "log10",
      breaks = breaks_vec,
      labels = c("0", formatC(breaks_vec[-1], format = "fg"))
    ) +
    
    scale_y_log10() +
    
    facet_wrap(~dataset, scales = "free_x") +
    
    scale_color_manual(values = factor_colors) +
    
    ggtitle(paste(gsub("_", " ", assay), title)) +
    xlab("Concentration (µM)") +
    ylab("Fold Change") +
    labs(color = color_col) +
    
    theme_classic()
  
  if (!is.null(shape_col)) {
    p = p + labs(shape = shape_col)
  }
  
  if (!is.null(linetype_col)) {
    p = p + scale_linetype_manual(values = linetypes) +
      labs(linetype = linetype_col)
  }
  
  # =========================
  # 8. SAVE
  # =========================
  
  if (save) {
    ggsave(
      filename = paste0(assay, "_combined_", title, ".png"),
      plot = p,
      width = width,
      height = height,
      dpi = 500
    )
  }
  
  return(p)
}



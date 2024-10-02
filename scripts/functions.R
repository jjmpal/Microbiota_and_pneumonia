#Function to extract tables
extr_coef_table <- function(fit, var_select=""){ 
  cbind(summary(fit)$conf.int, pval = summary(fit)$coefficients[,"Pr(>|z|)"]) %>%
    as.data.frame %>%
    rownames_to_column("Variable") %>%
    rename(est ="exp(coef)", low = "lower .95", high="upper .95") %>%
    filter(str_detect(Variable, var_select))
}  
  

my.kable <- function(my.table) {
  my.table %>%  
    knitr::kable(escape = F, sign=2) %>%
    kable_classic(bootstrap_options = c("hover"), html_font="Arial") %>%
    row_spec(0,bold=TRUE) %>%
    #kable_classic() borders are not generated properly at atlas, extra_css needed
    row_spec(0,bold=TRUE, extra_css = "border-bottom: solid; border-top: solid;
                                     border-bottom-width: thin; border-top-width: thin;") %>%
    row_spec(nrow(my.table), extra_css = "border-bottom: solid; border-bottom-width: thin;" ) #%>%
    #column_spec(1:ncol(my.table), extra_css = "column-cap: 10px;")
  
}

cox_wrapper <- function(data,
                        predictors,
                        covariates,
                        status,
                        time_to_event,
                        alpha_level,
                        normalize,
                        test_ph_assumption) {
  if(normalize) {  
    if(class(data[, predictors]) == "numeric") {
      x <- data[, predictors]
      data[, predictors] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
    } else {
      data[, predictors] <- apply(data[, predictors], 2, FUN = function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T) })
    }
  }
  ## Formulas ***************************
  linear_formulas <- lapply(predictors, function(x) {
    formula_data <- deparse(substitute(data))
    formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(c(covariates,x), collapse = "+"))
    return(formula)
  }) %>% 
    set_names(predictors)
  ## Cox regression *********************
  print("Cox")
  linear_cox_fit <- lapply(linear_formulas, function(x) {
    coxph(as.formula(x), data=data, x=TRUE)
  })
  ## Check PH assumptions ****************
  if(test_ph_assumption) {
    print("PH assumptions")
    ph_assumption <-  lapply(predictors, function(m) {
      ph_test <- cox.zph(linear_cox_fit[[m]])
      p_values <- ph_test$table[, "p"]
      # significant cases
      x <- which(p_values < 1)
      if(length(x) == 0) {
        return(NULL)
      }
      df <- data.frame(feature = m, variable_not_ph = names(x), p_value = p_values[x])
    }) %>% 
      do.call(rbind, .) %>%
      mutate(p_adj = p.adjust(p_value, "BH")) %>% 
      filter(p_value < alpha_level)
  }
  ## Results *****************************
  print("Results")
  results <- lapply(predictors, function(x) {
    df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
    df <- df[nrow(df), ] %>% 
      select(coef, "se(coef)", "z", "Pr(>|z|)") %>% 
      setNames(c("coef", "se_coef", "Z", "p")) 
    df <- df %>% 
      mutate(predictor = x)
  }) %>% 
    do.call(rbind, .) 
  # Multiple testing correction
  results <- results %>%
    mutate(P_adjusted = p.adjust(p, "fdr")) %>% 
    ungroup() %>% 
    group_by(predictor) 
  # Results in neat form for presentation
  neat_results <- results %>% 
    # filter(p == min(p)) %>% 
    ungroup() %>% 
    mutate(HR = round(exp(coef),3)) %>% 
    mutate(HR_lower_95 = round(exp(coef - 1.96*se_coef), 3),
           HR_upper_95 = round(exp(coef + 1.96*se_coef), 3),
           P = round(p, 5),
           Coefficient = round(coef, 3),
           "Coefficient SE" = round(se_coef, 3)) %>% 
    mutate(HR = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>% 
    select(Predictor = predictor, Coefficient, HR, "p","P_adjusted") %>%
    mutate(HR = ifelse(is.na(Coefficient), NA, HR))  %>% 
    filter(P_adjusted < alpha_level) %>% 
    arrange(P_adjusted) %>% 
    setNames(c("Predictor", "Coefficient", "HR","P-value" ,"P (adjusted)"))
  # Results in a form more convenient for further manipulations
  results <- results %>% 
    ungroup %>% 
    mutate(PH = exp(coef)) %>% 
    mutate(p_adj = P_adjusted) %>% 
    mutate(direction = ifelse(coef < 0, "negative", "positive"))
  if(nrow(neat_results) == 0) {
    return(list(results = results))
  }
  if(test_ph_assumption) {
    if(nrow(neat_results) == 0) {
      return(list(results = results, ph_assumption = ph_assumption))
    }
    return(list(neat_results = neat_results, 
                results = results,
                ph_assumption = ph_assumption))
  }
  return(list(neat_results = neat_results, results = results))
}


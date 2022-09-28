#' Estimate complete prevalence
#'
#' @description The function \code{prevEst} is used to estimate complete prevalence using survival and incomplete incidence data
#'
#' @param incidence Incomplete incidence data.frame containing age at diagnosis (named "ageDiag"), year diagnosed (named "yrDiag"),  count of new cases (named "count)
#' @param survival Survival data.frame for each age at prevalence year containing:  age at diagnosis (named "ageDiag"), year diagnosed (named "yrDiag"), observed survival rate (named "survival")
#' @param year Year for which prevalence is to be estimated
#' @param years Years to use for calculating prevalence. Defaults to all years included in incidence dataframe.
#' @param adjust logical. If TRUE, rates will be age adjusted
#' @param grouped_ages logical. If TRUE, multi-year age gruops are used for adjustment.
#' @param groups If grouped ages are use, vector of grouped age levels contain the lowest group in each agre groupe is selected
#' @param incidence_est Estimated incidence data frame for all years containing  age at diagnosis (named "ageDiag"), year diagnosed (named "yrDiag"),  count of new cases (named "count)
#' @param sex_specific character string containing sex ("Male" or "Female" if requesting sex-specific rates)
#'
#' @return Data.frame object
#'
#' @examples
#'
#' prevEst(incidence = incidence.reformat,
#'        survival = survival.nosurvival,
#'        year = 2018,
#'        adjust = T,
#'        standardpopulation = standard.population,
#'        censuspopulation = census.population,
#'        grouped_ages = T,
#'        groups = seq(0,85,5))
#'
#'
#' @export


prevEst <- function (
  incidence,                                               # Incomplete incidence dataframe containing the variables: 1) age at diagnosis, 2) year diagnosed, 3) count
  survival,                                                # Survival dataframe for each age at prevalence year containing: 1) age at diagnosis, 2) year diagnosed, 3) observed survival proportion
  year,                                                    # Complete prevalence year,
  years = NULL,
  adjust = F,                                              # Age-adjust results? If no, ignore subsequent arguments
  grouped_ages = F,                                        # If ages are grouped (i.e. 5-year age bands), then this should be TRUE
  groups = NULL,                                           # Vector of grouped age levels containing the lowest age in each group,
                                                           # (e.g. if using 5-year age bands starting at age = 0, the vector would look like c(0, 5, 10, 15, etc.))
  sex_specific = "Both sexes",                             # Use sex-specific statistics for rates? If so, choose "Male" or "Female"
  incidence_est = NULL                                     # Estimated incidence dataframe for years (very optional)
)

{
  # require(dplyr, tidyr)
  options(dplyr.summarise.inform = FALSE)
  # load(file = "sysdata.rda")                               # Data for age-adjusting
  
  if (is.null(incidence)){
    stop("Incidence dataframe is NULL")
  } else if (is.null(survival)){
    stop("Survival dataframe is NULL")
  } else if (is.null(year)){
    stop("Prevalence year not specified")
  } else if (is.null(years)){
    years <- unique(incidence$yrDiag)
  } 
  
  if(adjust == F) {
    inc <- incidence %>% mutate_all(as.numeric)  %>%
    mutate(yrPrev = year,
           agePrev = (year-yrDiag) + ageDiag)
      
    surv <- survival %>% mutate_all(as.numeric) %>%
      mutate(yrPrev = year,
             agePrev = (year-yrDiag) + ageDiag,
             survival = ifelse(agePrev >=100, 0, survival))
    
    
    prevest <- inc %>%
      dplyr::filter(yrDiag %in% years) %>%
      dplyr::left_join(surv, by = c("ageDiag", "yrDiag", "agePrev")) %>%
      dplyr::mutate(yrPrev = year,
                    final = count*survival) %>%
      dplyr::group_by(agePrev) %>%
      dplyr::summarise(prevalence = round(sum(final, na.rm=TRUE))) %>%
      dplyr::ungroup() %>%
      dplyr::arrange()
    
  } else {
    
    if(grouped_ages == T){
      
      prevest<- prevest %>%
        dplyr::mutate(agePrev = as.character(cut(as.numeric(agePrev), c(groups, max(groups)*2), include.lowest = F, right = F, labels = groups))) %>%
        dplyr::group_by(agePrev) %>%
        dplyr::summarise(prevalence  = sum(prevalence)) %>%
        ungroup()
      
      if(sex_specific != "Both sexes") {
        census.population <- census.population %>%
          dplyr::mutate(sex = tolower(sex)) %>%
          dplyr::filter(sex == tolower(sex_specific))
      }
      
      census.population <- census.population %>%
        dplyr::mutate(age = as.character(cut(as.numeric(age), c(groups, max(groups)*500), include.lowest = F, right = F, labels = groups))) %>%
        dplyr::group_by(age) %>%
        dplyr::summarise(census_population = sum(as.numeric(census_population))) %>%
        dplyr::ungroup() %>%
        arrange(as.numeric(age))
      
      
      standard.population <- standard.population %>%
        dplyr::mutate(age = as.character(cut(as.numeric(age), c(groups, max(groups)*500), include.lowest = F, right = F, labels = groups))) %>%
        dplyr::group_by(age) %>%
        dplyr::summarise(standard_population = sum(as.numeric(standard_population))) %>%
        dplyr::ungroup() %>%
        arrange(as.numeric(age))
    }
    
    prevest <- prevest %>%
      mutate(prevalence = ifelse(agePrev >= 85, sum(.$prevalence[which(.$agePrev >= 85)]), prevalence)) %>%
      filter(agePrev <= 85) %>%
      left_join(census.population, by = c("agePrev" = "age")) %>%
      left_join(standard.population, by = c("agePrev" = "age")) %>%
      group_by(agePrev) %>%
      mutate(crude_rate = prevalence/census_population*100000,
             weights = standard_population/sum(standard.population$standard_population),
             adjusted_rate = crude_rate*weights,
             adjusted_lci = asht::wspoissonTest(prevalence, w = as.numeric(weights)/as.numeric(census_population), wmtype = "tcz", mult=100000)$conf.int[[1]],
             adjusted_uci = asht::wspoissonTest(prevalence, w = as.numeric(weights)/as.numeric(census_population), wmtype = "tcz", mult=100000)$conf.int[[2]]) %>%
      ungroup() %>%
      arrange(as.numeric(agePrev))
  }
  
  return(prevest %>% mutate_all(as.numeric))  
}

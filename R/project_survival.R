#' Project survival data for use in prevEst
#'
#' @param data Survival dataframe for which prevalence is to be projected
#' @param ages Ages to be included
#' @param observation.years A vector of years to build the model on
#' @param projection.years A vector of years to predict survival for
#' @param Observed.Surv Logical
#' @param Expected.Surv Logical
#' @param assumption Character, either "population" or "nosurvival"
#' @param life.table SEER life table
#' @param names Named vector pointing to column names
#' @param keepExtraCols Logical
#' @return Predicted survival including \code{projection.years}
#' @examples
#'
#' project_survival(data=raw.survival,
#'                  ages= 0:85,
#'                  years=1975:2018,
#'                  final.year=2030,
#'                  observation.years = NULL,               # A vector of years to build the model on
#'                  projection.years = NULL,
#'                  assumption = "population",
#'                  life.table = life.tables,
#'                  years.observed.surv=20,
#'                  names = c("ageDiag"="Age_recode_with_single_ages_and_85",
#'                            "yrPrev"="yrPrev",
#'                            "yrDiag"="Year_of_diagnosis",
#'                            "Observed"="Observed"))
#'
#' @seealso [prevEst::project_incidence()]
#' @export

project_survival <- function(data,                # Case listing survival data, having either 0/1 or Dead/Alive as indicators
                             ages,
                             observation.years = NULL,               # A vector of years to build the model on
                             projection.years = NULL,
                             Observed.Surv=TRUE,
                             Expected.Surv=FALSE,
                             assumption="nosurvival",
                             life.table=NULL,
                             names=c("ageDiag"="ageDiag",
                                     "yrPrev"="yrPrev",
                                     "yrDiag"="yrDiag",
                                     "Observed"="observed"), # A list of names containing 1) age, 2) year, and 3) survival, of the form list("age" = ..., "year" = ..., etc.)
                             keepExtraCols=FALSE
){
  require(dplyr)
  options(dplyr.summarise.inform = FALSE)
  
  if(is.null(observation.years)) {
    stop("Please specify years for observations.")
  }
  if(is.null(projection.years)) {
    stop("Please specify years for projections.")
  }
  
  new <- data.frame(ageDiag  = as.numeric(gsub("\\D", "", data[[names[["ageDiag"]]]])),
                    yrPrev = as.numeric(gsub("\\D", "", data[[names[["yrPrev"]]]])),
                    yrDiag = as.numeric(gsub("\\D", "", data[[names[["yrDiag"]]]])),
                    survival = as.numeric(data[[names[["Observed"]]]]))
  
  first.year <- min(new$yrDiag)
  final.year <- as.numeric(max(projection.years))
  years.observed.surv = length(observation.years)
  
  if(keepExtraCols==TRUE) {
    new <- new %>%
      bind_cols(data %>% select(-names))
  }
 if(all(is.na(new$survival)) & is.null(life.table)){
   stop("No survival data provided")
 }
 if(all(is.na(new$survival)) & !is.null(life.table)){
    new <- new %>% 
      mutate(period = as.character(yrPrev - yrDiag)) %>%
      left_join(life.table, by = c("ageDiag", "period")) %>% 
      filter(period <= length(observation.years) & 
             ageDiag %in% ages) %>%
      fill(expected, .direction = "downup") %>%
      group_by(ageDiag, period) %>%
      arrange(ageDiag, period) %>% 
      mutate(survival = cumprod(expected)) %>%
      ungroup() %>%
      mutate_all(as.numeric) %>%
      select(-expected)
}
  cat("Projecting ", length(projection.years), " years of survival for ",min(projection.years),"-",max(projection.years), "\n", sep = "")
  
  new <- new  %>%
    drop_na() %>%
    filter(yrDiag %in% observation.years) %>%
    dplyr::mutate( period=yrPrev-yrDiag,
                   agePrev = ageDiag+period,
                   survival = case_when(survival  > 1 ~ survival /100,
                                        TRUE ~ survival))
  proj.surv <- new %>%
    dplyr::group_by(period) %>%
    tidyr::nest() %>%
    dplyr::mutate(predicted_survival=purrr::map(data, ~modelr::add_predictions(data=expand.grid(yrDiag = projection.years,
                                                                                                ageDiag = ages),
                                                                               lm(survival ~ as.numeric(yrDiag) + as.numeric(ageDiag), data = .x),
                                                                               var="survival"))) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols=predicted_survival) %>%
    dplyr::mutate(yrPrev=yrDiag+period,
                  agePrev=ageDiag+period,
                  survival=case_when(as.numeric(survival)<0~0,
                                     as.numeric(survival)>1~1,
                                     TRUE~as.numeric(survival)))
  
  survival.merged <- dplyr::bind_rows(new, proj.surv) %>%
    # dplyr::mutate(survival=case_when(period==0~1,
    #                                   TRUE~survival)) %>%
    dplyr::filter(yrDiag >= first.year) %>% 
    dplyr::mutate_all(as.numeric) %>%
    dplyr::arrange(ageDiag, yrDiag, period)
  
  # Create skeleton dataframe and left join to "new" dataframe to handle missing values. Will default to survival = 1 if data mising.
  skeleton <- tidyr::expand_grid(ageDiag  = ages,
                                 yrDiag = first.year:final.year,
                                 yrPrev = first.year:final.year) %>%
    dplyr::mutate(period=yrPrev-yrDiag,
                   agePrev=ageDiag+period) %>%
    dplyr::arrange(ageDiag , yrDiag)  %>%
    dplyr::mutate_all(as.numeric) %>%
    dplyr::filter(period >= 0)
  
  
  
  if (assumption=="population") {
    if (is.null(life.table)) {
      stop("Life table must be provided for population survival")
    }
    else {
      
      message("Applying population-level survival \n")
      
      if(any(!names(life.table) %in% c("period", "ageDiag", "expected"))){
        stop("Life tables must contain 'period', 'ageDiag', and 'expected' columns \n")
      }
      
      
      life.table <- life.table %>% 
        dplyr::mutate_all(as.numeric) %>%
        dplyr::select(period, ageDiag, expected) %>%
        dplyr::arrange(desc(period)) %>%
        dplyr::filter(period <= years.observed.surv & 
               ageDiag %in% ages)
      
      full.survival.temp1 <- skeleton  %>%
        dplyr::left_join(new, by=names(skeleton)) %>%
        dplyr::mutate_all(as.numeric) %>%
        dplyr::left_join(life.table %>% mutate_all(as.numeric), by = c("ageDiag", "period"))  %>%
        dplyr::mutate(expected = case_when(agePrev>=100 ~ 0,
                                           TRUE~expected),
                      survival = case_when(!is.na(survival)~survival,
                                           TRUE~expected)) 
      
      full.survival.temp2 <- full.survival.temp1  %>%
        filter(period >= (years.observed.surv)) %>%
        arrange(period) %>%
        dplyr::group_by(ageDiag) %>%
        dplyr::mutate(survival = cumprod(survival)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(ageDiag,yrDiag) 
      
      full.survival <- full.survival.temp1 %>%
        dplyr::filter(period <= years.observed.surv) %>%
        dplyr::bind_rows(full.survival.temp2 %>% filter(period > years.observed.surv)) %>%
        dplyr::arrange(ageDiag, period) %>%
        dplyr::select(-expected)
      
    }
  }
  if (assumption=="nosurvival") {
    cat("Applying no survival assumption")
    full.survival <- skeleton  %>%
      dplyr::left_join(new,by=names(skeleton)) %>%
      dplyr::mutate(survival=case_when(period > years.observed.surv ~ 0, TRUE ~ survival))
  }
  
  final <- full.survival %>%
    dplyr::mutate(survival=case_when(survival>1~survival/100,
                                     survival<=0~0,
                                     agePrev>=100~0,
                                     TRUE~survival))  %>%
    dplyr::group_by(ageDiag, yrPrev) %>%
    dplyr::arrange(ageDiag, period, yrPrev) %>%
    tidyr::fill(survival, .direction = "downup") %>%
    dplyr::select(ageDiag,agePrev, yrDiag, yrPrev, period , survival) %>%
    dplyr::filter(period >= 0)
  
  return(as.data.frame(final))
}
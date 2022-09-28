#' Estimate complete prevalence
#'
#' @param incidence Incomplete incidence data.frame containing age at diagnosis (named "ageDiag"), year diagnosed (named "yrDiag"),  count of new cases (named "count)
#' @param incidence_est Estimated incidence data frame for smaller area containing  all years containing  age at diagnosis (named "ageDiag"), year diagnosed (named "yrDiag"),  count of new cases (named "count)
#' @param regYr Years to do regression on ie. c(2000:2010)
#' @param durationYr Years that require estimates from regression
#' @return An incidence dataframe
#' @examples
#' regprev(incidence = incomplete_incidence,
#         incidence_est = complete_incidence,
#'        regYr = c(2001:2017),
#'        durationYr = c(1975:2000))
#'
#' @export


regPrev <- function(
  incidence,
  incidence_est,
  regYr = NULL,
  durationYr = NULL
){
  require(dplyr,tidyr,purrr,modelr)

  incidence <- incidence %>%
    select(c("ageDiag", "yrDiag","count")) %>%
    dplyr::mutate_all(as.numeric)

  incidence_est <- incidence_est %>%
    select(c("ageDiag", "yrDiag","count")) %>%
    dplyr::mutate_all(as.numeric)

  idf <- full_join(incidence, incidence_est, by = c("ageDiag", "yrDiag"))
  
  if(is.null(durationYr) & is.null(regYr)){
    message("Guessing durationYr and regYr as they are not specified. \n")
    # Guess years for regression
    # Linear regression and resulting diagnostic statistics
    regprev <- idf %>%
      arrange(ageDiag, yrDiag) %>%
      dplyr::group_by(ageDiag) %>%
      tidyr::nest() %>%
      dplyr::mutate(data = purrr::map(data, function(x) x %>% mutate(count.y = case_when(all(is.na(count.y)) ~ 0,
                                                                                         T ~ count.y),
                                                                     count.x = case_when(all(is.na(count.x)) & all(is.na(count.y)) ~ 0,
                                                                                         all(is.na(count.x)) ~ count.y,
                                                                                         T ~ count.x))),
                    model_data = purrr::map(data, function(x) drop_na(x)),
                    model = purrr::map(model_data, function(x) lm(count.x ~ count.y, data = x)),
                    predicted_incidence = purrr::map2(data,
                                               model,
                                               ~suppressWarnings(modelr::add_predictions(data = as.data.frame(.x  %>% filter(is.na(count.x))),
                                                                        model = .y,
                                                                        var = "count_pred"))),
                    complete_data = purrr::map2(.x = model_data, .y = predicted_incidence,
                                                ~mutate(bind_rows(.x, .y)))) %>%
      dplyr::select(-c(data, model_data, model, predicted_incidence)) %>%
      tidyr::unnest(cols=c(ageDiag, complete_data)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(count = case_when(is.na(count.x) & (count_pred <=0 | is.na(count_pred)) ~ 0, 
                                      is.na(count.x) ~ round(count_pred),
                                      T ~ round(count.x))) %>%
      dplyr::select(ageDiag, yrDiag, count) 

  } else {

    # Linear regression and resulting diagnostic statistics
    regprev <- idf %>%
      arrange(ageDiag, yrDiag) %>%
      dplyr::group_by(ageDiag) %>%
      tidyr::nest() %>%
      dplyr::mutate(data = purrr::map(data, function(x) x %>% mutate(count.y = case_when(all(is.na(count.y)) ~ 0,
                                                                                         T ~ count.y),
                                                                     count.x = case_when(all(is.na(count.x)) & all(is.na(count.y)) ~ 0,
                                                                                         all(is.na(count.x)) ~ count.y,
                                                                                         T ~ count.x))),
                    model = purrr::map(data, function(x) lm(count.x ~ count.y, data = x %>% filter(yrDiag %in% regYr))),
                    predicted_incidence = map2(data,
                                               model,
                                               ~modelr::add_predictions(data = as.data.frame(.x  %>% filter(yrDiag %in% durationYr)),
                                                                        model = .y,
                                                                        var = "count_pred"))) %>%
      dplyr::select(-c(data, model)) %>%
      tidyr::unnest(cols=predicted_incidence) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(count = ifelse(is.na(count.x), count_pred, count.x)) %>%
      dplyr::select(ageDiag, yrDiag, count) %>%
      dplyr::bind_rows(idf  %>% filter(yrDiag %in% regYr) %>% dplyr::rename(c("count"="count.x")) %>%  dplyr::select(ageDiag, yrDiag, count)) %>%
      dplyr::arrange(ageDiag, yrDiag) %>%
      dplyr::mutate(count = ifelse(as.numeric(count) <= 0, 0, round(count)))
    
}
  return(regprev)

  }

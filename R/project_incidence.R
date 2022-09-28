#' Project incidence for unobserved year
#'
#' @param incidence Incidence dataframe for which prevalence is to be projected
#' @param projection.years Incidence dataframe for which prevalence is to be projected
#' @param method GLM family to be base predictions off of
#' @return Predicted survival proportions for \code{years}
#' @return Projected incidence dataframe including projection years
#' @examples
#' project_incidence(data=complete_incidence,
#'                   final.year=2030)
#' @seealso [project_survival()]
#' @export


project_incidence <- function(data,
                              projection.years = NULL,
                              method = "poisson") {

  require(dplyr)
  options(dplyr.summarise.inform = FALSE)


  # Searches for columns containing age, year, incidence, then renames them to be used later.
  incidence <- data %>%
    dplyr::select(c("ageDiag", "yrDiag","count"))

  if(is.null(projection.years)) {
    stop("Please specify years for projections. \n")
  }

  cat("Projecting incidence for",min(projection.years),"-",max(projection.years), "\n")

  proj.inc <- incidence %>%
    dplyr::group_by(ageDiag  ) %>%
    tidyr::nest() %>%
    dplyr::mutate(predicted_incidence=purrr::map(data, ~modelr::add_predictions(data=data.frame(yrDiag  = projection.years),
                                                                                glm(count ~ as.numeric(yrDiag) , family=method, data = .x),
                                                                                var="count",type="response"))) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols=predicted_incidence) %>%
    dplyr::mutate(count=case_when(as.numeric(count)<0~0,
                                  TRUE~round(as.numeric(count),0))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(yrDiag %in% projection.years & !(yrDiag %in% incidence$yrDiag))

  incidence.merged <-  dplyr::bind_rows(incidence,proj.inc)
}

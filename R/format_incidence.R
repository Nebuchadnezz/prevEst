#' Format incidence data
#'
#' @description This function formats incidence data for easier/more streamlined use with functions included in this package. It takes new casesand summarizes them by age and year of diagnosis.
#' @param data Incidence dataframe
#' @param ages Numeric vector
#' @param years Numeric vector
#' @param names Named character vector
#' @param keepExtraCols Logical
#' @details The prevEst function requires both a properly formatted incidence and survival data. This function, the counterpart to [format_survival()], is designed 
#' to take SEER-like incidence data and format it to work more easily with the [prevEst()] function. 3 columns are necessary: a column for 1) age at diagnosis, 2) year of diagnosis, 
#' and 3) the reported incidence for that combination of the two. While these functions aren't necessary, they help wrap some simple transformation steps.
#' @return
#' @examples
#'
#' format_survival(incDf,
#'                 ages = c(0:85),
#'                 years = c(2004:2018),
#'                 names = c("ageDiag" = "age",
#'                           "yrDiag" = "year",
#'                           "incidence" = "count")
#'                 keepExtraCols = FALSE)
#'
#'
#' @seealso [format_survival()] The sister function that formats survival data
#' @export

format_incidence <- function(data,                 # A dataframe of counts for each unique combination of age and year
                             ages,                # A vector of ages that you'd like included in the output
                             years,               # A vector of years that you'd like included in the output
                             names = c("ageDiag" = "age",
                                       "yrDiag" = "year",
                                       "incidence" = "count"), # A vector of names containing 1) age, 2) year, and 3) counts, of the form list("age" = ..., "year" = ..., etc.)
                             keepExtraCols=FALSE
)  {
  require(dplyr,tidyr)
  options(dplyr.summarise.inform = FALSE)
  new <- data.frame(ageDiag = data[[names[["ageDiag"]]]],
                    yrDiag = as.numeric(data[[names[["yrDiag"]]]]),
                    inc = as.numeric(data[[names[["incidence"]]]])) %>%
    dplyr::mutate(ageDiag = as.numeric(gsub("\\D", "", ageDiag)))

  if( keepExtraCols==TRUE) {
    new <- new %>%
      bind_cols(data %>% select(-names))
  }

skeleton <- tidyr::expand_grid(ageDiag = ages,
                               yrDiag =  as.numeric(years),
                        inc = as.numeric(0)) %>%
  dplyr::arrange(ageDiag, yrDiag)

final <-  dplyr::left_join(skeleton, new, by = c("ageDiag", "yrDiag")) %>%
  dplyr::mutate(count = ifelse(inc.y == 0, inc.x, inc.y)) %>%
  dplyr::select(-c(inc.y, inc.x)) %>%
  dplyr:: arrange(ageDiag, yrDiag)

return(as.data.frame(final))
}

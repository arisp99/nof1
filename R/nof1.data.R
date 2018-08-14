#' Make an N of 1 object containing data, priors, and a jags model file
#'
#' @param Y Outcome of the study. This should be a vector with length of total number of observations.
#' @param Treat Treatment indicator vector with same length as the outcome.
#' @param baseline baseline Treatment name. This serves as a baseline/placebo when comparing different treatments.
#' @param ncat Number of categories. Used in ordinal models.
#' @param response Type of outcome. Can be normal, binomial, poisson or ordinal.
#' @param Time parameter used for modelling splines. Still under development.
#' @param knots parameter used for modelling splines. Still under development.
#' @param alpha.prior Prior for the intercept of the model.
#' @param beta.prior Prior for the treatment coefficient.
#' @param gamma.prior Prior for modelling splines. Still under development.
#' @param dc.prior Prior for the length between cutpoints. Used only for ordinal logistic models.
#' @param c1.prior Prior for the first cutpoint. Used only for ordinal logistic models.
#' @param rho.prior Prior for the correlated error model. Still under development.
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response and wishart for multinomial response. It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm, dwish) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter. For wishart distribution, the last two parameter would be the scale matrix and the degrees of freedom.
#' @return Creates list of variables that are used to run the model using \code{\link{nof1.run}}
#' \item{Y}{Outcome}
#' \item{Treat}{Treatment}
#' \item{baseline}{Baseline variable}
#' \item{ncat}{Number of categories for ordinal response}
#' \item{nobs}{Total number of observations in a study}
#' \item{Treat.name}{Treatment name besides baseline treatment}
#' \item{response}{The type of response variable}
#' \item{priors}{Priors that the code will be using. Default priors are used if prior was not specified}
#' \item{code}{Rjags model file code that is generated using information provided by the user. To view model file inside R, use \code{cat(nof1$code).}}
#' @examples
#' ###Blocker data example
#' laughter
#' Y <- laughter$Y
#' Treat <- laughter$Treat
#' nof1 <- nof1.data(Y, Treat, ncat = 11, baseline = "Usual Routine", response = "ordinal")
#' str(nof1)
#' cat(nof1$code)
#' @export

nof1.data <- function(Y, Treat, baseline = "baseline", ncat = NULL, response = NULL, Time=NULL, knots = NULL,
                      alpha.prior = NULL, beta.prior = NULL, gamma.prior = NULL, dc.prior = NULL, c1.prior = NULL,
                      rho.prior = NULL, hy.prior = NULL){

  if(response == "ordinal"){
    if(is.null(ncat)){
      stop("ncat (number of categories) must be entered for ordinal response")
    }
  }
  nobs <- length(Y)

  if(!baseline %in% Treat){
    stop("baseline treatment name is not in Treat")
  }

  Treat <- gsub(" ", "\\_", Treat)
  baseline <- gsub(" ", "\\_", baseline)

  Treat.name <- unique(Treat)
  Treat.name <- Treat.name[Treat.name != baseline]

  nof1 = list(Y = Y, Treat = Treat, baseline = baseline, ncat = ncat, nobs = nobs, Treat.name = Treat.name, response = response)

  # for correlated model, not used
  if(!is.null(Time)){
    cen.Time <- (Time - mean(Time, na.rm = TRUE)) / sd(Time, na.rm = TRUE)
    nof1$Time = cen.Time
  }

  # for splines, not used
  if(!is.null(knots)){
    cen.knots <- (knots - mean(Time, na.rm = TRUE))/ sd(Time, na.rm = TRUE)
    BS <- bs(cen.Time, knots = cen.knots)
    nof1$BS <- BS
    nof1$knots <- knots
  }

  for(i in Treat.name){
    nam <- paste("Treat_", i, sep = "")
    nam <- assign(nam, as.numeric(Treat == i))
    nof1[[ paste("Treat_", i, sep = "")]] <- nam
  }

  prior.param <- list(response = response, dc.prior = dc.prior, c1.prior = c1.prior, alpha.prior = alpha.prior, beta.prior = beta.prior, gamma.prior = gamma.prior, hy.prior = hy.prior, rho.prior = rho.prior)
  prior.data <- nof1.prior.default(prior.param)

  nof1 <- c(nof1, prior.data)

  code <- nof1.rjags(nof1)
  nof1$code <- code

  class(nof1) <- "nof1.data"
  nof1
}

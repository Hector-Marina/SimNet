#' Estimate transmission parameters from infection-status data
#'
#' @description
#' `TPest()` estimates the infection (`beta`) and recovery (`gamma`) rates from
#' the temporal infection-status information stored in a `simPOP` object.
#'
#' For each simulation replicate, the function identifies new infections and
#' recoveries between consecutive simulation steps. Transmission parameters are
#' then estimated independently using intercept-only negative binomial
#' regression models.
#'
#' The infection rate is estimated using `log(S*I/N)` as an offset, where `S`
#' and `I` are the numbers of susceptible and infectious individuals,
#' respectively, and `N=S+I`. The recovery rate is estimated using `log(I)` as
#' an offset.
#'
#' The expected infection duration is calculated as the inverse of the
#' estimated recovery rate.
#'
#' @param POP A `simPOP` object containing temporal infection-status
#' information.
#'
#' @param verbose Logical indicating whether model-fitting errors are printed.
#' Default is `TRUE`.
#'
#' @return A data frame containing the estimated infection rate
#' (`beta_estimate`), recovery rate (`gamma_estimate`), and expected infection
#' duration (`Inf_duration`) for each simulation replicate.
#'
#' @examples
#' \dontrun{
#' TPest(POP=POP)
#' }
#'
#' @export
TPest=function(POP=NULL, verbose=TRUE) {

  # Validate simPOP object
  if(!inherits(POP, "simPOP")) stop("`POP` must be class `simPOP`.")
  if(!validate(POP)) stop("Invalid `simPOP` object.")

  # Validate verbose
  if(!is.logical(verbose) || length(verbose)!=1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.")

  # Output table
  TPtable=data.frame(beta_estimate=rep(NA_real_, POP@rep),
                    gamma_estimate=rep(NA_real_, POP@rep),
                    Inf_duration  =rep(NA_real_, POP@rep))

  # Estimate transmission parameters per replicate
  for(i in seq_along(POP@statusM)){

    Status=POP@statusM[[i]]

    # Status at consecutive simulation steps
    Status1=Status[,-ncol(Status),drop=FALSE]
    Status2=Status[,-1,drop=FALSE]

    # Number of susceptible and infectious individuals
    S_time=colSums(Status1==0)
    I_time=colSums(Status1==1)

    # Population at risk
    N_time=S_time+I_time

    df=data.frame(
      S_time=S_time,
      I_time=I_time,
      N_time=N_time,
      new_infections=colSums(Status1==0 & Status2==1),
      new_recoveries=colSums(Status1==1 & Status2!=1)
    )

    beta_hat=NA_real_
    gamma_hat=NA_real_

    #- Estimate infection rate
    df_beta=df[df$S_time>0 &df$I_time>0 & df$N_time>0,]

    if(nrow(df_beta)>0 && sum(df_beta$new_infections)>0){

      fit_beta=function(){
        MASS::glm.nb(
          new_infections ~ 1 +
            offset(log(S_time) + log(I_time) - log(N_time)),
          data=df_beta#,
          #control=stats::glm.control(maxit=500)
        )
      }

      if(verbose){
        model_foi=try(fit_beta(), silent=FALSE)
      }else{
        model_foi=try(
          suppressWarnings(fit_beta()),
          silent=TRUE
        )
      }

      if(!inherits(model_foi,"try-error") &&
         isTRUE(model_foi$converged))
        beta_hat=exp(stats::coef(model_foi)[1])
    }

    #- Estimate recovery rate
    df_gamma=df[df$I_time>0,]

    if(nrow(df_gamma)>0 && sum(df_gamma$new_recoveries)>0){

      fit_gamma=function(){
        MASS::glm.nb(
          new_recoveries ~ 1 +
            offset(log(I_time)),
          data=df_gamma#,
          #control=stats::glm.control(maxit=500)
        )
      }

      if(verbose){
        model_cure=try(fit_gamma(), silent=FALSE)
      }else{
        model_cure=try(
          suppressWarnings(fit_gamma()),
          silent=TRUE
        )
      }

      if(!inherits(model_cure,"try-error") &&
         isTRUE(model_cure$converged))
        gamma_hat=exp(stats::coef(model_cure)[1])
    }

    # Expected infection duration
    Inf_duration=ifelse(!is.na(gamma_hat) && gamma_hat>0,
                        1/gamma_hat,NA_real_)

    # Store estimates
    TPtable[i,]=c(beta_hat,gamma_hat,Inf_duration)
  }

  return(TPtable)
}

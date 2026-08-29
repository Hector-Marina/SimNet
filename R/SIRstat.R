#' Simulate a stochastic status transition
#'
#' @description
#' `SIRstat()` simulates the stochastic transition of individuals between two
#' states during a single simulation step.
#'
#' Individuals with status `stat1` at the previous simulation step are sampled
#' independently with probability `par`. Individuals for which a transition
#' occurs are assigned status `stat2` at the current step, while the remaining
#' individuals retain status `stat1`.
#'
#' @param Status A numeric matrix containing the status of each individual
#' across simulation steps. Rows represent individuals and columns represent
#' simulation steps.
#'
#' @param step A positive integer indicating the current simulation step.
#' Transitions are sampled using the individual statuses at `step - 1`.
#'
#' @param stat1 A numeric value indicating the status from which individuals
#' may transition.
#'
#' @param stat2 A numeric value indicating the status to which individuals
#' may transition.
#'
#' @param par A numeric value between `0` and `1` indicating the probability
#' of transition from `stat1` to `stat2`.
#'
#' @return A numeric status matrix with the simulated transitions recorded at
#' the specified simulation step.
#'
#' @keywords internal
#'
#' @examples
#' set.seed(123)
#'
#' Status=matrix(0, nrow=10, ncol=5)
#' Status[1:5,1]=1
#'
#' # Simulate recovery from infected (1) to susceptible (0)
#' Status=SIRstat(Status=Status, step=2, stat1=1, stat2=0, par=0.2)
#' Status
#'
#' @export
SIRstat=function(Status=NULL, step=NULL, stat1=NULL, stat2=NULL, par=NULL) {

  # Validate Status
  if(!is.matrix(Status) || !is.numeric(Status))
    stop("`Status` must be a numeric matrix.")

  # Validate step
  if(!is.numeric(step) || length(step)!=1L || !is.finite(step) ||
     step%%1!=0 || step<2 || step>ncol(Status))
    stop("`step` must be an integer between 2 and the number of columns in `Status`.")

  # Validate statuses
  if(!is.numeric(stat1) || length(stat1)!=1L || !is.finite(stat1))
    stop("`stat1` must be a single numeric value.")

  if(!is.numeric(stat2) || length(stat2)!=1L || !is.finite(stat2))
    stop("`stat2` must be a single numeric value.")

  if(stat1==stat2)
    stop("`stat1` and `stat2` must be different.")

  # Validate transition probability
  if(!is.numeric(par) || length(par)!=1L || !is.finite(par) ||
     par<=0 || par>=1)
    stop("`par` must be a numeric value between 0 and 1.")

  # Individuals eligible for transition
  IDsStat=which(Status[,step-1]==stat1)

  if(length(IDsStat)>0){

    # Sample stochastic transition
    Rstat=stats::rbinom(length(IDsStat), size=1, prob=par)

    # Individuals remaining in stat1
    if(any(Rstat==0))
      Status[IDsStat[Rstat==0],step]=stat1

    # Individuals transitioning to stat2
    if(any(Rstat==1))
      Status[IDsStat[Rstat==1],step]=stat2
  }

  return(Status)
}

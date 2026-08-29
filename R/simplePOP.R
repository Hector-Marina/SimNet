#' Simulate a simple population
#'
#' @description
#' `simplePOP` creates a simple population for stochastic transmission
#' simulations. The population contains a specified number of individuals,
#' simulation steps, and replicates, without additional individual-level
#' information.
#'
#' All individuals are initially assigned as susceptible. Contact networks and
#' transmission parameters can subsequently be added using other `SimNet`
#' functions.
#'
#' @param NIDs A positive integer indicating the number of individuals in the
#' simulated population (Default: *NIDs=200*).
#' @param steps A positive integer indicating the number of transmission steps
#' to simulate (Default: *steps=100*).
#' @param rep A positive integer indicating the number of simulation replicates
#' (Default: *rep=100*).
#'
#' @return A `simPOP` object containing the simulated population. The
#' `statusM` slot contains one status matrix per replicate, with individuals in
#' rows and simulation steps in columns. All individuals are initially coded as
#' susceptible (`0`).
#'
#' @examples
#' # Create a population of 200 individuals
#' POP=simplePOP(NIDs = 200, steps = 100, rep = 100)
#' POP
#'
#' @export

simplePOP=function(NIDs=200, steps=100, rep=100) {

  if (!is.numeric(NIDs) || length(NIDs) != 1 || !is.finite(NIDs) ||
      NIDs < 1 || NIDs %% 1 != 0) {
    stop("`NIDs` must be a positive integer.")
  }

  if (!is.numeric(steps) || length(steps) != 1 || !is.finite(steps) ||
      steps < 1 || steps %% 1 != 0) {
    stop("`steps` must be a positive integer.")
  }

  if(steps < 2) stop("`POP` must contain at least two simulation steps.")

  if (!is.numeric(rep) || length(rep) != 1 || !is.finite(rep) ||
      rep < 1 || rep %% 1 != 0) {
    stop("`rep` must be a positive integer.")
  }

  POP=methods::new("simPOP",
          method = "simplePOP",
          IDs=1:NIDs,
          steps = steps,
          rep=rep,
          statusM = lapply(1:rep, function(i) {
            matrix(0, nrow = NIDs, ncol = steps)})
  )

  # Validate resulting object
  if(!validate(POP))
    stop("The simulation produced an invalid `simPOP` object.")

  return(POP)
}

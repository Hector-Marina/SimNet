#' Simulate a homogeneous-mixing contact network
#'
#' @description
#' `uniCONT()` generates a temporal uniform contact network representing
#' homogeneous mixing within the simulated population.
#'
#' In a homogeneous-mixing network, every individual is connected to every
#' other individual with the same edge weight, defined by `val`. Contact
#' matrices are identical across simulation steps, resulting in a fully
#' connected population without temporal or individual heterogeneity in
#' contact structure.
#'
#' This network can be used as a naive transmission scenario in which
#' transmission is not constrained by an underlying contact structure.
#'
#' @param POP A `simPOP` object containing the simulated population.
#'
#' @param NameCont A character string specifying the name assigned to the
#' simulated contact network. If a network with the same name already exists
#' in `POP`, it is replaced with a warning.
#'
#' @param val A positive numeric value specifying the uniform edge weight
#' assigned to every pair of individuals. Default is `1`.
#'
#' @param verbose Logical indicating whether informative messages are printed.
#' Default is `TRUE`.
#'
#' @return A `simPOP` object with the homogeneous-mixing temporal contact
#' matrices added to the `ContactM` slot under `NameCont`.
#'
#' @examples
#' # Create a simple population
#' POP=simplePOP(NIDs = 200, steps = 100, rep = 100)
#'
#' # Simulate a homogeneous-mixing contact network
#' POP=uniCONT(POP=POP, NameCont="Homogeneous", val=1)
#'
#' POP
#'
#' @export
uniCONT=function(POP=NULL, NameCont="Homogeneous", val=1, verbose=TRUE) {

  # Validate simPOP object
  if(!inherits(POP, "simPOP")) stop("`POP` must be class `simPOP`.")
  if(!validate(POP))           stop("Invalid `simPOP` object.")

  # Validate NameCont
  if(!is.character(NameCont) || length(NameCont)!=1L ||
     is.na(NameCont) || !nzchar(NameCont))
    stop("`NameCont` must be a single non-empty character string.")

  if(NameCont %in% names(POP@ContactM))
    warning("Contact network `", NameCont,
            "` already exists in `POP` and has been replaced.")

  # Validate val
  if(!is.numeric(val) || length(val)!=1L ||
     !is.finite(val) || val<=0)
    stop("`val` must be a positive numeric value.")

  # Create homogeneous-mixing contact matrix
  mat=matrix(
    val,
    nrow=length(POP@IDs),
    ncol=length(POP@IDs),
    dimnames=list(POP@IDs, POP@IDs)
  )

  # Remove self-contacts
  diag(mat)=0

  # Repeat the homogeneous network across simulation steps
  Mlist=rep(list(mat), POP@steps-1)

  POP@ContactM[[NameCont]]=Mlist

  # Validate resulting object
  if(!validate(POP))
    stop("The simulation produced an invalid `simPOP` object.")

  return(POP)
}

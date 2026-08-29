#' Simulate temporally correlated contact networks
#'
#' @description
#' `simCONT()` generates temporal undirected contact networks with
#' gamma-distributed edge weights. Contact weights between the same pair of
#' individuals can be correlated across consecutive simulation steps,
#' representing temporal stability in the contact structure.
#'
#' For each pair of individuals, contact weights are generated from a gamma
#' distribution defined by `gshape` and `gscale`. Temporal correlation between
#' contact weights is introduced using an autoregressive correlation structure
#' controlled by `rho`. The resulting contact matrices are symmetric, with
#' zero values on the diagonal.
#'
#' Edge weights can subsequently be rescaled using `val`, and contacts below
#' `min` are set to zero.
#'
#' @param POP A `simPOP` object containing the simulated population.
#'
#' @param NameCont A character string specifying the name assigned to the
#' simulated contact network. If a network with the same name already exists
#' in `POP`, it is replaced with a warning.
#'
#' @param gshape A positive numeric value specifying the shape parameter of
#' the gamma distribution used to generate contact weights. Default is `1`.
#'
#' @param gscale A positive numeric value specifying the scale parameter of
#' the gamma distribution used to generate contact weights. Default is `5`.
#'
#' @param val A non-negative numeric value used as a multiplicative scaling
#' factor for the simulated edge weights. Default is `1`.
#'
#' @param min A non-negative numeric value specifying the minimum edge weight
#' retained in the contact matrices. Values below this threshold are set to
#' zero. Default is `1e-7`.
#'
#' @param rho A numeric value specifying the temporal correlation between
#' contact weights across consecutive simulation steps. Default is `0.25`.
#'
#' @param verbose Logical indicating whether informative messages are printed.
#' Default is `TRUE`.
#'
#' @return A `simPOP` object with the simulated temporal contact matrices added
#' to the `ContactM` slot under `NameCont`.
#'
#' @examples
#' set.seed(123)
#'
#' # Create a simple population
#' POP=simplePOP(NIDs=20, steps=100, rep=100)
#'
#' # Simulate a temporally correlated contact network
#' POP=simCONT(POP=POP,  NameCont="SpatialInteraction",
#'             gshape=1, gscale=5, rho=0.25)
#'
#' POP
#'
#' @export

simCONT=function(POP=NULL, NameCont="Unknown",
                 gshape=1, gscale=5, val=1,
                 min=1e-7, rho=0.25, verbose=TRUE) {

  # Control parameters
  if (!inherits(POP, "simPOP")) stop("Object must be class 'simPOP'.")
  if( !validate(POP))           stop("Invalid `simPOP` object.")

  # Validate NameCont
  if(!is.character(NameCont) || length(NameCont)!=1L ||
     is.na(NameCont) || !nzchar(NameCont))
    stop("`NameCont` must be a single non-empty character string.")

  if(NameCont %in% names(POP@ContactM))
    warning("Contact network `", NameCont,
            "` already exists in `POP` and has been replaced.")

  # Validate gamma parameters
  if(!is.numeric(gshape) || length(gshape)!=1L ||
     !is.finite(gshape)  || gshape<=0)
    stop("`gshape` must be a positive numeric value.")

  if(!is.numeric(gscale) || length(gscale)!=1L ||
     !is.finite(gscale)  || gscale<=0)
    stop("`gscale` must be a positive numeric value.")

  # Validate network scaling
  if(!is.numeric(val) || length(val)!=1L ||
     !is.finite(val)  || val<0)
    stop("`val` must be a non-negative numeric value.")

  if(!is.numeric(min) || length(min)!=1L ||
     !is.finite(min)  || min<0)
    stop("`min` must be a non-negative numeric value.")

  # Validate temporal correlation
  if(!is.numeric(rho) || length(rho)!=1L ||
     !is.finite(rho)  || rho<=-1 || rho>=1)
    stop("`rho` must be a numeric value between -1 and 1.")

  # Validate verbose
  if(!is.logical(verbose) || length(verbose)!=1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.")

  # Simulate temporally correlated gamma-distributed edge weights
  gamma_gcd=simstudy::genCorGen(choose(length(POP@IDs),2),
                          nvars = POP@steps-1, params1 = gshape, params2 = gscale,  rho = rho,
                          dist = "gamma", corstr="ar1", wide = TRUE)
  gamma_gcd=as.matrix(gamma_gcd[,-1,drop=FALSE])

  # Scale edge weights
  gamma_gcd=gamma_gcd*val

  # Create temporal contact matrices
  Mlist=vector("list", POP@steps-1)
  for(i in seq_len(ncol(gamma_gcd))){
    mat=matrix(0, nrow=length(POP@IDs), ncol=length(POP@IDs),
               dimnames=list(POP@IDs,POP@IDs))

    mat[upper.tri(mat)]=gamma_gcd[,i]
    mat[lower.tri(mat)]=t(mat)[lower.tri(mat)]

    # Remove contacts below the minimum threshold
    mat[mat<min]=0

    Mlist[[i]]=mat
  }

  POP@ContactM[[NameCont]]=Mlist

  # Validate simPOP object
  VAL=validate(POP)
  if(!VAL) stop("The simulation produced an invalid `simPOP` object.")

  return(POP)
}

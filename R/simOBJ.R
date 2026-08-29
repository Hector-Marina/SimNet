#' Simulate object-mediated contact networks
#'
#' @description
#' `simOBJ()` generates temporal directed contact networks based on the
#' sequential use of shared objects by individuals in a simulated population.
#'
#' At each simulation step, every individual is assigned `freq` object-use
#' events among `pos` available objects. For individuals using the same object,
#' directed contacts are generated from earlier users to subsequent users.
#' Multiple shared-object encounters between the same pair of individuals are
#' accumulated in the corresponding edge weight.
#'
#' The resulting contact matrices can be used to represent indirect
#' transmission pathways through shared environmental resources, such as
#' milking stations, feeding locations, or resting areas.
#'
#' @param POP A `simPOP` object containing the simulated population.
#'
#' @param NameCont A character string specifying the name assigned to the
#' simulated contact network. If a network with the same name already exists
#' in `POP`, it is replaced with a warning.
#'
#' @param pos A positive integer indicating the number of available objects
#' at each simulation step. Default is `24`.
#'
#' @param freq A positive integer indicating the number of object-use events
#' assigned to each individual per simulation step. Default is `1`.
#'
#' @param val A non-negative numeric value used to multiply the edge weights
#' of the simulated contact matrices. Default is `1`.
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
#' POP=simplePOP(NIDs = 200, steps = 100, rep = 100)
#'
#' # Simulate an object-mediated contact network
#' POP=simOBJ(POP = POP, NameCont = "ObjectNetwork", pos = 24, freq = 1)
#' POP
#'
#' @export

simOBJ=function(POP=NULL, NameCont="Unknown", pos=24, freq=1, val=1, verbose=TRUE) {

  # Control parameters
  if (!inherits(POP, "simPOP")) stop("Object must be class 'simPOP'.")
  if( !validate(POP))           stop("Invalid `simPOP` object.")

  # Validate NameCont
  if (!is.character(NameCont) || length(NameCont) != 1L ||
      is.na(NameCont) || !nzchar(NameCont))
    stop("`NameCont` must be a single non-empty character string.")

  if (NameCont %in% names(POP@ContactM))
    warning("Contact network `", NameCont,
            "` already exists in `POP` and has been replaced.")

  # Validate pos
  if (!is.numeric(pos) || length(pos) != 1L || !is.finite(pos) ||
      pos < 1 || pos %% 1 != 0)
    stop("`pos` must be a positive integer.")

  # Validate freq
  if (!is.numeric(freq) || length(freq) != 1L || !is.finite(freq) ||
      freq < 1 || freq %% 1 != 0)
    stop("`freq` must be a positive integer.")

  # Validate val
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val < 0)
    stop("`val` must be a non-negative numeric value.")

  # Validate verbose
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.")

  # Individuals involved in object-use events
  IDs=rep(POP@IDs, times=freq)

  # Create a direct contact matrix using object time interaction
  Mlist=list()
  for( i in 1:(POP@steps-1)){
    mat=matrix(0, length(POP@IDs), length(POP@IDs))

    # Simulate interaction with the object
    ObjC= cbind(sample(IDs,length(POP@IDs)*freq, replace=F),sample(1:pos,length(POP@IDs)*freq, replace=T))

    for(j in sort(unique(ObjC[,2]))){
      if(sum(ObjC[,2]==j)>1){
        ObjT=ObjC[ObjC[,2]==j,]
        for (k in 1:(nrow(ObjT)-1)){
          mat[ObjT[k,1],ObjT[(k+1):nrow(ObjT),1]]=mat[ObjT[k,1],ObjT[(k+1):nrow(ObjT),1]]+1
        }
      }
    }

    # Remove self-contacts
    diag(mat)=0

    # Matrix weighting
    mat=mat*val

    Mlist[[i]]=mat
  }

  POP@ContactM[[NameCont]]=Mlist

  # Validate simPOP object
  VAL=validate(POP)
  if(!VAL) stop("The simulation produced an invalid `simPOP` object.")

  return(POP)
}

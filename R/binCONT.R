#' Transform contact networks into binary networks
#'
#' @description
#' `binCONT()` transforms a temporal weighted contact network stored in a
#' `simPOP` object into a binary contact network using a specified threshold.
#'
#' For each contact matrix, edge weights below `thr` are set to `0`, whereas
#' edge weights equal to or greater than `thr` are replaced by `val`. The
#' temporal structure and directionality of the original contact network are
#' preserved.
#'
#' If `NameContOut=NULL`, the original contact network is replaced by its
#' binary version. Alternatively, a new network can be created by specifying
#' a different name in `NameContOut`.
#'
#' @param POP A `simPOP` object containing the contact network to transform.
#'
#' @param NameCont A character string specifying the name of the contact
#' network to transform.
#'
#' @param NameContOut A character string specifying the name assigned to the
#' resulting binary contact network. If `NULL`, `NameCont` is overwritten.
#' Default is `NULL`.
#'
#' @param thr A non-negative numeric value specifying the threshold used to
#' classify contacts. Edge weights below `thr` are set to `0`, whereas values
#' equal to or greater than `thr` are replaced by `val`.
#'
#' @param val A positive numeric value assigned to contacts equal to or above
#' the threshold. Default is `1`.
#'
#' @return A `simPOP` object with the binary temporal contact matrices added
#' to the `ContactM` slot under `NameContOut`.
#'
#' @examples
#' set.seed(123)
#'
#' # Create and simulate a contact network
#' POP=simplePOP(NIDs=200, steps=100, rep=100)
#' POP=simCONT(POP=POP, NameCont="SpatialInteraction")
#'
#' # Transform the weighted network into a binary network
#' POP=binCONT(POP=POP, NameCont="SpatialInteraction",
#'             NameContOut="BinaryInteraction", thr=0.10, val=1)
#'
#' POP
#'
#' @export
binCONT=function(POP=NULL, NameCont=NULL, NameContOut=NULL,
                 thr=NULL, val=1) {

  # Validate simPOP object
  if(!inherits(POP, "simPOP"))
    stop("`POP` must be class `simPOP`.")

  if(!validate(POP))
    stop("Invalid `simPOP` object.")

  # Validate NameCont
  if(!is.character(NameCont) || length(NameCont)!=1L ||
     is.na(NameCont) || !nzchar(NameCont))
    stop("`NameCont` must be a single non-empty character string.")

  if(!NameCont %in% names(POP@ContactM))
    stop("Contact network `",NameCont,"` was not found in `POP`.")

  # Validate NameContOut
  if(is.null(NameContOut)){

    NameContOut=NameCont

  }else{

    if(!is.character(NameContOut) || length(NameContOut)!=1L ||
       is.na(NameContOut) || !nzchar(NameContOut))
      stop("`NameContOut` must be a single non-empty character string.")

    if(NameContOut %in% names(POP@ContactM) &&
       NameContOut!=NameCont)
      warning("Contact network `",NameContOut,
              "` already exists in `POP` and has been replaced.")
  }

  # Validate threshold
  if(!is.numeric(thr) || length(thr)!=1L ||
     !is.finite(thr) || thr<0)
    stop("`thr` must be a non-negative numeric value.")

  # Validate contact value
  if(!is.numeric(val) || length(val)!=1L ||
     !is.finite(val) || val<=0)
    stop("`val` must be a positive numeric value.")

  # Transform temporal contact matrices
  Mlist=vector("list",length(POP@ContactM[[NameCont]]))

  for(i in seq_along(POP@ContactM[[NameCont]])){

    mat=POP@ContactM[[NameCont]][[i]]

    mat[mat<thr]=0
    mat[mat>=thr]=val

    # Remove self-contacts
    diag(mat)=0

    Mlist[[i]]=mat
  }

  POP@ContactM[[NameContOut]]=Mlist

  # Validate resulting object
  if(!validate(POP))
    stop("Contact-network transformation produced an invalid `simPOP` object.")

  return(POP)
}

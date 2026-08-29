#' Estimate topological parameters from contact networks
#'
#' @description
#' `parCONT()` estimates topological parameters describing the structure of
#' temporal contact networks stored in a `simPOP` object.
#'
#' For each contact network and simulation step, the function estimates
#' weighted degree, betweenness, and closeness centrality. The mean value of
#' each parameter across individuals is calculated for each temporal network,
#' and subsequently summarised across simulation steps using the mean and
#' interquartile range.
#'
#' @param POP A `simPOP` object containing temporal contact networks.
#'
#' @param NameCont A character vector specifying the contact networks to
#' analyse. If `NULL`, all contact networks stored in `POP` are analysed.
#'
#' @param directed A logical value or vector indicating whether each selected
#' contact network is directed. If `NULL`, network directionality is inferred
#' from the symmetry of the contact matrices.
#'
#' @param verbose Logical indicating whether informative messages are printed.
#' Default is `TRUE`.
#'
#' @return A character matrix containing the mean and interquartile range
#' of weighted degree, betweenness, and closeness for each selected contact
#' network.
#'
#' @examples
#' set.seed(123)
#'
#' # Create a simple population
#' POP=simplePOP(NIDs=200, steps=100, rep=100)
#'
#' # Simulate a temporally correlated contact network
#' POP=simCONT(POP=POP, NameCont="SpatialInteraction",
#'             gshape=1, gscale=5, rho=0.25)
#'
#' # Estimate network topological parameters
#' parCONT(POP=POP, NameCont="SpatialInteraction")
#'
#' @export

parCONT=function(POP=NULL, NameCont=NULL, directed=NULL, verbose=TRUE) {
  # Validate simPOP object
  if(!inherits(POP, "simPOP")) stop("`POP` must be class `simPOP`.")
  if(!validate(POP))           stop("Invalid `simPOP` object.")

  if(length(POP@ContactM)==0) stop("No contact networks found in `POP`.")

  # validate contact networks
  if(is.null(NameCont)) NameCont=names(POP@ContactM)

  if(!is.character(NameCont) || length(NameCont)==0)
    stop("`NameCont` must be a character vector.")

  if(!all(NameCont %in% names(POP@ContactM)))
    stop("Contact network(s) not found in `POP`. Available contact networks: ",
         paste(names(POP@ContactM), collapse=", "),".")

  # Determine network directionality
  if(is.null(directed)){
    directed=vapply(NameCont,function(nm){!all(vapply(
      POP@ContactM[[nm]],isSymmetric,logical(1)))},logical(1))
  }else{
    if(!is.logical(directed) || any(is.na(directed)))
      stop("`directed` must be logical.")
    if(length(directed)==1)
      directed=rep(directed, length(NameCont))
    if(length(directed)!=length(NameCont))
      stop("`directed` must have length 1 or the same length as `NameCont`.")
  }

  # Output table
  tbl=matrix(NA_character_, nrow=3, ncol=length(NameCont),
    dimnames=list(c("Degree","Betweenness","Closeness"),NameCont))

  # Estimate topological parameters
  for(c in seq_along(NameCont)){

    Mlist=POP@ContactM[[NameCont[c]]]

    tmp=matrix(NA_real_, nrow=length(Mlist), ncol=3,
      dimnames=list(NULL,c("Degree","Betweenness","Closeness")))

    for(i in seq_along(Mlist)){

      g=igraph::graph_from_adjacency_matrix(
        adjmatrix=as.matrix(Mlist[[i]]),
        mode=if(directed[c]) "directed" else "undirected",
        weighted=TRUE,
        diag=FALSE)

      # Topological parameters
      Wdegree=igraph::strength(g,
        loops=FALSE,
        mode="all",
        weights=igraph::E(g)$weight)

      Wbetweenness=igraph::betweenness(g,
        weights=igraph::E(g)$weight,
        normalized=FALSE)

      Wcloseness=igraph::closeness(g,
        weights=igraph::E(g)$weight,
        normalized=FALSE)

      tmp[i,]=c(mean(Wdegree, na.rm=TRUE),
                mean(Wbetweenness, na.rm=TRUE),
                mean(Wcloseness, na.rm=TRUE))
    }

    # Summarise across temporal networks
    for(i in seq_len(ncol(tmp))){
      tbl[i,c]=paste0(
        round(mean(tmp[,i], na.rm=TRUE),2),
        " (",
        round(stats::quantile(tmp[,i], probs=0.25, na.rm=TRUE),2),
        "-",
        round(stats::quantile(tmp[,i], probs=0.75, na.rm=TRUE),2),
        ")"
      )
    }
  }

  return(tbl)
}

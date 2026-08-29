#' Infer the role of candidate contact networks in disease transmission
#'
#' @description
#' `RISKcont()` evaluates the probability that observed infection events could
#' have occurred through different temporal contact networks stored in a
#' `simPOP` object.
#'
#' For each newly infected individual, the function identifies the individuals
#' that were infectious at the previous sampling step and extracts their
#' contacts with the newly infected individual during the interval between
#' consecutive sampling steps. The probability of transmission is then
#' calculated independently for each candidate contact network.
#'
#' A transmission probability of zero indicates an impossible transmission
#' event under the corresponding contact network, i.e. the newly infected
#' individual had no effective contact with an infectious individual during
#' the evaluated interval.
#'
#' The resulting event-level transmission probabilities can subsequently be
#' combined across infection events to compare the likelihood of alternative
#' contact networks as potential transmission pathways.
#'
#' @param POP A `simPOP` object containing infection-status information,
#' temporal contact networks, and transmission parameters.
#'
#' @param rep A numeric vector indicating the simulation replicates to analyse.
#' If `NULL`, all replicates are included. Default is `NULL`.
#'
#' @param steps A numeric vector indicating the sampling steps at which
#' infection status is assumed to be observed. Transmission probabilities are
#' evaluated using all contact matrices between consecutive sampling steps.
#' If `NULL`, all simulation steps are used. Default is `NULL`.
#'
#' @param NameCont A character vector specifying the contact networks to
#' evaluate. If `NULL`, all contact networks stored in `POP` are included.
#' Default is `NULL`.
#'
#' @param ref An optional character string specifying a reference contact
#' network used to standardise the likelihood summary printed when
#' `verbose=TRUE`. The reference does not affect the transmission probabilities
#' stored in the returned object. Default is `NULL`.
#'
#' @param verbose Logical indicating whether summaries of transmission
#' probabilities and impossible transmission events are printed.
#' Default is `TRUE`.
#'
#' @return A `simPOP` object with event-level transmission probabilities for
#' each candidate contact network stored in the `RISKcont` slot.
#'
#' @examples
#' \dontrun{
#' POP=RISKcont(POP=POP)
#' }
#'
#' @export

RISKcont=function(POP, rep=NULL, steps=NULL, NameCont=NULL, ref=NULL, verbose=TRUE) {

  # Validate simPOP object
  if(!inherits(POP, "simPOP")) stop("`POP` must be class `simPOP`.")
  if(!validate(POP))           stop("Invalid `simPOP` object.")

  # Validate replicates
  if(is.null(rep)) rep=seq_len(POP@rep)
  if(!is.numeric(rep)     || length(rep)==0L ||
     any(!is.finite(rep)) || any(rep%%1!=0)  ||
     any(!rep %in% seq_len(POP@rep)))
    stop("`rep` contains invalid replicate numbers.")

  # Validate sampling steps
  if(is.null(steps)) steps=seq_len(POP@steps)
  if(!is.numeric(steps)     || length(steps)<2L ||
     any(!is.finite(steps)) || any(steps%%1!=0) ||
     any(!steps %in% seq_len(POP@steps)))
    stop("`steps` must contain at least two valid simulation steps.")
  if(is.unsorted(steps, strictly=TRUE))
    stop("`steps` must be provided in increasing order without duplicates.")

  # Select contact networks
  if(is.null(NameCont)) NameCont=names(POP@ContactM)
  if(!is.character(NameCont) || length(NameCont)<2L ||
     any(is.na(NameCont))    || any(!nzchar(NameCont)))
    stop("`NameCont` must contain at least two contact-network names.")
  if(!all(NameCont %in% names(POP@ContactM)))
    stop("Contact network(s) not found in `POP`. Available contact networks: ",
          paste(names(POP@ContactM), collapse=", "),".")

  # Validate reference network
  if(!is.null(ref)){
    if(!is.character(ref) || length(ref)!=1L ||
       is.na(ref)         || !nzchar(ref))
      stop("`ref` must be a single non-empty character string.")
    if(!ref %in% NameCont)
      stop("`ref` must correspond to one of the contact networks in `NameCont`.")
  }

  # Validate transmission parameter
  if(is.null(POP@TransPar$beta)    || !is.numeric(POP@TransPar$beta) ||
     length(POP@TransPar$beta)!=1L || !is.finite(POP@TransPar$beta))
    stop("A valid transmission parameter `beta` was not found in `POP@TransPar`.")

  # Validate verbose
  if(!is.logical(verbose) || length(verbose)!=1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.")


  # Store transmission probabilities
  TblList=list()
  nTbl=0L

  # Estimate transmission probabilities
  for(r in rep){
    # Infection-status matrix
    S=POP@statusM[[r]]

    for(i in seq.int(2L, length(steps))){
      # Individuals infectious at previous sampling step
      Ainf=sort(which(S[,steps[i-1]]==1))

      # Newly infected individuals
      Ninf=sort(which(paste0(S[,steps[i-1]],S[,steps[i]])=="01"))

      if(length(Ninf)>0){
        for(j in seq_along(Ninf)){
          TProb=base::rep(0, length(NameCont))

          for(k in seq_along(NameCont)){

            x=sapply(steps[i-1]:(steps[i]-1),function(tt){
                POP@ContactM[[NameCont[k]]][[tt]][Ainf,Ninf[j]]})
            x=as.matrix(x)

            if(ncol(x)>1){
              TProb[k]=1-prod(1-(1-exp(-POP@TransPar$beta*x)))
            }else{
              TProb[k]=1-prod(1-(1-exp(-rowSums(POP@TransPar$beta*x))))
            }
          }

          nTbl=nTbl+1L

          TblList[[nTbl]]=data.frame(
            Replicate=r,
            Step=steps[i],
            ContactM=NameCont,
            Susceptible=Ninf[j],
            Ninf=length(Ainf),
            TProb=TProb
          )
        }
      }
    }
  }

  # Combine results
  if(length(TblList)>0){
    Tbl=do.call(rbind,TblList)
    rownames(Tbl)=NULL
  }else{
    Tbl=data.frame(
      Replicate=numeric(),
      Step=numeric(),
      ContactM=character(),
      Susceptible=numeric(),
      Ninf=numeric(),
      TProb=numeric()
    )
  }

  # Store results
  POP@RISKcont=Tbl

  # Print summaries
  if(verbose){
    if(nrow(Tbl)==0){
      message("No new infection events were detected for the selected replicates and steps.")
    }else{
      if(is.null(ref)){
        message("Average transmission probability per contact network:")
        print(tapply(Tbl$TProb,Tbl$ContactM,function(x) mean(x,na.rm=TRUE)))
      }else{
        message(
          "Mean replicate-level likelihood relative to reference contact network `",
          ref,"`:")

        TProbAgg=stats::aggregate(TProb ~ Replicate + ContactM,data=Tbl,FUN=prod)

        L=tapply(TProbAgg$TProb,TProbAgg$ContactM,function(x) mean(x,na.rm=TRUE))

        print(L/L[ref])
      }

      message("Proportion of impossible events per contact network:")
      print(tapply(Tbl$TProb,Tbl$ContactM,function(x) sum(x==0,na.rm=TRUE)
                  )/tapply(Tbl$TProb,Tbl$ContactM,function(x) sum(!is.na(x))))
    }
  }

  # Validate resulting object
  if(!validate(POP))
    stop("Transmission-network inference produced an invalid `simPOP` object.")

  return(POP)
}

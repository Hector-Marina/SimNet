#' Infer potential sources of infection
#'
#' @description
#' `RISKest()` evaluates potential infection sources for newly infected
#' individuals using a selected temporal contact network.
#'
#' For each newly infected individual, the function identifies individuals
#' that were infectious at the previous sampling step and extracts their
#' contacts with the newly infected individual during the interval between
#' consecutive sampling steps. Contact weights are accumulated across the
#' interval and converted into individual transmission risks.
#'
#' Potential infection sources are ranked according to their accumulated
#' contact weights. Shannon entropy and Shannon evenness are calculated from
#' the transmission risks of individuals with non-zero transmission
#' probability to quantify the uncertainty associated with identifying the
#' infection source.
#'
#' Because simulated effective transmission contacts are stored in
#' `POP@TransDF`, the function also compares the inferred potential sources
#' with the effective source generated during the simulation. These results
#' can be used to evaluate the accuracy of infection-source inference.
#'
#' @param POP A `simPOP` object containing infection-status information,
#' temporal contact networks, transmission parameters, and effective
#' transmission contacts.
#'
#' @param rep A numeric vector indicating the simulation replicates to analyse.
#' If `NULL`, all replicates are included. Default is `NULL`.
#'
#' @param steps A numeric vector indicating the sampling steps at which
#' infection status is assumed to be observed. Contacts between consecutive
#' sampling steps are combined when evaluating potential infection sources.
#' If `NULL`, all simulation steps are used. Default is `NULL`.
#'
#' @param NameCont A character string specifying the contact network used to
#' infer potential infection sources. If `NULL`, the contact network recorded
#' in `POP@TransDF` is used when a single network is present. Default is
#' `NULL`.
#'
#' @param thr A numeric value between `0` and `1` specifying the proportional
#' difference from the maximum accumulated contact weight used to identify
#' candidate infection sources. For example, `thr=0.2` retains individuals
#' with contact weights at least 80\% of the maximum. Default is `0.2`.
#'
#' @param verbose Logical indicating whether informative messages are printed.
#' Default is `TRUE`.
#'
#' @return A `simPOP` object with infection-source inference results stored in
#' the `RISKest` slot.
#'
#' @examples
#' \dontrun{
#' POP=RISKest(POP=POP, rep=1)
#' }
#'
#' @export
RISKest=function(POP, rep=NULL, steps=NULL, NameCont=NULL,thr=0.2, verbose=TRUE) {

  # Validate simPOP object
  if(!inherits(POP, "simPOP"))    stop("`POP` must be class `simPOP`.")
  if(!validate(POP))              stop("Invalid `simPOP` object.")

  # Check transmission information
  if(nrow(POP@TransDF)==0)
    stop("No effective transmission contacts found in `POP@TransDF`.")

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

  # Select contact network
  if(is.null(NameCont)){
    TransNetworks=unique(POP@TransDF$ContactM)
    if(length(TransNetworks)==1L){
      NameCont=TransNetworks
    }else{
      stop("`NameCont` must be specified when multiple contact networks are present in `POP@TransDF`.")
    }
  }
  if(!is.character(NameCont) || length(NameCont)!=1L ||
     is.na(NameCont)         || !nzchar(NameCont))
    stop("`NameCont` must be a single non-empty character string.")

  if(!NameCont %in% names(POP@ContactM))
    stop("Contact network `",NameCont,"` was not found in `POP`.")

  if(!NameCont %in% POP@TransDF$ContactM)
    stop("No effective transmission contacts were simulated through `",
         NameCont,"`.")

  # Validate threshold
  if(!is.numeric(thr) || length(thr)!=1L ||
     !is.finite(thr)  || thr<0   ||  thr>1)
    stop("`thr` must be a numeric value between 0 and 1.")

  # Validate verbose
  if(!is.logical(verbose) || length(verbose)!=1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.")

  # Effective transmission contacts
  Tbl=POP@TransDF[POP@TransDF$Replicate %in% rep &
                  POP@TransDF$ContactM==NameCont &
                  POP@TransDF$Step>min(steps)    &
                  POP@TransDF$Step<=max(steps),  ]

  if(nrow(Tbl)==0)
    stop("No effective contacts found with the selected parameters.")

  # Infection-source information
  Tbl$Ninf=NA_real_                         # Number of infected individuals in previous step
  Tbl$indID=NA_character_                   # Individual with the highest contact value
  Tbl$infCont=NA_real_                      # Highest contact value

  Tbl$indIDsL=vector("list",nrow(Tbl))      # List of infected individuals contacted
  Tbl$infContL=vector("list",nrow(Tbl))     # List of contact values for those individuals
  Tbl$infRiskL=vector("list",nrow(Tbl))     # List of transmission risks from those individuals

  Tbl$SEntropy=NA_real_                     # Shannon's entropy value (standarized)
  Tbl$SEvenness=NA_real_                    # Shannon's evenness (Table 2: https://doi.org/10.1111/j.1574-6941.2003.tb01040.x)
  Tbl$NpropPos=NA_real_                     # Number of possible events

  Tbl$PpropPos=NA_real_                     # Position of the proposed causative infected ID
  Tbl$PpropDiff=vector("list", nrow(Tbl))   # Proportion difference with the maximum value
  Tbl$NpropThr=vector("list", nrow(Tbl))    # Number of infected individuals above the threshold
  Tbl$PpropThrIDs=vector("list", nrow(Tbl)) # Names of the infected individuals above the threshold
  Tbl$PpropThrMatch=FALSE                   # Is the causative ID capture by the threshold?

  # Estimate infection-source risk
  for(r in rep){
    # Infection-status matrix
    S=POP@statusM[[r]]

    for(i in seq.int(2L,length(steps))){
      # Infectious individuals at previous sampling step
      Ainf=sort(which(S[,steps[i-1]]==1))
      # Newly infected individuals
      Ninf=sort(which(paste0(S[,steps[i-1]],S[,steps[i]])=="01"))

      if(length(Ninf)>0){
        for(j in seq_along(Ninf)){

          # Contacts between infectious and newly infected individuals
          x=sapply(steps[i-1]:(steps[i]-1),function(tt){
              POP@ContactM[[NameCont]][[tt]][Ainf,Ninf[j]]},simplify="array")


          # Preserve dimensions when only one infectious individual exists
          if(length(x)==(steps[i]-steps[i-1])){
            x=matrix(x,nrow=1)
          }else{
            x=as.matrix(x)
          }
          rownames(x)=Ainf

          # Accumulated contact and transmission risk
          if(ncol(x)==1){
            x_contact=as.vector(x)
            names(x_contact)=rownames(x)
            x_contact=x_contact[order(x_contact,decreasing = T)]

            x_risk=1-exp(-POP@TransPar$beta*x_contact)
          }else{
            x_contact=rowSums(x)
            names(x_contact)=rownames(x)
            x_contact=x_contact[order(x_contact,decreasing = T)]

            x_risk=1-exp(-rowSums(POP@TransPar$beta*x))
          }

          # Find corresponding effective transmission
          k=which(Tbl$Replicate==r & (Tbl$Step>steps[i-1] & Tbl$Step<=steps[i]) & Tbl$Susceptible==Ninf[j])
          if(length(k)==0) next

          # If an individual was infected multiple times between sampling
          # points, only the last infection is evaluated
          Tbl[k,"Ninf"]=0
          if(length(k)>1) k=k[length(k)]

          # Number of infected individuals
          Tbl[k,"Ninf"]=length(Ainf)
          Tbl[k,"indID"]=names(which.max(x_contact))
          Tbl[k,"infCont"]=round(as.numeric(x_contact[which.max(x_contact)]),4)

          # Potential infectious contacts
          Tbl$indIDsL[[k]]=as.vector(names(x_contact))
          Tbl$infContL[[k]]=round(unname(x_contact),4)
          Tbl$infRiskL[[k]]=round(unname(x_risk),4)

          # Shannon entropy and evenness
          p=x_risk[x_risk>0]
          p_se=p/sum(p)

          Tbl[k,"SEntropy"]=-sum(p_se*log2(p_se))
          Tbl[k,"SEvenness"]=-sum(p_se*log2(p_se))/log2(length(p_se))
          Tbl[k,"NpropPos"]=length(p)

          # Compare inferred sources with the simulated effective source
          if(any(names(x_contact)==Tbl$Infected[[k]])){

            Tbl[k,"PpropPos"] =which(Tbl$indIDsL[[k]]==Tbl$Infected[[k]])

            Tbl$PpropDiff[[k]]=x_contact[names(x_contact)==Tbl$Infected[[k]]]/
                               x_contact[which.max(x_contact)]

            candidateIDs=names(x_contact)[
                               x_contact>=x_contact[which.max(x_contact)]*(1-thr)]

            Tbl$NpropThr[[k]]=length(candidateIDs)
            Tbl$PpropThrIDs[[k]]=candidateIDs
            Tbl[k,"PpropThrMatch"]=any(candidateIDs==Tbl$Infected[[k]])
          }
        }
      }
    }
  }

  # Store infection-source inference
  POP@RISKest=Tbl

  # Validate resulting object
  if(!validate(POP))
    stop("Infection-source inference produced an invalid `simPOP` object.")

  return(POP)
}

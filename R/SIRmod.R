#' Simulate stochastic transmission through a contact network
#'
#' @description
#' `SIRmod()` simulates individual-level stochastic disease transmission through
#' a temporal contact network stored in a `simPOP` object.
#'
#' At each simulation step, infectious individuals are evaluated in random
#' order. For each infectious individual, susceptible individuals connected
#' through the selected contact network are exposed to infection according to
#' the corresponding edge weight and the infection parameter `beta`. Once an
#' individual becomes infected during a simulation step, it cannot be infected
#' again by another infectious individual during the same step. Successful
#' transmission events are stored as effective transmission contacts.
#'
#' Infectious individuals recover according to `gamma`. When `epsilon=NULL`,
#' recovered individuals immediately return to the susceptible state, resulting
#' in an SIS model. When `epsilon=0`, recovered individuals remain recovered,
#' resulting in an SIR model. Positive values of `epsilon` allow recovered
#' individuals to return stochastically to the susceptible state, resulting in
#' an SIRS model.
#'
#' The rate parameters `beta`, `gamma`, and `epsilon` are internally converted
#' to transition probabilities using `1-exp(-rate)`.
#'
#' @param POP A `simPOP` object containing the simulated population and at least
#' one temporal contact network.
#'
#' @param NameCont A character string specifying the contact network used to
#' simulate transmission.
#'
#' @param beta A non-negative numeric value specifying the infection rate
#' (`S -> I`). Default is `0.2`.
#'
#' @param gamma A non-negative numeric value specifying the recovery rate
#' (`I -> R`). Default is `0.2`.
#'
#' @param epsilon A non-negative numeric value specifying the rate at which
#' recovered individuals return to the susceptible state (`R -> S`).
#' Set to `NULL` to simulate an SIS model. Default is `0.2`.
#'
#' @param IDseed A numeric vector containing the identifiers of individuals
#' assigned as infectious at the first simulation step. If `NULL`,
#' approximately 1\% of the population is randomly selected.
#'
#' @param verbose Logical indicating whether informative messages are printed.
#' Default is `TRUE`.
#'
#' @return A `simPOP` object containing the simulated infection-status matrices
#' in `statusM`, the effective transmission contacts in `TransDF`, and the
#' transformed transmission probabilities in `TransPar`.
#'
#' @examples
#' set.seed(123)
#'
#' # Create a simple population
#' POP=simplePOP(NIDs=200, steps=100, rep=100)
#'
#' # Simulate a temporal contact network
#' POP=simCONT(POP=POP, NameCont="SpatialInteraction",
#'             gshape=0.01, gscale=0.05, rho=0.25)
#'
#' # Simulate transmission using an SIS model
#' POP=SIRmod(POP=POP, NameCont="SpatialInteraction",
#'            beta=0.015, gamma=0.015, epsilon=NULL, IDseed=c(1:7))
#' POP
#'
#' @export
#'
SIRmod=function(POP=NULL, NameCont=NULL,
                beta=0.2, gamma=0.2, epsilon=0.2,
                IDseed=NULL, verbose=TRUE) {


  # Validate simPOP object
  if(!inherits(POP, "simPOP")) stop("`POP` must be class `simPOP`.")
  if(!validate(POP))           stop("Invalid `simPOP` object.")

  # Validate contact network
  if(!is.character(NameCont) || length(NameCont)!=1L ||
     is.na(NameCont) || !nzchar(NameCont))
    stop("`NameCont` must be a single non-empty character string.")

  if(!NameCont %in% names(POP@ContactM))
    stop("Contact network `", NameCont, "` not found in `POP`.")

  # Validate transmission parameters
  if(!is.numeric(beta) || length(beta)!=1L ||
     !is.finite(beta)  || beta<0   || beta>1)
    stop("`beta` must be a non-negative numeric value between 0 and 1.")

  if(!is.numeric(gamma) || length(gamma)!=1L ||
     !is.finite(gamma)  || gamma<0  || gamma>1)
    stop("`gamma` must be a non-negative numeric value between 0 and 1.")

  if(!is.null(epsilon) &&
     (!is.numeric(epsilon) || length(epsilon)!=1L   ||
      !is.finite(epsilon)  || epsilon<0 || epsilon>1))
    stop("`epsilon` must be a non-negative numeric value or NULL.")

  # Initial infected individuals
  if(is.null(IDseed))
    IDseed=sample(POP@IDs, ceiling(length(POP@IDs)/100))

  if(!is.numeric(IDseed) || length(IDseed)==0L ||
     any(!is.finite(IDseed)))
    stop("`IDseed` must contain valid individual identifiers.")

  if(anyDuplicated(IDseed))
    stop("`IDseed` must contain unique individual identifiers.")

  if(!all(IDseed %in% POP@IDs))
    stop("All `IDseed` values must correspond to individuals in `POP`.")

  # Validate verbose
  if(!is.logical(verbose) || length(verbose)!=1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.")


  # Convert transmission rates to probabilities
  p.Beta=1-exp(-beta)
  p.gamma=1-exp(-gamma)

  if(is.null(epsilon)){
    p.epsilon=NULL
  }else{
    p.epsilon=1-exp(-epsilon)
  }


  # Store effective transmission contacts
  TransList=list()
  nTrans=0L

  # Simulate transmission across replicates
  for(rep in seq_len(POP@rep)){

    Status=POP@statusM[[rep]]

    # Reset previous infection information
    if(any(Status!=0)){
      Status[]=0

      if(verbose && rep==1)
        message("Infection information found in `statusM` and reset before simulation.")
    }

    # Assign initially infected individuals
    Status[IDseed,1]=1
      # being 0=Susceptible; 1=Infectious; 2=Recovered

    # Simulate transmission across steps
    for(step in 2:POP@steps){

      # Contact network corresponding to the transition
      net=POP@ContactM[[NameCont]][[step-1]]

      # Individuals infected at the previous step
      InfIDs=which(Status[,step-1]==1)

      if(length(InfIDs)>0){

        # Randomise order of infectious individuals
        InfIDs=InfIDs[sample(seq_along(InfIDs), size=length(InfIDs),
                 replace=FALSE)]

        for(id in InfIDs){
          # Individuals susceptible at the previous and current step
          SIds=which(Status[,step-1]==0 & Status[,step]==0)

          if(length(SIds)>0){

            # Probability of transmission
            trans_prob=net[id,SIds]*p.Beta
            trans_prob[trans_prob>1]=1

            # Sample transmission events
            RInf=SIds[stats::rbinom(length(trans_prob),size=1,
                                    prob=trans_prob)==1]
            if(length(RInf)>0){
              Status[RInf,step]=1

              nTrans=nTrans+1L
              TransList[[nTrans]]=data.frame(
                Replicate=rep,
                Step=step,
                Infected=id,
                Susceptible=RInf,
                ContactM=NameCont
              )
            }
          }
        }
      }

      # Infected to recovered
      Status=SIRstat(Status=Status, step=step, stat1=1, stat2=2, par=p.gamma)


      # Recovered to susceptible
      if(is.null(epsilon)){
        Status[Status[,step]==2,step]=0
      }else{
        Status=SIRstat(Status=Status, step=step, stat1=2, stat2=0, par=p.epsilon)
      }
    }
    POP@statusM[[rep]]=Status
  }

  # Effective transmission contacts
  if(length(TransList)>0){
    TransDF=do.call(rbind, TransList)
    rownames(TransDF)=NULL
  }else{
    TransDF=data.frame(
      Replicate=numeric(),
      Step=numeric(),
      Infected=numeric(),
      Susceptible=numeric(),
      ContactM=character()
    )
  }

  POP@TransDF=TransDF

  # Transmission probabilities
  POP@TransPar=list(
    beta=p.Beta,
    gamma=p.gamma,
    epsilon=p.epsilon
  )

  # Validate resulting object
  if(!validate(POP))
    stop("Transmission simulation produced an invalid `simPOP` object.")

  return(POP)
}

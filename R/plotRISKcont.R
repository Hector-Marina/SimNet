#' Plot transmission probabilities across candidate contact networks
#'
#' @description
#' `plotRiskcont()` visualises and summarises the transmission probabilities
#' estimated for candidate contact networks using `RISKcont()`.
#'
#' The function produces a horizontal boxplot showing the distribution of
#' transmission probabilities (`TProb`) across observed infection events for
#' each candidate contact network.
#'
#' In addition, transmission probabilities are multiplied across infection
#' events within each replicate and contact network to obtain replicate-level
#' likelihoods. These likelihood distributions are summarised by their mean,
#' median, and 2.5th and 97.5th percentiles.
#'
#' The function also estimates the pairwise overlap between likelihood
#' distributions. Overlap is calculated as the shared area under the estimated
#' probability-density curves after restricting each distribution to its
#' 2.5th–97.5th percentile interval. Larger values indicate greater similarity
#' between the likelihood distributions of two candidate contact networks.
#'
#' @param x A `simPOP` object containing contact-network inference results in
#' the `RISKcont` slot, generated using `RISKcont()`.
#'
#' @param ... Additional graphical parameters passed to `boxplot()`.
#'
#' @return The function is called primarily for its graphical side effect.
#' Descriptive statistics of the replicate-level likelihood distributions and
#' their pairwise overlap are printed to the console.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' plotRiskcont(POP)
#' }
#'
#' @export
plotRiskcont=function(x,...) {
    # Validate simPOP object
    if(!inherits(x, "simPOP")) stop("`x` must be class `simPOP`.")
    if(!validate(x))           stop("Invalid `simPOP` object.")

    if(nrow(x@RISKcont)==0)
      stop("No contact-network inference results found in `x@RISKcont`.")

    if(!all(c("Replicate","ContactM","TProb") %in% names(x@RISKcont)))
      stop("`x@RISKcont` does not contain the required information.")

    # Preserve graphical parameters
    oldpar=graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(oldpar), add=TRUE)

    # Transmission probabilities per infection event
    TProb=x@RISKcont

    area_groups=split(TProb$TProb, TProb$ContactM)
    colors=grDevices::rainbow(length(area_groups))

    # Plot transmission probabilities
    graphics::par(bg="white", mfrow=c(1,1), mar=c(4,8,2,1))

    graphics::boxplot(
      TProb ~ ContactM,
      data=TProb,
      horizontal=TRUE,
      col=grDevices::adjustcolor(colors, alpha.f=0.3),
      border=colors,
      xlab="Transmission probability",
      ylab="",
      las=1,
      ...
    )

    # Calculate replicate-level likelihoods
    TProbAgg=stats::aggregate(
      TProb ~ Replicate + ContactM,
      data=TProb,
      FUN=prod
    )

    ContactM_groups=split(TProbAgg$TProb, TProbAgg$ContactM)

    # Descriptive statistics
    Q=matrix(
      NA_real_,
      nrow=4,
      ncol=length(ContactM_groups),
      dimnames=list(
        c("mean","median","2.5p","97.5p"),
        names(ContactM_groups)
      )
    )

    # Pairwise overlap between likelihood distributions
    O=matrix(
      0,
      nrow=length(ContactM_groups),
      ncol=length(ContactM_groups),
      dimnames=list(
        names(ContactM_groups),
        names(ContactM_groups)
      )
    )

    for(i in seq_along(ContactM_groups)){

      Q[,i]=c(
        mean(ContactM_groups[[i]], na.rm=TRUE),
        stats::median(ContactM_groups[[i]], na.rm=TRUE),
        stats::quantile(
          ContactM_groups[[i]],
          probs=c(0.025,0.975),
          na.rm=TRUE
        )
      )
    }

    # Estimate overlap between distributions
    for(i in seq_along(ContactM_groups)){
      for(j in seq_along(ContactM_groups)){

        if(i!=j){

          xi=ContactM_groups[[i]]
          xi=xi[xi>=Q[3,i] & xi<=Q[4,i]]
          xi=xi[!is.na(xi)]

          yj=ContactM_groups[[j]]
          yj=yj[yj>=Q[3,j] & yj<=Q[4,j]]
          yj=yj[!is.na(yj)]

          # Estimate probability-density curves
          dx=stats::density(xi)
          dy=stats::density(yj)

          # Interpolate second density to the first density grid
          dy_interp=stats::approx(
            dy$x,
            dy$y,
            xout=dx$x
          )$y

          # Shared density
          min_density=pmin(dx$y, dy_interp)

          # Approximate shared area under the curves
          grid_step=diff(dx$x[1:2])

          O[i,j]=sum(min_density, na.rm=TRUE)*grid_step
        }
      }
    }

    message("Descriptive statistics:")
    print(Q)

    message("Shared ContactM under the curves:")
    print(O)

    invisible(list(
      statistics=Q,
      overlap=O
    ))
  }

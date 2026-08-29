#' Plot relative transmission probabilities across contact networks
#'
#' @description
#' `plotRISKcont2()` compares the transmission probabilities estimated for
#' candidate contact networks relative to a reference network.
#'
#' For each observed infection event, the transmission probability estimated
#' for each candidate contact network is divided by the probability obtained
#' from the reference network. Values greater than `1` indicate that the
#' candidate network provides a higher transmission probability than the
#' reference network, whereas values below `1` indicate a lower probability.
#'
#' @param x A `simPOP` object containing contact-network inference results in
#' the `RISKcont` slot, generated using `RISKcont()`.
#'
#' @param ref A character string specifying the contact network used as
#' reference. Default is last.
#'
#' @param Nsim An optional character string used as the plot title.
#' Default is `NULL`.
#'
#' @param colors A named character vector specifying the colours used for the
#' candidate contact networks. Names must correspond to the contact-network
#' names displayed in the plot.
#'
#' @return A `ggplot` object showing the distribution of transmission
#' probabilities relative to the selected reference network.
#'
#' @examples
#' \dontrun{
#' colors=c("SIN"="magenta3", "ASN"="red2",
#'          "MPN"="dodgerblue", "RCN"="green3")
#'
#' plotRISKcont2(POP, ref="Uniform", colors=colors)
#' }
#' @export
plotRISKcont2=function(x, ref=NULL, Nsim=NULL, colors=c(
      "SIN"="magenta3","ASN"="red2","MPN"="dodgerblue","RCN"="green3")){

  # Validate simPOP object
  if(!inherits(x, "simPOP")) stop("`x` must be class `simPOP`.")
  if(!validate(x))           stop("Invalid `simPOP` object.")

  if(nrow(x@RISKcont)==0)
    stop("No contact-network inference results found in `x@RISKcont`.")

  # Validate reference network
  if(is.null(ref))
    ref=names(x@ContactM)[length(x@ContactM)]

  # Validate reference network
  if(!is.character(ref) || length(ref)!=1L ||
     is.na(ref) || !nzchar(ref))
    stop("`ref` must be a single non-empty character string.")

  if(!ref %in% x@RISKcont$ContactM)
    stop("Reference contact network `",ref,"` was not found in `x@RISKcont`.")

  # Transmission probabilities
  TProb=x@RISKcont

  # Extract reference probabilities per infection event
  ref_vals=TProb[TProb$ContactM==ref,c("Replicate","Step","Susceptible","TProb")]
  names(ref_vals)[4]="TProb_ref"

  # Match candidate and reference probabilities
  TProb2=merge(TProb, ref_vals, by=c("Replicate","Step","Susceptible"))

  # Relative transmission probability
  TProb2$Prop=TProb2$TProb/TProb2$TProb_ref

  # Remove reference network
  TProb2=TProb2[TProb2$ContactM!=ref,]

  # Preserve colour-defined plotting order
  if(!is.null(names(colors))){
    levels_plot=c(
      names(colors)[names(colors) %in% unique(TProb2$ContactM)],
      setdiff(unique(TProb2$ContactM),names(colors))
    )

    TProb2$ContactM=factor(
      TProb2$ContactM,
      levels=levels_plot
    )
  }

  # Validate colours
  if(!is.null(colors)){
    missing_colors=setdiff(
      unique(as.character(TProb2$ContactM)),
      names(colors)
    )

    if(length(missing_colors)>0)
      stop(
        "Colours are missing for the following contact networks: ",
        paste(missing_colors,collapse=", "),
        "."
      )
  }

  # Plot
  ggplot2::ggplot(TProb2,
                  ggplot2::aes(x=ContactM, y=Prop,
                               colour=ContactM, fill=ContactM)) +
    ggplot2::geom_boxplot(width=0.5, alpha=0.25, outlier.shape=NA) +
    ggplot2::geom_hline( yintercept=1, colour="grey40",
                         linetype="dashed", linewidth=1) +
    ggplot2::coord_cartesian(ylim=c(0,3)) +
    ggplot2::scale_colour_manual(values=colors) +
    ggplot2::scale_fill_manual(values=colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none",
                   plot.title=ggplot2::element_text(hjust=0.5)) +
    ggplot2::labs( title=Nsim, x=NULL, y=NULL)
}

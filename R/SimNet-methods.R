#' Show method for simPOP
#'
#' @description
#' Displays a summary of the information contained in an `simPOP` object.
#'
#' @param object An `simPOP` object.
#'
#' @importMethodsFrom methods show
#'
#' @export
#'
setMethod("show", "simPOP", function(object) {
  cat("-------simulated POPulation object-------\n")
  cat("Number of individuals:", length(object@IDs), "\n")
  cat("Individuals information:", nrow(object@IDsinfo)!=0, "\n")
  cat("Number of steps:", object@steps, "\n")
  cat("Number of replicates:", object@rep, "\n")
  if(length(object@ContactM)==0){
    cat("Contact matrices:", length(object@ContactM)!=0, "\n")
  }else{
    cat("Contact matrices:", c(names(object@ContactM)), "\n")
  }
  cat("Transmission information:", nrow(object@TransDF)!=0, "\n")
  cat("Transmission pathway information:", nrow(object@RISKest)!=0, "\n")
})



#' Plot simPOP objects
#'
#' @description
#' `plot()` provides graphical summaries of the population, contact networks,
#' and disease-transmission dynamics stored in a `simPOP` object.
#'
#' Four plot types are available:
#'
#' * `"inf.mat"`: infection-status matrix showing the status of each individual
#'   across simulation steps.
#' * `"inf.curv"`: temporal trajectories of the number of susceptible,
#'   infected, and recovered individuals. For multiple replicates, the mean
#'   trajectory and its variability are displayed.
#' * `"cont.net"`: aggregated contact network obtained by summing the selected
#'   temporal contact matrices.
#' * `"inf.net"`: realised transmission network containing the effective
#'   transmission contacts recorded during the stochastic simulation.
#'
#' @param x A `simPOP` object containing population, contact-network, and
#' transmission information.
#'
#' @param rep Numeric scalar or vector indicating the simulation replicate(s)
#' to plot. For `"inf.mat"` and `"inf.net"`, a single replicate is used.
#' For `"inf.curv"`, multiple replicates can be supplied to summarise
#' transmission dynamics across simulations. If `NULL`, the default depends
#' on the selected plot type.
#'
#' @param steps Numeric vector indicating the simulation steps to include.
#' If `NULL`, all available simulation steps are used.
#'
#' @param ContactM Character string indicating the name of the contact network to
#' plot. Required for `type = "cont.net"`. For `type = "inf.net"`, the first
#' contact network stored in the object is used when `ContactM = NULL`.
#'
#' @param type Character string indicating the plot to produce. Available
#' options are `"inf.mat"`, `"inf.curv"`, `"cont.net"`, and `"inf.net"`.
#' The default is `"inf.mat"`
#'
#' @param ... Additional graphical parameters passed to the corresponding
#' plotting function where applicable.
#'
#' @return The function is called for its graphical side effect and produces
#' the selected plot.
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
#'
#' # Plot infection-status matrix for the first replicate
#' plot(POP, rep = 1, type = "inf.mat")
#'
#' # Plot infection trajectories across replicates
#' plot(POP, type = "inf.curv")
#'
#' @export
methods::setMethod("plot", "simPOP",
                   function(x, rep=NULL, steps=NULL, ContactM=NULL,
                            type="inf.mat",...) {

  # Control parameters
  if (!inherits(x, "simPOP")) stop("Object must be class 'simPOP'.")
  if( !validate(x))           stop("Invalid `simPOP` object.")

  type=match.arg(type,
                 choices = c("inf.mat", "inf.curv", "cont.net", "inf.net"))
  if(is.null(steps)) steps=1:x@steps

  if(type=="inf.mat"){
    if(is.null(rep)) rep=1
    Status=x@statusM[[rep]]
    Status=Status[,steps]

    cols=c(grDevices::rgb(0, 1, 0, alpha=0.4),"orange","blue")
    if(length(unique(as.vector(Status)))==2) cols=cols[1:2]
    if(length(unique(as.vector(Status)))==1) cols=cols[1]

    graphics::image(t(apply(Status, 2, rev)), col = cols, axes = FALSE,main = "Matrix Plot")
    graphics::axis(1, at=seq(0, 1, length.out=ncol(Status)), labels=1:ncol(Status))
    graphics::axis(2, at=seq(0, 1, length.out=nrow(Status)), labels=nrow(Status):1, las=1)
  }

  if(type=="inf.curv"){
    if(is.null(rep)) rep=1:x@rep


    if (length(rep)==1){
      Status=x@statusM[[rep]]
      Status=Status[,steps]

      plot(colSums(Status == 0), type = "l", col = "darkgreen", ylim = c(0, length(x@IDs)), xlab = "Time step", ylab = "Count", ...)
      graphics::lines(colSums(Status == 1), col = "orange")
      graphics::lines(colSums(Status == 2), col = "blue")
      graphics::legend("right", legend = c("Susceptible", "Infected", "Recovered"), col = c("darkgreen", "orange", "blue"), lty = 1, bty = "n")
    }else{
      Status=lapply(x@statusM[rep], function(M) M[, steps, drop = FALSE])

      Status_stats=function(val) {counts=sapply(Status, function(M) {colSums(M == val)})
        list(min = apply(counts, 1, min),        mean = apply(counts, 1,  mean),
             sd  =apply(counts, 1, stats::sd),   max  = apply(counts, 1, max))}
      stats0=Status_stats(0)
      stats1=Status_stats(1)
      stats2=Status_stats(2)
      statsL=list(stats0,stats1,stats2)

      col=c("darkgreen","orange","blue")
      for (i in 1:3){
        if(i==1){
          plot(statsL[[i]]$mean,type="l", col=col[i], ylim = c(0,length(x@IDs)),
               xlab = "Time step", ylab = "Count")#, ...)
        }else{
          graphics::lines(statsL[[i]]$mean, col=col[i])
        }
        x_axis=c(steps, rev(steps))
        y_axis=c(statsL[[i]]$mean-statsL[[i]]$sd, rev(statsL[[i]]$mean+statsL[[i]]$sd))
        polygon(x_axis, y_axis, col = grDevices::adjustcolor(col[i], alpha.f = 0.1), border = NA)
        rm(x_axis,y_axis)
      }
      graphics::legend("right", legend = c("Susceptible", "Infected", "Recovered"), col = c("darkgreen", "orange", "blue"), lty = 1, bty = "n")
    }
  }

  if(type=="cont.net"){
    m=Reduce(`+`, x@ContactM[[ContactM]][steps])
    colnames(m)=rownames(m)=1:nrow(m)

    edge_list=as.data.frame(as.table(m))
    edge_list=subset(edge_list, Var1 != Var2 & Freq != 0)
    names(edge_list)=c("from", "to", "weight")

    graph=igraph::graph_from_data_frame(edge_list, directed = TRUE)

    igraph::V(graph)$name  =igraph::V(graph)$name  # ensure node names match 1:5
    igraph::E(graph)$weight=edge_list$weight


    plot(graph,
         edge.arrow.size = 0.1,
         vertex.size = 8,
         vertex.label.cex = 0.8,
         edge.width = igraph::E(graph)$weight / max(igraph::E(graph)$weight) * 5,  # scale width 1–5
         main = paste0("Transmission Network (ContactM:",ContactM,";steps:",length(steps),")"),...)
  }

  if(type=="inf.net"){
    if(is.null(rep)) rep=1

    TransDF=x@TransDF
    TransDF=TransDF[TransDF$Replicate==rep,]

    if(!is.null(steps))TransDF=TransDF[TransDF$Step %in% steps,]
    if(is.null(ContactM))  ContactM=names(x@ContactM)[1]
    TransDF=TransDF[TransDF$ContactM %in% ContactM,]

    graph=igraph::graph_from_edgelist(as.matrix(TransDF[,c(3,4)]), directed = TRUE)

    targets=setdiff(unique(TransDF[,4]), unique(TransDF[,3]))
    both= intersect(unique(TransDF[,4]), unique(TransDF[,3]))

    igraph::V(graph)$name=igraph::V(graph)$name  # ensure node names match 1:5
    igraph::V(graph)$color=ifelse(igraph::V(graph)$name %in% targets, "darkgreen",
                           ifelse(igraph::V(graph)$name %in% both   , "darkorange2","red3"))

    # Plot
    plot(graph,
         vertex.color = igraph::V(graph)$color,
         edge.width=0.1,
         edge.arrow.size=0.1,
         main = paste0("Transmission Network (rep:",rep,";ContactM:",ContactM,";steps:",steps,")"))
    message("Plot legend: Red nodes represent original infections, orange nodes represent susceptible animals that become infected and infect another, and green nodes represent susceptible individuals that become infected but do not infect anyone else.")
  }

})


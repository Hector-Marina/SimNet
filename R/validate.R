#' Validate simPOP object's information
#'
#' @description
#' `validate()` performs a series of structural and consistency checks on a
#' `simPOP` object. The function verifies population information, infection
#' status matrices, temporal contact networks, transmission parameters, and
#' transmission-inference results when available.
#'
#' @param POP A `simPOP` object.
#'
#' @return A logical value. `TRUE` indicates that no inconsistencies were
#' detected and `FALSE` indicates that one or more validation checks failed.
#'
#' @keywords internal
#'
#' @examples
#' POP=simplePOP(NIDs = 200, steps = 100, rep = 100)
#' validate(POP)
#'
#' @export
validate=function(POP){
  VAL=TRUE

  #-- Check object class
  if (!methods::is(POP, "simPOP")) stop("`POP` must be a `simPOP` object.")

  #-- Method
  if (!is.character(POP@method) || length(POP@method) != 1L ||
      is.na(POP@method) || !nzchar(POP@method)) {
    message("POP validation: Population method information is incorrect.")
    VAL=FALSE
  }

  #-- IDs
  if (!is.numeric(POP@IDs)     || length(POP@IDs) == 0L  ||
      any(!is.finite(POP@IDs)) || any(POP@IDs %% 1 != 0) ||
      any(POP@IDs <= 0)        || anyDuplicated(POP@IDs)) {
    message("POP validation: Individual identifiers are incorrect.")
    VAL=FALSE
  }

  #-- IDsinfo
  if (!is.data.frame(POP@IDsinfo)) {
    message("POP validation: Individual information must be a data frame.")
    VAL=FALSE
  }

  if (is.data.frame(POP@IDsinfo) && nrow(POP@IDsinfo) > 0L &&
      nrow(POP@IDsinfo) != length(POP@IDs)) {
    message("POP validation: Individual identifiers are inconsistent with ",
            "individual-level information.")
    VAL=FALSE
  }

  #-- steps
  if (!is.numeric(POP@steps) || length(POP@steps) != 1L ||
      !is.finite(POP@steps)  || POP@steps < 1    || POP@steps %% 1 != 0) {
    message("POP validation: Number of simulation steps is incorrect.")
    VAL=FALSE
  }

  #-- rep
  if (!is.numeric(POP@rep) || length(POP@rep) != 1L ||
      !is.finite(POP@rep)  || POP@rep < 1    || POP@rep %% 1 != 0) {
    message("POP validation: Number of simulation replicates is incorrect.")
    VAL=FALSE
  }

  #-- statusM
  if (!is.list(POP@statusM)) {
    message("POP validation: Infection-status information must be a list.")
    VAL=FALSE
  }

  if (is.list(POP@statusM) && length(POP@statusM) != POP@rep) {
    message("POP validation: Number of status matrices does not correspond ",
            "to the number of replicates.")
    VAL=FALSE
  } else {
    ok=vapply(POP@statusM, function(M) {
          is.matrix(M)               && is.numeric(M)        &&
          nrow(M) == length(POP@IDs) && ncol(M) == POP@steps &&
          all(is.finite(M))          && all(M %in% c(0, 1, 2))
      },logical(1))

    if (!all(ok)) {
      message("POP validation: Infection-status matrices are incorrect. ",
              "Matrices must contain one row per individual, one column per ",
              "simulation step, and status values 0, 1, or 2.")
      VAL=FALSE
    }
  }

  #-- ContactM
  if (!is.list(POP@ContactM)) {
    message("POP validation: Contact-network information must be a list.")
    VAL=FALSE
  }

  if (is.list(POP@ContactM) && length(POP@ContactM) > 0L) {

    # Network names
    if (is.null(names(POP@ContactM)) || any(names(POP@ContactM) == "") ||
        anyDuplicated(names(POP@ContactM))) {
      message(
        "POP validation: Contact networks must have unique, non-empty names."
      )
      VAL=FALSE
    }

    for (ii in seq_along(POP@ContactM)) {
      net=POP@ContactM[[ii]]
      net_name=names(POP@ContactM)[ii]
      if (!is.list(net)) {
        message("POP validation: Contact network `", net_name,
          "` must contain a list of matrices.")
        VAL=FALSE
        next
      }

      if (!length(net) == (POP@steps - 1L)) {
        message("POP validation: Contact network `", net_name,
                "` contains an incorrect number of temporal matrices.")
        VAL=FALSE
      }

      if (length(net) > 0L) {
        ok=vapply(net, function(M) {is.matrix(M) && is.numeric(M) &&
              nrow(M) == length(POP@IDs) &&
              ncol(M) == length(POP@IDs) &&
              all(is.finite(M)) && all(M >= 0)},logical(1))

        if (!all(ok)) {
          message("POP validation: Contact network `", net_name,
                  "` contains invalid contact matrices.")
          VAL=FALSE
        }

        # Self-contact should not occur
        ok_diag=vapply(net,function(M) {
          is.matrix(M) && nrow(M) == ncol(M) && all(diag(M) == 0)
          },logical(1))

        if (!all(ok_diag)) {
          message("POP validation: Contact network `", net_name,
                  "` contains non-zero diagonal values.")
          VAL=FALSE
        }
      }
    }
  }


  #-- TransPar
  if (!is.list(POP@TransPar)) {
    message("POP validation: Transmission parameters must be stored as a list.")
    VAL=FALSE
  }

  if (is.list(POP@TransPar) && length(POP@TransPar) > 0L) {

    if (!all(c("beta", "gamma", "epsilon") %in% names(POP@TransPar))) {
      message(
        "POP validation: Transmission parameters must contain `beta`, ",
        "`gamma`, and `epsilon`."
      )
      VAL=FALSE
    } else {

      if (!is.numeric(POP@TransPar$beta) ||
          length(POP@TransPar$beta) != 1L ||
          !is.finite(POP@TransPar$beta) ||
          POP@TransPar$beta < 0 ||
          POP@TransPar$beta > 1) {
        message("POP validation: Transmission parameter `beta` is incorrect.")
        VAL=FALSE
      }

      if (!is.numeric(POP@TransPar$gamma) ||
          length(POP@TransPar$gamma) != 1L ||
          !is.finite(POP@TransPar$gamma) ||
          POP@TransPar$gamma < 0 ||
          POP@TransPar$gamma > 1) {
        message("POP validation: Transmission parameter `gamma` is incorrect.")
        VAL=FALSE
      }

      if (!is.null(POP@TransPar$epsilon) &&
          (!is.numeric(POP@TransPar$epsilon) ||
           length(POP@TransPar$epsilon) != 1L ||
           !is.finite(POP@TransPar$epsilon) ||
           POP@TransPar$epsilon < 0 ||
           POP@TransPar$epsilon > 1)) {
        message("POP validation: Transmission parameter `epsilon` is incorrect.")
        VAL=FALSE
      }
    }
  }

  #-- TransDF
  if (!is.data.frame(POP@TransDF)) {
    message("POP validation: Transmission information must be a data frame.")
    VAL=FALSE
  }

  if (is.data.frame(POP@TransDF) && nrow(POP@TransDF) > 0L) {
    req=c("Replicate", "Step", "Infected", "Susceptible", "ContactM")

    if (!all(req %in% names(POP@TransDF))) {
      message("POP validation: Transmission data are missing required columns.")
      VAL=FALSE
    } else {

      if (any(!POP@TransDF$Replicate %in% seq_len(POP@rep))) {
        message(
          "POP validation: Transmission data contain invalid replicate numbers."
        )
        VAL=FALSE
      }

      if (any(!POP@TransDF$Step %in% seq_len(POP@steps))) {
        message(
          "POP validation: Transmission data contain invalid simulation steps."
        )
        VAL=FALSE
      }

      if (any(!POP@TransDF$Infected    %in% POP@IDs) ||
          any(!POP@TransDF$Susceptible %in% POP@IDs)) {
        message(
          "POP validation: Transmission data contain unknown individual IDs."
        )
        VAL=FALSE
      }

      if (length(POP@ContactM) > 0L &&
          any(!POP@TransDF$Area %in% names(POP@ContactM))) {
        message(
          "POP validation: Transmission data contain unknown contact networks."
        )
        VAL=FALSE
      }
    }
  }


  #-- RISKest
  if (!is.data.frame(POP@RISKest)) {
    message(
      "POP validation: Infection-source inference information must be a data frame."
    )
    VAL=FALSE
  }

  if (is.data.frame(POP@RISKest) && nrow(POP@RISKest) > 0L) {
    req=c("Replicate","Step","Infected","Susceptible","ContactM","Ninf","SEntropy",
          "SEvenness","NpropPos")

    if (!all(req %in% names(POP@RISKest))) {
      message(
        "POP validation: Infection-source inference results are missing ",
        "required columns."
      )
      VAL=FALSE
    }
  }

  #-- RISKcont
  if (!is.data.frame(POP@RISKcont)) {
    message(
      "POP validation: Contact-network inference information must be a data frame."
    )
    VAL=FALSE
  }

  if (is.data.frame(POP@RISKcont) && nrow(POP@RISKcont) > 0L) {
    req=c("Replicate","Step","ContactM","Susceptible","Ninf","TProb")

    if (!all(req %in% names(POP@RISKcont))) {
      message(
        "POP validation: Contact-network inference results are missing ",
        "required columns."
      )
      VAL=FALSE
    }
  }

  return(VAL)
}

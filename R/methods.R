#' @title Methods for flipscores objects
#'
#' @description Methods for \code{flipscores} objects.
#' The following are methods to extract and manipulate relevant information from
#' a \code{flipscores} object.
#'
#' @name flipscores-method
#' @docType methods

NULL


#' Update method for flipscores
#'
#' Ensures that \code{update()} on a \code{flipscores} object returns
#' a proper \code{flipscores} object rather than a plain \code{glm}.
#'
#' @param object A \code{flipscores} object.
#' @param ... Additional arguments to update in the call, such as
#'   \code{formula}, \code{data}, \code{family}, \code{score_type}, etc.
#'   Supports the \code{.} shorthand in formula updates, e.g.
#'   \code{formula = . ~ . + offset(.OFFSET___)}.
#' @return A \code{flipscores} object.
#' @method update flipscores
#' @export
update.flipscores <- function(object, ...) {

  call <- object$flipscores_call
  if (is.null(call))
    call <- object$call

  extras <- match.call(expand.dots = FALSE)$...

  if (length(extras) > 0) {
    # handle formula update separately using update.formula
    # to correctly resolve '.' shorthand
    if ("formula" %in% names(extras)) {
      call$formula <- stats::update.formula(formula(object),
                                            eval(extras[["formula"]],
                                                 parent.frame()))
      extras[["formula"]] <- NULL
    }
    # override remaining arguments
    for (a in names(extras))
      call[[a]] <- extras[[a]]
  }

  # evaluate in the environment of the original model's terms
  # so that variables referenced in the formula are found correctly
  env <- tryCatch(
    attr(terms(object), ".Environment"),
    error = function(e) NULL
  )
  if (is.null(env)) env <- parent.frame()

  eval(call, envir = env)
}

#' print.flipscores print method for a flipscores object.
#' @param x a flipscores object
#' @method print flipscores
#' @docType methods
#' @rdname flipscores-method
#' @export


print.flipscores <- function(x, ...) {
  cat(get_head_flip_out(x))
  cat("Call: ")
  print(x$flipscores_call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  # print.default(x)
}

#' summary.flipscores summary method for a flipscores object.
#' @rdname flipscores-method
#' @param object a flipscores object
#' @param ... additional arguments to be passed
#' @method  summary flipscores
#' @docType methods
#' @export

summary.flipscores <- function(object, ...) {
  sum_model <- summary.glm(object = object)

  display_call <- object$flipscores_call
  if (!is.null(display_call)) {
    # replace family with a clean deparsed version
    fam <- object$family
    if (!is.null(fam)) {
      display_call$family <- if (is.character(fam)) {
        fam
      } else {
        str2lang(paste0(fam$family, "(link='", fam$link, "')"))
      }
    }
    # only shorten data if it was NOT passed as a named object
    # i.e. if it is not a simple symbol like `df` or `mydata`
    data_arg <- object$flipscores_call$data
    if (!is.null(data_arg) && !is.symbol(data_arg)) {
      display_call$data <- as.symbol(
        paste0("data.frame_", nrow(object$model), "x", ncol(object$model))
      )
    }
  }

  sum_model$coefficients <- sum_model$coefficients[, c(1,1:4,4), drop=FALSE]
  sum_model$coefficients[, -1] <- NA
  sum_model$coefficients[names(object$p.values), -1] <- NA
  sum_model$coefficients[names(object$p.values), 2] <- unlist(colSums(object$scores))
  sum_model$coefficients[names(object$p.values), 3] <- attributes(object$scores)$sd
  sum_model$coefficients[, 4] <- sum_model$coefficients[, 2] / sum_model$coefficients[, 3]
  sum_model$coefficients[names(object$p.values), 5] <-
    (sum_model$coefficients[names(object$p.values), 2] / attributes(object$scores)$nrm)[]
  sum_model$coefficients[names(object$p.values), 6] <- object$p.values
  colnames(sum_model$coefficients)[c(1,2,3,4,5,6)] <- c("estimate","score", "se","z",
                                                        "pcor", "p")
  sum_model$aliased <- rep(FALSE, length(sum_model$aliased))
  sum_model$call <- display_call
  sum_model
}

#' @description \code{p.adjust} method for class "\code{joint_flipscores}" and "\code{flipscores}".
#' Add adjusted p-values into the \code{joint_flipscores}and \code{flipscores} object.
#' @rdname joint_flipscores-methods
#' @param object an object of class \code{joint_flipscores} and \code{flipscores}.
#' @param method any method implemented in \code{flip.adjust} or
#' a custom function. In the last case it must be a function that uses a matrix
#' as input and returns a vector of adjusted p-values equal to the number of columns of the inputed matrix.
#' @param tail argument: expresses the tail direction of the alternative hypothesis.
#' It can be "two.sided" (or 0, the default), "less" (or -1) or "greater" (or +1).
#' @param ... additional arguments to be passed
#' @docType methods
#' @export

p.adjust <- function (object, method = "maxT", tail = 0, ...)
{
  if(is.character(method)){
    if(method=="maxT"){
      #      if("alphas"%in%names(as.list(match.call())))
      p.adj=maxT.light(.set_tail(object$Tspace, tail = tail),...)
    } else
      if(method%in%c("minp","minP","Tippet","Tippett") ) {
        p.adj=maxT.light(-.t2p(object$Tspace, tail = tail),...)
      } else
        p.adj = flip.adjust(.set_tail(object$Tspace, tail = tail),
                            method = method)
  } else if(is.function(method)){
    p.adj = method(.set_tail(object$Tspace, tail = tail))
  }
  object$summary_table$p.adj <- p.adj
  object$p.adjust.method <- method
  object
}

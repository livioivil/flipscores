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
  colnames(sum_model$coefficients)[c(2,4,5,6)] <- c("Score", "z value",
                                                    "Part. Cor", "Pr(>|z|)")
  sum_model$aliased <- rep(FALSE, length(sum_model$aliased))
  sum_model$call <- display_call
  sum_model
}



###########
get_head_flip_out <- function(x){
  if(is.null(dim(x$n_flips)))
    n_flips=x$n_flips else
      n_flips=nrow(x$n_flips)
  if(length(grep("Negative Binomial",x$family$family))==0)
  {paste("Flip Score Test:
         score_type =",x$score_type,
         ", n_flips =",n_flips,"\n")}
  else
    paste("Flip Score Test:
          score_type =",x$score_type,
          ", n_flips =",n_flips,
          ", theta =",round(x$theta,digits=5),"\n")
}

.intercept_in_Z_and_count_family <- function(mod){
  has_intercept=TRUE
  # check if family is count data
  if(mod$family$family%in% c("poisson", "quasipoisson", "Negative Binomial"))
    {
    Z=model.matrix(mod)
    if(ncol(Z)==0) {
      has_intercept=FALSE
      } else {
        P=Z%*%solve(t(Z)%*%Z)%*%t(Z)
        ones=matrix(1,nrow(Z),1)
        has_intercept=!(sum(abs(P%*%ones-ones))>1E-10)
      }
    if(!has_intercept) { #if intercept is the space spanned by Z
    warning("The Normalized Generalized Partial Correlation (Determination) Coefficient for Count families without interncept in the null model has not implemented, yet. NA will be returned.")
    }
  }
  has_intercept
}


.set_tail <- function(Tspace,tail=0){
  if((tail==0)|(tail=="two.sided"))
    Tspace=abs(Tspace) else
      if((tail<0)|(tail=="less"))
        Tspace=-Tspace

      Tspace
}

.t2p_only_first <- function(Tspace,tail=0){
  Tspace=.set_tail(Tspace,tail=tail)
  if(is.vector(Tspace)){
    P = mean(Tspace>=Tspace[1])
  } else if(ncol(Tspace)==1){
    P = mean(Tspace[,]>=Tspace[1,1])
  } else    {
    P = apply(Tspace,2, function(Tsp)mean(Tspace>=Tspace[1]))
  }
  P
}

.find_common_pattern <- function(vettore) {
  if (length(vettore) < 2) {
    return(vettore)
  }

  if(length(grep(":",vettore[1]))>0){
    vettore_splt=strsplit(vettore,":")
    pttrns=sapply(1:length(vettore_splt[[1]]),function(i){
      .find_common_pattern (sapply(vettore_splt,function(x) x[i]) )
    })
    return(paste(pttrns,collapse=":"))
  }

  chars=sapply(vettore,strsplit,"")
  #matrix of all chracters
  charsMat=suppressWarnings(do.call(cbind,chars))
  # ask row-wise if they are all equals
  all_eqs=apply(charsMat,1,function(x)length(unique(x))==1)
  # this is the first different character
  if(all(all_eqs)) {
    common_pattern=vettore[1]
  } else {
    common_pattern=substr(vettore[1],1,max(which.min(all_eqs)-1,1))
  }


  return(common_pattern)
}

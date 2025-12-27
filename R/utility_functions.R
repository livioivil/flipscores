#' @title Methods for flipscores objects
#'
#' @description Methods for \code{flipscores} objects.
#' The following are methods to extract and manipulate relevant information from
#' a \code{flipscores} object.
#'
#' @name flipscores-method
#' @docType methods

NULL



#' print.flipscores print method for a flipscores object.
#' @param x a flipscores object
#' @method print flipscores
#' @docType methods
#' @rdname flipscores-method
#' @export


print.flipscores <- function(x, ...) {
  cat(get_head_flip_out(x))
  cat("Call: ")
  print(x$call)
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

summary.flipscores <- function (object, ...) {
  sum_model=summary.glm(object = object)
  sum_model$coefficients=sum_model$coefficients[,c(1,1:4,4),drop=FALSE]
  sum_model$coefficients[,-1]=NA
  # temp=sum_model$coefficients
  # matrix(NA,length(sum_model$coefficients),6)
  # rownames(temp)=rownames(sum_model$coefficients)
  # colnames(temp)=colnames(sum_model$coefficients)
  # temp[rownames(sum_model$coefficients),]=sum_model$coefficients
  #sum_model$coefficients=temp
  sum_model$coefficients[names(object$p.values),-1]=NA
  sum_model$coefficients[names(object$p.values),2]=unlist(object$Tspace[1,,drop=TRUE])
  sum_model$coefficients[names(object$p.values),3]=attributes(object$scores)$sd#unlist(sapply(object$scores,sd)*sqrt(nrow(object$scores)))
  sum_model$coefficients[,4]=sum_model$coefficients[,2]/sum_model$coefficients[,3]
  sum_model$coefficients[names(object$p.values),5]=(sum_model$coefficients[names(object$p.values),2]/attributes(object$scores)$nrm)[]
  sum_model$coefficients[names(object$p.values),6]=object$p.values
  # sum_model$coefficients=sum_model$coefficients[,c(1,4)]
  colnames(sum_model$coefficients)[c(2,4,5,6)]=c("Score","z value","Part. Cor","Pr(>|z|)")

  sum_model$aliased=rep(FALSE,length(sum_model$aliased))
  structure(sum_model, heading = get_head_flip_out(object), class = c("data.frame"))
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

#' @export
plot.flipscores <- function (x, ...) {
  summary.flipscores(x, ...)
  # object <- x
  # dots <- list(...)
  # a <- object$plotting
  # b <- object$Descriptive
  # fac.names <- a$fac_names
  # exist <- hasArg(factor)
  #
  # if(length(fac.names) != 1){
  #   if(!exist){
  #     print("Please choose the factor you wish to plot (for interaction type something like group1:group2) and confirm by pressing 'Enter'")
  #     Faktor <- scan("", what="character")
  #     while(length(Faktor)==0){
  #       print("Please enter the name of the factor you wish to plot!")
  #       Faktor <- scan("", what="character")
  #     }
  #     } else {
  #      Faktor <- dots$factor
  #     }
  #   x.label <- ""
  # } else {
  #   Faktor <- fac.names
  #   x.label <- fac.names
  # }
  #
  # match.arg(Faktor, fac.names)
  #
  # # default values
  # args <- list(plot.object = a, descr.object = b, factor = Faktor,
  #              lwd =2, ylab = "Means", xlab = x.label, col = 1:length(fac.names), pch = 1:18, legendpos = "topright")
  #
  # args[names(dots)] <- dots
  #
  # do.call(plotting, args = args)
}

get_head_flip_out <- function(x){
  if(stringr::word(x$family$family, 1)!="Negative")
  {paste("Flip Score Test: 
         score_type =",x$score_type,
         ", n_flips =",x$n_flips,"\n")}
  else 
    paste("Flip Score Test: 
          score_type =",x$score_type,
          ", n_flips =",x$n_flips,
          ", theta =",round(x$theta,digits=5),"\n")
}

#' @export
print.flipscores <- function(x, ...) {
  cat(get_head_flip_out(x))
  cat("Call: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  # print.default(x)
}

#' @export
summary.flipscores <- function (object, ...) {
  sum_model=summary.glm(object = object)
  sum_model$coefficients=sum_model$coefficients[,c(1,1:4)]
  sum_model$coefficients[,5]=object$p.values
  sum_model$coefficients[,2]=object$Tspace[1,]
  sum_model$coefficients[,3]=apply(object$scores,2,sd)/sqrt(nrow(object$scores))
  sum_model$coefficients[,4]=sum_model$coefficients[,2]/sum_model$coefficients[,3]
  # sum_model$coefficients=sum_model$coefficients[,c(1,4)]
  colnames(sum_model$coefficients)[c(2,4)]=c("Score","z value")
  cat(get_head_flip_out(object))
  
  sum_model
}

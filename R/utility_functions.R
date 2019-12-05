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


#' print.flipscores

print.flipscores <- function(x, ...) {
  cat(get_head_flip_out(x))
  cat("Call: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  # print.default(x)
}

#' summary.flipscores
#' @export summary.flipscores
summary.flipscores <- function (object, ...) {
  sum_model=summary.glm(object = object)
  sum_model$coefficients=sum_model$coefficients[,c(1,1:4)]
  sum_model$coefficients[,5]=object$p.values
  sum_model$coefficients[,2]=object$Tspace[1,]
  sum_model$coefficients[,3]=apply(object$scores,2,sd)*sqrt(nrow(object$scores))
  sum_model$coefficients[,4]=sum_model$coefficients[,2]/sum_model$coefficients[,3]
  # sum_model$coefficients=sum_model$coefficients[,c(1,4)]
  colnames(sum_model$coefficients)[c(2,4)]=c("Score","z value")
  
  structure(sum_model, heading = get_head_flip_out(object), class = c("data.frame"))
  sum_model
}

#' anova.flipscores
#' @description This is the \code{anova} method for \code{flipscores} object. Remark: it performs type III deviance decomposition as in \code{car::Anova}.
#' @param object (the object) \code{glm} (or \code{flipscores}) object with the model under the null hypothesis (i.e. the covariates, the nuisance parameters).
#' @param model1 a \code{glm} (or \code{flipscores}) or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The variables in \code{model1} are the same variables in \code{object} plus one or more variables to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested variables column-wise.
#' @param score_type The type of score that is computed. It can be "orthogonalized", "effective" or "basic".
#' Default is "orthogonalized". "effective" and "orthogonalized" take into account the nuisance estimation. The default is \code{NULL}, in this case the value is taken from \code{object}.
#' @param n_flips The number of random flips of the score contributions.
#' When \code{n_flips} is equal or larger than the maximum number of possible flips (i.e. n^2), all possible flips are performed. 
#' Default is 5000.
#' @param id a \code{vector} identifying the clustered observations. If \code{NULL} (default) observations are assumed to be independent. NOTE: if \code{object} is a \code{flipscores} and \code{model$flip_param_call$id} is not \code{NULL}, this is considered in the inference.
#' @param ... other parameters allowed in \code{stats::anova}.
#' @examples
#' set.seed(1)
#' dt=data.frame(X=scale(rnorm(50)),
#'    Z=factor(rep(LETTERS[1:3],length.out=50)))
#' dt$Y=rpois(n=nrow(dt),lambda=exp(dt$X*(dt$Z=="C")))
#' mod0=flipscores(Y~Z+X,data=dt,family="poisson",score_type = "standardized")
#' summary(mod0)
#' anova(mod0)
#'
#' mod1=flipscores(Y~Z*X,data=dt,family="poisson",score_type = "standardized")
#' summary(mod1)
#' anova(mod0,model1 = mod1)
#'
#' @export


anova.flipscores <- function(object, model1=NULL,
                           score_type=NULL,
                           n_flips=5000, id=NULL,
                           ...){
  
  if(is.null(object$x)||(length(object$x)==0)) object=update(object,x=TRUE)
  
  if(is.null(score_type)) score_type = object$score_type
  
  anova_temp=get("anova.glm", envir = asNamespace("stats"),
                 inherits = FALSE)
  
  ## comparison of 2 nested models
  if(!is.null(model1)){ 
    scores=compute_scores(model0 = object,
                          model1 = model1,
                          score_type = score_type)
    mf <- match.call(expand.dots = TRUE)
    if(!is.null(mf$flip_param_call$id))
      scores=rowsum(scores,group = id)
    
    Tspace=sapply(1:ncol(scores), function(id_col){
      score1=scores[,id_col,drop=FALSE]
      attributes(score1)$scale_objects=attributes(scores)$scale_objects[[id_col]]
      as.matrix(.flip_test(score1,perms=n_flips,
                                        score_type = score_type)$Tspace)
      })
    dst=mahalanobis_npc(Tspace)
    
    out_param=anova_temp(object,model1,test="Rao")
    heading2=attributes(out_param)$heading[2]
    out_param = out_param[-1,]
    out_param = out_param[,-c(1:2,4)]
    names(out_param)[2]="Score"
    out_param[[2]]=dst[1]
    names(out_param)[3]="Pr(>Score)"
    out_param[[3]]=.t2p(dst)
    rownames(out_param)[1]="Model 2 vs Model 1"
    
    } else   { ## one anova for all variables
    
    varlist <- attr(object$terms, "variables")
    if (!is.matrix(object[["x"]])) 
      object$x = model.matrix(object)
    varseq <- attr(object$x, "assign")
    
    ############################
    subsets_npc=lapply(unique(varseq[varseq!=0]),function(i)which(varseq==i))

    res=mahalanobis_npc_multi(ids_list = subsets_npc,permT = as.matrix(object$Tspace))
    # flip::npc(ps@permT,comb.funct = "mahalanobist",subsets = subsets_npc)
    # res@res[, 3]=res@res[, 3]*nrow(ps@permT)
    
    out_param = anova_temp(object,test="Rao")
    heading2=paste0("Model: ",deparse1(formula(object)))
    out_param = out_param[-1,]
    out_param = out_param[,-(2:4)]
    names(out_param)[2]="Score"
    out_param[[2]]=res[1,]
    names(out_param)[3]="Pr(>Score)"
    out_param[[3]]=apply(res,2,.t2p)
  } # closes if
  #make up
  title <- paste0("Analysis of Deviance Table (Type III test)", 
                  "\n\nModel: ", 
                  object$family$family, ", link: ", object$family$link,"\n") 
  attr(out_param,"heading")=c(title,heading2)
  
  
  attr(out_param,"heading")[[1]]=  paste(attr(out_param,"heading")[[1]],sep="",
          "\nInference is provided by FlipScores approach (",object$n_flips," sign flips).\n")
  return(out_param)
}
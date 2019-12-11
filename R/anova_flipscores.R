#' anova.flipscores
#' @param model0 a \code{glm} (or \code{flipscores}) object with the model under the null hypothesis (i.e. the covariates, the nuisances parameters).
#' @param model1 a \code{glm} (or \code{flipscores}) or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The variables in \code{model1} are the same variables in \code{model0} plus one or more variables to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested variables column-wise.
#' @param score_type The type of score that is computed, either "orthogonalized", "effective" or "basic". 
#' By default is "orthogonalized". "effective" and "orthogonalized" takes into account nuisance estimation.
#' @param n_flips The number of random flips of the scores contributions
#' When \code{n_flips} is equal or larger than the maximum number of possible flips (i.e. n^2), all possible flips are performed. 
#' Default is 5000.
#' @param id a \code{vector} identifying the clustered observations. If \code{NULL} (default) observations are assumed to be independent. NOTE: if \code{model0} is a \code{flipscores} and \code{model$id} not \code{NULL}, this is considered in the inference.
#' @param type type of test, "I", "III", 1, or 3. Roman numerals are equivalent to the corresponding Arabic numerals.
#' @examples
#'
#' set.seed(1)
#' Z=rnorm(20)
#' X=factor(rep(LETTERS[1:3],length.out=20))
#' Y=rpois(n=20,lambda=exp(Z*(X=="C")*.5))
#' mod0=flipscores(Y~Z+X,family="poisson")
#' summary(mod0)
#' anova(model0 = mod0)
#'
#' mod1=flipscores(Y~Z*X,family="poisson")
#' anova(model0 = mod0,mod1)
#'
#' @export

anova.flipscores <- function(model0, model1=NULL,
                           score_type="orthogonalized",
                           n_flips=5000, type=3,id=NULL,
                           ...){
  
  if(type=="I") type=1
  if(type=="III") type=3
  if(!(type %in%c(1,3))) {
    warning("type=",type," is not implemented. Forcing to type=3")
    type=3
  }
  
  if(is.null(model0$x)||(length(model0$x)==0)) model0=update(model0,x=TRUE)
  
  ## comparison of 2 nested models
  if(!is.null(model1)){ 
    scores=compute_scores(model0=model0,
                          model1=model1,
                          score_type=score_type)
    mf <- match.call(expand.dots = TRUE)
    if(!is.null(mf$id))
      scores=rowsum(scores,group = id)
    ps=flip::flip(scores,perms=n_flips)
    res=flip::npc(ps@permT,comb.funct = "mahalanobist")
    out_param=stats:::anova.glm(model0,model1,test="Rao")
    out_param$Rao[2]=res@res[1,3]
    out_param$`Pr(>Chi)`[2]=res@res[1,4]
    
    } else   { ## type I or III deviance decomposition
    
    varlist <- attr(model0$terms, "variables")
    if (!is.matrix(model0$x)) 
      model0$x = model.matrix(model0)
    varseq <- attr(model0$x, "assign")
    nvars <- max(0, varseq)
    resdev <- resdf <- NULL
  
    mf <- match.call(expand.dots = TRUE)
    if(!is.null(mf$id)) 
      model0$id=mf$id
    
    scores=lapply(1:nvars, function(var_i) {
      tested=which(varseq==var_i)
      # print(tested)
      excluded <- if(type==2) which(varseq>var_i) else c()
      # print(excluded)
      flipscores:::socket_compute_scores(tested,model0,exclude=excluded,score_type=score_type)
      })
    temp=sapply(scores,ncol)
    subsets_npc=list(1:temp[1])
    if(length(temp)>1)
      for(i in 2:length(temp))
        subsets_npc=c(subsets_npc,list(sum(temp[1:(i-1)])+(1:temp[i])))
        
    scores = as.data.frame(scores)

  
    ps=flip::flip(as.matrix(scores),perms=n_flips,tail=1)
    res=flip::npc(ps@permT,comb.funct = "mahalanobist",subsets = subsets_npc)
    # res@res[, 3]=res@res[, 3]*nrow(ps@permT)
    out_param=car::Anova(model0, test = "LR",type=type) 
    #stats:::anova.glm(model0,test="Rao")
    # out_param$Rao[-1]=res@res[,3]
    # out_param$`Pr(>Chi)`[-1]=res@res[,4]
    names(out_param)[1]="Score"
    out_param[[1]]=res@res[,3]
    names(out_param)[3]="Pr(>Score)"
    out_param[[3]]=res@res[,4]
    if(type==1){
      type_test=": Type I test (i.e. terms added sequentially, first to last)"
    } else if (type==3) {
      type_test=": Type III test"
    }
    title <- paste0("Analysis of Deviance Table", type_test,
                    "\nModel: ", 
                    model0$family$family, ", link: ", model0$family$link, 
                    "\nResponse: ", as.character(varlist[-1L])[1L])
    title=gsub("Terms added sequentially .first to last.\\n\\n",
                                        "",title)
    attr(out_param,"heading")[[1]]=title
    
  } # closes if
  #make up
  
  attr(out_param,"heading")[[1]]=  paste(attr(out_param,"heading")[[1]],sep="",
          "Inference is provided by FlipScores approach (",model0$n_flips," sign flips).\n")
  return(out_param)
}
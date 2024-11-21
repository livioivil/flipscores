#' compute_scores
#' @param model0 a \code{glm} object with the model under the null hypothesis (i.e. the covariates, the nuisance parameters).
#' @param model1 a \code{glm} or a \code{matrix} (or \code{vector}). If it is a \code{glm} object, it has the model under the alternative hypothesis. The variables in \code{model1} are the same variables in \code{model0} plus one or more variables to be tested.  Alternatively, if
#' \code{model1} is a \code{matrix}, it contains the tested variables column-wise.
#' @param score_type The type of score that is computed. It is "orthogonalized", "effective" or "basic".
#' "effective" and "orthogonalized" take into account the nuisance estimation.
#' @param ... other arguments.
#' @description 
#' Same usage as \code{anova.glm}. 
#' The parameter \code{id}  is used too, 
#' if present in \code{model0} (with priority) or in \code{model1}.
#'
#' @author Jesse Hemerik, Riccardo De Santis, Vittorio Giatti, Jelle Goeman and Livio Finos
#' @examples
#' set.seed(1)
#' Z=rnorm(20)
#' X=Z+rnorm(20)
#' Y=rpois(n=20,lambda=exp(Z+X))
#' mod0=glm(Y~Z,family="poisson")
#' X=data.frame(X=X)
#' scr0=compute_scores(model0 = mod0, model1 = X)
#' head(scr0)
#' @export

compute_scores <- function(model0, model1,score_type = "standardized",...){
  score_type=match.arg(score_type,c("orthogonalized","standardized","effective","basic","my_lab"))
  if(missing(score_type))
    stop("test type is not specified or recognized")
  X=get_X(model0,model1)
  if(!is.null(model0))  
  {id_model1=if(is.list(model1)) eval(model1$flip_param_call$id) 
  else 
    NULL
  } else id_model1=NULL
  
  rm(model1)
  
  ###############
  if(is.null(model0$x)||is.null(ncol(model0$x))) {
    call <- getCall(model0)
    call$x=TRUE
    term <- terms(model0)
    env <- attr(term, ".Environment")
    model0=eval(call, env, parent.frame())
  }
  
  
  # no variables in the null model
  if(ncol(model0$x)==0){
    scores=matrix(model0$y)
    rownames(scores)=names(model0$fitted.values)
    scale_objects=list(list(nrm = 1))
    score_type="basic"
    sqrtinvV_vect_times_residuals = rep(1,nrow(model0$x))
  } else { # at least one covariate
    if(is.null(model0$y)) model0$y=model0$model[,1]
    Z=model0$x
    residuals=(model0$y-model0$fitted.values)
    if(is.null(model0$weights))  sqrtW=rep(1,length(residuals)) else sqrtW=(as.numeric(model0$weights)**0.5)
    pars_score<-get_par_expo_fam(model0)
    D_vect<-pars_score$D
    V_vect<-pars_score$V
    sqrtinvV_vect<-V_vect**(-0.5)
    sqrtinvV_vect_times_residuals=sqrtinvV_vect*residuals
    
    #BASIC SCORE
    if (score_type == "basic") {
      .get_1score_basic <- function(X){
        B=t(t(X * D_vect) %*% (diag(sqrtinvV_vect**2, nrow = length(sqrtW))))
        scores = B * (sqrtinvV_vect_times_residuals) #* (1/sum(!is.na(model0$y))^0.5)
        scale_objects=list(nrm=sqrt(sum(B^2)*sum((sqrtinvV_vect_times_residuals)^2)))
        list(scores=scores, scale_objects=scale_objects)
      }
  
      temp=apply(X,2,.get_1score_basic)
      scores=sapply(temp,function(obj) obj$scores)
      # print(names(scores))
      scale_objects=lapply(temp,function(obj) obj$scale_objects)
      # print(names(scale_objects))
    } else
      
      ###############################
      ######################## not a documented score type, just a place for experiments
      if (score_type == "my_lab") {
      .get_1score_orthogonalized <- function(X){
          B=(t(X*sqrtW)%*%OneMinusH*(sqrtinvV_vect))
          dinv=deco$d^-1
          Xt=(B%*%deco$u%*%diag(deco$d))
          scores=(diag(dinv)%*%t(deco$u)%*%(residuals))[,]#*(1/sum(!is.na(model0$y))**0.5)
          nrm=sqrt(sum(B^2)*sum(residuals^2))
          scale_objects=list(U=deco$u,B=B,nrm=nrm,Xt=Xt)
         
          list(scores=scores, scale_objects=scale_objects)
        }
        
        OneMinusH = diag(nrow(Z)) - ((sqrtW)* Z) %*% solve(t(Z) %*% ((sqrtW**2) * Z)) %*% t(Z * (sqrtW))
        deco=svd((V_vect^0.5)*OneMinusH,nv = 0)
        deco$u=deco$u[,deco$d>1E-12]
        deco$d=deco$d[deco$d>1E-12]
        temp=apply(X,2,.get_1score_orthogonalized)
        scores=sapply(temp,function(obj) obj$scores)
        # print(names(scores))
        scale_objects=lapply(temp,function(obj) obj$scale_objects)
        # print(names(scale_objects))
  
        
        attr(scores,"scale_objects")=scale_objects
        rownames(scores)=1:length(scores)
        return(scores)
        
        } else
      ##  EFFECTIVE SCORE 
      if(score_type=="effective"){
        .get_1score_effective <- function(X){
          B<-X*(sqrtW)-t(crossprod(crossprod(A,X*(sqrtW)),solve(crossprod(A),t(A))))
          scores=B*sqrtinvV_vect_times_residuals#/(sum(!is.na(model0$y))**0.5)
          scale_objects=list(nrm=sqrt(sum(B^2)*sum((sqrtinvV_vect_times_residuals)^2)))
          list(scores=scores, scale_objects=scale_objects)
        }
        
        A<-(sqrtW)*Z
        temp=apply(X,2,.get_1score_effective)
        scores=sapply(temp,function(obj) obj$scores)
        # print(names(scores))
        scale_objects=lapply(temp,function(obj) obj$scale_objects)
        # print(names(scale_objects))
        
      } else
        ##  STANDARDIZED SCORE
        if(score_type=="standardized"){
          .get_1score_standardized <- function(X,U){
            b=crossprod(diag(nrow(Z))-tcrossprod(U),X*sqrtW)
            m = sum(b^2)
            # we divide it by sqrt(m) which is the sd scaling factor of the observed test stat (i.e. effective and standardized have the same observed test stat)
            A=b[,]*U/sqrt(m)
            scores=b*sqrtinvV_vect_times_residuals#/(sum(!is.na(model0$y))**0.5)
            nrm=sqrt(sum(b^2)*sum((sqrtinvV_vect_times_residuals)^2))
            scale_objects=list(A=A,nrm=nrm)
            list(scores=scores, scale_objects=scale_objects)
          }
          
          U=svd((sqrtW*Z),nv=0)$u
          temp=apply(X,2,.get_1score_standardized,U)
          scores=sapply(temp,function(obj) obj$scores)
          # print(names(scores))
          scale_objects=lapply(temp,function(obj) obj$scale_objects)
          # print(names(scale_objects))
          
        } else  
          #ORTHO EFFECTIVE SCORE
          if(score_type=="orthogonalized"){
            .get_1score_orthogonalized <- function(X){
              B=(t(X*sqrtW)%*%OneMinusH*(sqrtinvV_vect))
              scores=t(B%*%deco$u)*(t(deco$u)%*%(residuals))[,]#*(1/sum(!is.na(model0$y))**0.5)
              nrm=sqrt(sum(B^2)*sum(residuals^2))
              scale_objects=list(U=deco$u,B=B,nrm=nrm)
              list(scores=scores, scale_objects=scale_objects)
            }
            
            OneMinusH = diag(nrow(Z)) - ((sqrtW)* Z) %*% solve(t(Z) %*% ((sqrtW**2) * Z)) %*% t(Z * (sqrtW))
            deco=svd((V_vect^0.5)*OneMinusH,nv = 0)
            deco$d[deco$d<1E-12]=0
            temp=apply(X,2,.get_1score_orthogonalized)
            scores=sapply(temp,function(obj) obj$scores)
            # print(names(scores))
            scale_objects=lapply(temp,function(obj) obj$scale_objects)
            # print(names(scale_objects))
            
   
          }
  }
  
  rownames(scores)=names(sqrtinvV_vect_times_residuals)
  
  nobservations=list(...)$nobservations
  if(!is.null(nobservations))
    if(nrow(scores)<nobservations){
      if(score_type=="standardized"){
        for(i in 1:length(scale_objects)){
          temp=matrix(0,nobservations,ncol(scale_objects[[i]]$A))
          temp[as.numeric(rownames(scores)),]=scale_objects[[i]]$A
          rownames(temp)=1:nobservations
          scale_objects[[i]]$A=temp
        }
      }
      
      temp=matrix(0,nobservations,ncol(scores))
      temp[as.numeric(rownames(scores)),]=scores
      scores=temp
      rownames(scores)=1:nobservations
      
      temp=rep(0,nobservations) 
      names(temp)=1:nobservations
      temp[as.numeric(names(sqrtinvV_vect_times_residuals))]=sqrtinvV_vect_times_residuals
      sqrtinvV_vect_times_residuals=temp
  
    }
  
  std_dev=get_std_dev_score(model0,X)
#  scale_objects$df.residual <- df.residual(model0)
  attr(scores,"sd")=std_dev
  attr(scores,"scale_objects")=scale_objects
  attr(scores,"score_type")=score_type
  attr(scores,"resid_std")=sqrtinvV_vect_times_residuals
  
  return(scores)
}

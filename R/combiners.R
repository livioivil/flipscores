#' @title Nonparametric combination of \code{joint_flipscores} objects
#' @description Methods for combining \code{joint_flipscores} objects.
#' @docType methods
#' @name combine-methods
#' @param mods a \code{joint_flipscores} object.
#' @param comb_funct  combining function to be used.
#' Several functions are implemented: "mean", "median", "Fisher", "Liptak", (equal to) "Stoufer", "Tippet", (equal to) "minp", "maxT", "Mahalanobis".
#' Alternatively it can be a custom function that has a Tspace matrix as input.
#' For \code{combine_tests} the default is \code{comb_funct="maxT"},
#' while for \code{combine_contrasts} the default is \code{comb_funct="Mahalanobis"}.
#' @returns The function returns a \code{joint_flipscores}-object.
#' @param by if \code{NULL} (default), it combines all test statistics.
#' If a characters, it refers to the column's name of \code{summary_table} (and printed by something like \code{summary(mods)}).
#' The elements with the same value will be combined. If \code{by} is a vector, the values are defined by row-wise concatenation of the values of the columns in \code{by}.
#' The argument is inactive if \code{by_list} is not \code{NULL}.
#' @param by_list NULL (default) or a list of vectors. For each vector of the list it combines test statistics with position given by the element of the vector.
#' If the vectors in the list are characters, these refer to names(mods$Tspace).
#' @param tail direction of the alternative hypothesis. It can be "two.sided" (or 0, the default), "less" (or -1) or "greater" (or +1).



#' @usage combine_tests(mods, comb_funct = "maxT", by = NULL, by_list=NULL, tail = 0)
#' @description \code{combine_tests} combines the tests derived from multiverse models.
#' @docType methods
#' @rdname combine_tests
#' @export
#' @examples
#' #First example
#' library(flipscores)
#' set.seed(123)
#'
#' #Simulate data
#' n=20
#' D=data.frame(X=rnorm(n),Z1=rnorm(n),Z2=rnorm(n))
#' D$Y=D$Z1+D$X+rnorm(n)
#'
#' # Run four glms abd combine it in a list
#' mod1=glm(Y~X+Z1+Z2,data=D)
#' mod2=glm(Y~X+poly(Z1,2)+Z2,data=D)
#' mod3=glm(Y~X+poly(Z1,2)+poly(Z2,2),data=D)
#' mod4=glm(Y~X+Z1+poly(Z2,2),data=D)
#' mods=list(mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
#'
#' # Let us analyze the tests related to coefficient "X" and combine them
#' res=flipscores(mods,n_flips = 5000, seed = 1, tested_coeffs = "X")
#' summary(combine_tests(res))


combine_tests <- function (mods, comb_funct = "maxT", by = NULL, by_list=NULL, tail = 0)
{

  # for(i in 1:length(mods))
  #     mods[[i]]$call$data=mods[[i]]$data

  if(!is.null(by_list)){
    combined=by_list
  } else
  { # by_list is not active, let's use by

    # names(mods) = .set_mods_names(mods)
    if(is.null(by)){ # combine overall
      combined = list(overall = 1:ncol(mods$Tspace))

    } else { # more complex combinations
      smr=apply(mods$summary_table[,by,drop=FALSE],1,paste,collapse=".")
      uniq_nm = unique(smr)
      combined = lapply(uniq_nm, function(nm) which(smr == nm))
      names(combined) = uniq_nm
    }
  }
  res = lapply(1:length(combined), .npc2jointest,
               mods = mods, combined = combined, tail = tail, comb_funct = comb_funct)
  names(res) = names(combined)
  res=list(Tspace=.get_all_Tspace(res),summary_table=.get_all_summary_table(res))
  class(res) <- unique(c("jcombined", "joint_flipscores", class(res)))
  res
}


#' @usage combine_contrasts(mods, comb_funct = "Mahalanobis", tail = 0)
#' @description \code{combine_contrasts} combines the tests derived from the contrasts of a factor variable to get a
#' global test for the factor (i.e. categorical predictor).
#' It has strong analogies with ANOVA test.
#' @docType methods
#' @rdname combine_tests
#' @export
#' @examples
#' # Second (continued) example
#' # flipscores jointly on all models and all coefficients
#' mod1=glm(Y~X+Z1+Z2,data=D)
#' mod2=glm(Y~X+poly(Z1,2)+Z2,data=D)
#' mod3=glm(Y~X+poly(Z1,2)+poly(Z2,2),data=D)
#' mod4=glm(Y~X+Z1+poly(Z2,2),data=D)
#' mods=list(mod1=mod1,mod2=mod2,mod3=mod3,mod4=mod4)
#' res=flipscores(mods,n_flips = 2000)
#' summary(combine_tests(res))
#' summary(combine_tests(res, by="model"))
#' summary(combine_tests(res, by="coefficient"))
#' res2=combine_contrasts(res)
#' summary(res2)
#' #custom combinations:
#' coeffs=c("(Intercept)","X","Z1","Z2")
#' coeffs_ids=lapply(coeffs,grep,res2$summary_table$coeff)
#' names(coeffs_ids)=coeffs
#' summary(combine_tests(res2,by_list =   coeffs_ids))


combine_contrasts <- function (mods, comb_funct = "Mahalanobis", tail = 0)
{
  names(mods) = .set_mods_names(mods)
   if(is.null(mods$summary_table$model)){
    suffix <-as.integer(factor(sub(".*_model\\.*", "", colnames(mods$Tspace))))
    mods$summary_table$model <- paste0("model", suffix)
  }
  smr=apply(mods$summary_table[,c("model",".assign"),drop=FALSE],1,paste,collapse=".")
  new_names=sapply(unique(smr),function(x).find_common_pattern(mods$summary_table$coefficient[smr==x]))
  res=combine_tests(mods,by=c("model",".assign"),comb_funct = comb_funct, tail = tail)
  res$summary_table$coefficient=new_names
  assigns=paste0(".",unique(mods$summary_table$.assign),"$")
  as_ids=lapply(assigns, function(as)  grep(as,res$summary_table$Model))
  for(i in 1:length(assigns)){
    res$summary_table$Model[as_ids[[i]]]=gsub(assigns[i],"",res$summary_table$Model[as_ids[[i]]])
  }
  res
}
#######################
.npc2jointest <- function (id, mods, combined, tail, comb_funct)
{
  comb_name = names(combined)[id]
  if (is.null(comb_name))
    comb_name = "combined"
  Tspace = npc(mods$Tspace[, combined[[id]],drop=FALSE], comb_funct = comb_funct,
               tail = tail)
  colnames(Tspace) = comb_name
  Coeff=mods$summary_table[combined[[id]],"coefficient"]
  Coeff=unique(Coeff)
  if(length(Coeff) > 1) Coeff="many"
  summary_table = data.frame(coefficient = Coeff, stat = comb_funct,
                             ntests = max(1, length(combined[[id]])), S = Tspace[1],
                             p = .t2p_only_first(Tspace, tail = 1))
  list(Tspace = Tspace, summary_table = summary_table)
}


#####################
#' @importFrom stats qnorm
#' @importFrom stats median
#' @importFrom stats as.formula
#' @importFrom stats coefficients
#' @importFrom stats formula
#' @importFrom stats glm
#' @importFrom stats model.matrix
#' @importFrom stats quantile
#' @importFrom stats reformulate
#' @importFrom stats terms
#' @importFrom stats update

npc <- function (Tspace, comb_funct = "Fisher", tail = 0)
{
  if(is.character(comb_funct)){
    implemented_comb_funcs=c("mean", "median","Fisher",
                             "Liptak", "Stoufer",
                             "Tippet", "minp","maxT","Mahalanobis")
    comb_funct=match.arg(comb_funct,implemented_comb_funcs)
    if (comb_funct %in% c("mean", "median","maxT"))
      Tspace = .set_tail(Tspace, tail = tail)

    if (comb_funct %in% c("Fisher",
                          "Liptak", "Stoufer",
                          "Tippet", "minp"))
      Tspace = .t2p(Tspace, tail = tail)

    if (comb_funct == "minp")
      Tspace = .rowMin(Tspace)
    else if (comb_funct %in% c("maxT")) {
      Tspace = .rowMax(Tspace)
    } else if (comb_funct == "mean") {
      Tspace = rowMeans(Tspace)
      Tspace = .set_tail(Tspace, tail = tail)
    }
    else if (comb_funct == "median") {
      Tspace = .rowMedians(Tspace)
      Tspace = .set_tail(Tspace, tail = tail)
    }
    else if (comb_funct == "Fisher") {
      Tspace = rowSums(.comb_funct_fisher(Tspace))
    }
    else if (comb_funct %in% c("Liptak", "Stoufer")) {
      Tspace = rowSums(.comb_funct_liptak(Tspace))
    } else if (comb_funct %in% c("Mahalanobis")) {
      Tspace = mahalanobis_npc(Tspace)
    }
  } else if (is.function(comb_funct)) {
    Tspace = comb_funct(Tspace)
    comb_funct = "custom"
  }
  Tspace = matrix(Tspace)
  Tspace
}


.rowMedians <- function(X,...)
  apply(X,1,median,...)
.rowMax <- function(X,...)
  apply(X,1,max,...)
.rowMin <- function(X,...)
  apply(X,1,min,...)

.comb_funct_fisher <- function(p)-log(p)
.comb_funct_liptak <- function(p)-qnorm(p)
#add here other functions..
################################

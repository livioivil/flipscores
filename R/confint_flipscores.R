#' Confidence intervals for flipscores
#' @description This function allows to apply the \code{confint} method to \code{flipscores} objects.
#' @param object a \code{flipscores}) model object.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, they are recovered from the input model (\code{to_be_tested}).
#' @param level the confidence level required.
#' @param type the type of confidence interval to build, either "equitailed" or "symmetric". Default is "equitailed".
#' @param score_type The type of score that is computed. It can be "standardized" or "effective".
#' Default is inherited by the input model.
#' @param alternative The direction of the test to invert, and corresponding interval with respect to the estimate. It can be "two.sided" (default), "greater", or "less".
#' @param flips the fixed flips to use in the interval search. If \code{NULL} (default), it tries to inherit it from the input model.
#' @param n_flips If \code{flips} are not given or retrieved by the model, the number of random flips of the score contributions. Default inherits them from the input model.
#' @param ... additional arguments for methods.
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%), or "Est." and level or (1-level) if a one-sided alternative is selected.
#'
#' @details All families allowed in \code{glm}, as well as Negative Binomial, are supported, and inherited here by the input model \code{object}. See the main function \code{\link{flipscores}} for further information.
#'
#' \code{flipscores} is based on randomly sign-flipped score contributions. Performing the test (and generating the confidence interval) different times might give slightly different results. Providing a flip matrix generated with \code{\link{make_flips}} is suggested to uniform the results.
#'
#' \code{to_be_tested} works differently for factors and character variables: the corresponding dummy variables are to be tested separately. So, to test a numeric variable \code{X}, one could write \code{to_be_tested = "X"}, but to test a factor \code{Z} with levels \code{A}, \code{B} and \code{C} one should write \code{to_be_tested = c("ZB", "ZC")} as they appear in the model summary.
#'
#'#' @seealso \code{\link{flipscores}}, \code{\link{make_flips}}
#'
#' @examples
#' set.seed(1)
#' dt <- data.frame(X=scale(rnorm(50)),
#'    Z=factor(rep(LETTERS[1:3],length.out=50)))
#' dt$Y=rpois(n=nrow(dt),lambda=exp(dt$X*(dt$Z=="C")))
#' mod0 <- flipscores(Y~Z+X,data=dt,family="poisson")
#' summary(mod0)
#' confint(mod0)
#'
#' xx <- rnorm(20)
#' zz <- rnorm(20, 0.2 * xx)
#' yy <- rnorm(20, 1 + 2 * xx - zz)
#' flps <- make_flips(20, 5000)
#' mod1 <- flipscores(yy ~ xx + zz)
#' summary(mod1)
#' confint(mod1, flips = flps)
#' @export

confint.flipscores <- function(object,
                               parm = NULL,
                               level = 0.95,
                               type = c("equitailed", "symmetric"),
                               score_type = c("standardized", "effective"),
                               alternative = c("two.sided", "greater", "less"),
                               # -1, 0, 1
                               flips = NULL,
                               n_flips = NULL,
                               ...) {
    {
        # robaccia per i default
        if (is.null(parm))
        {
            parm <- names(object$p.values)
            if (is.null(parm))
            {
                parm <- 1:length(coef(object))
            }

        }

        if (all(type == c("equitailed", "symmetric")) ||
            (substr("equitailed", 0, nchar(type)) == type))
        {
            type <- "equitailed"
        } else
        {
            if ((substr("symmetric", 0, nchar(type)) == type))
            {
                type <- "symmetric"
            } else
            {
                stop("Invalid type.")
            }
        }

        if (all(score_type == c("standardized", "effective")))
        {
            score_type <- object$score_type
        } else
        {
            if (substr("standardized", 0, nchar(score_type)) == score_type)
            {
                score_type <- "standardized"
            } else
            {
                if ((substr("effective", 0, nchar(score_type)) == score_type))
                {
                    score_type <- "effective"
                } else
                {
                    stop("Invalid score type.")
                }
            }
        }
        # matlab convention
        if(length(alternative) == 1)
        {
            if(alternative == -1)
                alternative <- "less"
            if(alternative == 0)
                alternative <- "two.sided"
            if(alternative == 1)
                alternative <- "greater"
        }
        if (all(alternative == c("two.sided", "greater", "less")))
        {
            alternative <- "two.sided" # default se non c'è modo di prendersi l'alternative
        } else
        {
            if (substr("two.sided", 0, nchar(alternative)) == alternative)
            {
                alternative <- "two.sided"
            } else
            {
                if ((substr("greater", 0, nchar(alternative)) == alternative))
                {
                    alternative <- "greater"
                } else
                    if ((substr("less", 0, nchar(alternative)) == alternative))
                    {
                        alternative <- "less"
                    } else
                    {
                        stop("Invalid alternative.")
                    }
            }
        }

        if (type == "symmetric" & alternative != "two.sided")
        {
            warning(
                "Symmetric confidence bounds make no sense with a one-sided alternative -- two.sided used instead"
            )
            alternative <- "two.sided"
        }

        if (is.null(flips))
            # if not given in input
        {
            flips <- tryCatch({
                # use the same flips from the model, if in the environment
                eval(object$call$flips, parent.frame()) # attenzione
            }, error = function(e)
            {
                NULL
            })
            if (is.null(flips))
                # if not given to the model, create new
            {
                if (is.null(n_flips))
                {
                    n_flips <- eval(object$flip_param_call$n_flips,
                                    parent.frame()) # try to get it from the model
                }
                if (is.environment(eval(object$data)) | is.null(object$data)) # formula with elements from the environment (i.e. "no prebuilt dataframe") or for negative binomial modeling
                {
                    .form___ <- all.vars(object$terms)[1]
                    n_obs <- length(eval(parse(text = .form___)))
                } else
                {
                    n_obs <- NROW(eval(object$data))
                }
                flips <- flipscores::make_flips(n_obs = n_obs, n_flips = n_flips)
            }
        }
    }

    alpha <- 1 - level

    ci <- lapply(
        parm,
        eval(parse(text = paste0(
            "conf_bound_", type
        ))),
        mod = object,
        flips = flips,
        alpha = alpha,
        score_type = score_type,
        alternative = alternative
    )
    ci <- t(do.call(cbind, ci))
    # 4. Format output as a matrix with appropriate column names
    # ci <- cbind(lb, ub)
    if (alternative == "two.sided")
    {
        colnames(ci) <- paste(format(
            100 * c(alpha / 2, 1 - alpha / 2),
            trim = TRUE,
            scientific = FALSE,
            digits = 3
        ),
        "%")
    }
    if (alternative == "greater")
    {
        colnames(ci) <- c("Est.", paste(
            format(
                100 * c(1 - alpha),
                trim = TRUE,
                scientific = FALSE,
                digits = 3
            ),
            "%"
        ))
    }
    if (alternative == "less")
    {
        colnames(ci) <- c(paste(
            format(
                100 * c(alpha),
                trim = TRUE,
                scientific = FALSE,
                digits = 3
            ),
            "%"
        ), "Est.")
    }

    # cerca
    # ci <- data.frame(ci)
    return(ci)
}

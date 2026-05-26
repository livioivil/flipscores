pkgname <- "flipscores"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('flipscores')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("anova.flipscores")
### * anova.flipscores

flush(stderr()); flush(stdout())

### Name: anova.flipscores
### Title: anova.flipscores
### Aliases: anova.flipscores

### ** Examples

set.seed(1)
dt=data.frame(X=scale(rnorm(50)),
   Z=factor(rep(LETTERS[1:3],length.out=50)))
dt$Y=rpois(n=nrow(dt),lambda=exp(dt$X*(dt$Z=="C")))
mod0=flipscores(Y~Z+X,data=dt,family="poisson")
summary(mod0)
anova(mod0)

mod1=flipscores(Y~Z*X,data=dt,family="poisson")
summary(mod1)
anova(mod0,model1 = mod1)




cleanEx()
nameEx("compute_scores")
### * compute_scores

flush(stderr()); flush(stdout())

### Name: compute_scores
### Title: compute_scores
### Aliases: compute_scores

### ** Examples

set.seed(1)
Z=rnorm(20)
X=Z+rnorm(20)
Y=rpois(n=20,lambda=exp(Z+X))
mod0=glm(Y~Z,family="poisson")
X=data.frame(X=X)
scr0=compute_scores(model0 = mod0, model1 = X)
head(scr0)



cleanEx()
nameEx("confint.flipscores")
### * confint.flipscores

flush(stderr()); flush(stdout())

### Name: confint.flipscores
### Title: Confidence intervals for flipscores
### Aliases: confint.flipscores

### ** Examples

set.seed(1)
dt <- data.frame(X=scale(rnorm(50)),
   Z=factor(rep(LETTERS[1:3],length.out=50)))
dt$Y=rpois(n=nrow(dt),lambda=exp(dt$X*(dt$Z=="C")))
mod0 <- flipscores(Y~Z+X,data=dt,family="poisson")
summary(mod0)
confint(mod0, parm="X")

xx <- rnorm(20)
zz <- rnorm(20, 0.2 * xx)
yy <- rnorm(20, 1 + 2 * xx - zz)
flps <- make_flips(20, 1000)
mod1 <- flipscores(yy ~ xx + zz)
summary(mod1)
confint(mod1, flips = flps)



cleanEx()
nameEx("custom_contrasts")
### * custom_contrasts

flush(stderr()); flush(stdout())

### Name: custom_contrasts
### Title: Build custom contrasts for flip-score tests
### Aliases: custom_contrasts apply_custom_contrasts

### ** Examples

set.seed(1)
dat <- data.frame(
  treatment = factor(rep(c("Control", "Drug"), each = 30)),
  sex = factor(rep(c("F", "M"), times = 30)),
  age_group = factor(rep(c("Young", "Old", "Middle"), length.out = 60)),
  severity = runif(60, 0, 1)
)
eta <- with(dat,
  -1 + 0.8 * (treatment == "Drug") + 0.4 * (sex == "M") +
    0.6 * (age_group == "Old") - 0.3 * severity
)
dat$y <- rbinom(nrow(dat), 1, plogis(eta))

fit <- glm(y ~ treatment * sex + age_group + severity,
           data = dat, family = binomial)

# Simple treatment contrast for one factor
custom_contrasts(fit, ~ treatment)

# Pairwise treatment contrasts inside each sex
custom_contrasts(fit, pairwise ~ treatment | sex)

# Interaction-cell contrasts, similar to emmeans pairwise comparisons
cc <- custom_contrasts(fit, pairwise ~ treatment:age_group | sex)
apply_custom_contrasts(cc)



cleanEx()
nameEx("flipscores-package")
### * flipscores-package

flush(stderr()); flush(stdout())

### Name: flipscores-package
### Title: flipscores: Robust Score Testing in GLMs, by Sign-Flip
###   Contributions
### Aliases: flipscores-package
### Keywords: internal

### ** Examples

set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp(dt$X))
mod=flipscores(Y~Z+X,data=dt,family="poisson",x=TRUE)
summary(mod)

# Anova test
anova(mod)
# or
mod0=flipscores(Y~Z,data=dt,family="poisson",x=TRUE)
anova(mod0,mod)
# and
mod0=flipscores(Y~X,data=dt,family="poisson")
anova(mod0,mod)




cleanEx()
nameEx("flipscores")
### * flipscores

flush(stderr()); flush(stdout())

### Name: flipscores
### Title: Robust testing in GLMs, by sign-flipping score contributions
### Aliases: flipscores

### ** Examples

set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp((dt$Z=="C") + 2*dt$X))
mod=flipscores(Y~Z+X,data=dt,family="poisson",n_flips=1000)
summary(mod)

# Equivalent to:
model=glm(Y~Z+X,data=dt,family="poisson")
mod2=flipscores(model)
summary(mod2)




cleanEx()
nameEx("flipscores_contrasts")
### * flipscores_contrasts

flush(stderr()); flush(stdout())

### Name: flipscores_contrasts
### Title: Flip-score tests for custom contrasts
### Aliases: flipscores_contrasts

### ** Examples

set.seed(1)
dat <- data.frame(
  treatment = factor(rep(c("Control", "Drug"), each = 30)),
  sex = factor(rep(c("F", "M"), times = 30)),
  age_group = factor(rep(c("Young", "Old", "Middle"), length.out = 60)),
  severity = runif(60, 0, 1)
)
eta <- with(dat,
  -1 + 0.8 * (treatment == "Drug") + 0.4 * (sex == "M") +
    0.6 * (age_group == "Old") - 0.3 * severity
)
dat$y <- rbinom(nrow(dat), 1, plogis(eta))

fit <- glm(y ~ treatment * sex + age_group + severity,
           data = dat, family = binomial, x = TRUE)

# Test the treatment effect separately for F and M.
flipscores_contrasts(fit, pairwise ~ treatment | sex,
                     n_flips = 100, seed = 1)

# Test selected interaction-cell contrasts. Increase n_flips in real analyses.
flipscores_contrasts(fit, ~ treatment:age_group | sex,
                     n_flips = 100, seed = 1)



cleanEx()
nameEx("gcor")
### * gcor

flush(stderr()); flush(stdout())

### Name: gcor
### Title: Compute Generalized Partial Correlations for GLM terms
### Aliases: gcor

### ** Examples

set.seed(1)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rpois(n=20,lambda=exp(dt$Z=="C"))
mod=flipscores(Y~Z+X,data=dt,family="poisson",n_flips=1000)
summary(mod)

# Compute generalized partial correlations for all terms
(results <- gcor(mod))

# Compute for specific terms only
gcor(mod, terms = c("X", "ZC"))

gcor(mod, terms = c("X", "ZC"),normalize=TRUE)


gcor(mod, intercept_too=TRUE, normalize=TRUE)
set.seed(123)
dt=data.frame(X=rnorm(20),
   Z=factor(rep(LETTERS[1:3],length.out=20)))
dt$Y=rbinom(n=20,prob=plogis((dt$Z=="C")*2),size=1)
mod=flipscores(Y~Z+X,data=dt,family="binomial",n_flips=1000)
summary(mod)

(results <- gcor(mod,normalize=TRUE))
# Compute for specific terms only
gcor(mod, terms = c("X", "ZC"),normalize=TRUE)





cleanEx()
nameEx("make_flips")
### * make_flips

flush(stderr()); flush(stdout())

### Name: make_flips
### Title: It creates a 'n_flips'x'n_obs' matrix of random +1 and -1.  The
###   first row is made by ones (i.e. the observed test statistic is
###   computed)
### Aliases: make_flips

### ** Examples

# example code
make_flips(n_obs=10,n_flips=20)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

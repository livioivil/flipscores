rm(list=ls())
setwd("C:\\Users\\livio\\Documents\\github")

library(devtools)
document("flipscores")

install("flipscores")
library(flipscores)

set.seed(1)
x=(rep(0:1,20))
D=data.frame(y=rbinom(40,1,.25+x*.5),x=x,
             z=rnorm(40),id=rep(1:20,each=2))

mod_par=glm(y~x*z,data=D,family = binomial)
summary(mod_par)
mod_par

mod=flipscores(y~x*z,data=D,family = binomial,score_type = "effect")
summary(mod)

print(mod)

mod0=flipscores(y~x+z,data=D,family = binomial,score_type = "effect")
# undebug(compute_scores)
# compute_scores(mod0,mod)

anova(mod0,mod)
stats:::anova.glm(mod0,mod)


#########
D2=D
D2$z=as.factor(round(D2$z))
mod=flipscores(y~x+z,data=D2,family = binomial,score_type = "effect")
summary(mod)


stats:::anova.glm(mod)
anova(mod)

## clustered units
mod=flipscores(y~x*z,data=D,family = binomial,id=D$id)
summary(mod)


##### poisson
D$y=rpois(40,exp(1+.25*x))
summary(glm(y~x*z,data=D,family = poisson))
mod=flipscores(y~x*z,data=D,family = poisson)
summary(mod)

#### negative binomial
summary(MASS::glm.nb(y~x*z,data=D))
mod=flipscores(y~x*z,data=D,family="negbinom")
summary(mod)


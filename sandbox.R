# rm(list=ls())
# setwd("C:\\Users\\livio\\Documents\\github")
# 
# library(devtools)
# document("flip")
# install("flip")
# library(flip)
# 
# res=flip(cbind(1:5,5:1))
# npc(res,subsets = c(list(1),list(2)))
# npc(res,subsets = c(list(1:2),list(2:1)))
# npc(res,subsets = c(list(1:2)))

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

modz=flipscores(y~z,data=D,family = binomial)
summary(modz)
anova(modz)
modz=update(modz,x=TRUE)
anova(modz)


mod0=flipscores(y~x+z,data=D,family = binomial,score_type = "effect")


modI=glm(y~1,data=D,family = binomial)
modI=flipscores(y~1,data=D,family = binomial,score_type = "effect")

# undebug(compute_scores)
# compute_scores(mod0,mod)

anova(mod,mod0)
stats:::anova.glm(mod,mod0)


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


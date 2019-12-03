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


mod=flipscores(y~x*z,data=D,family = binomial)
summary(mod)

anova(mod)

modz=flipscores(y~z,data=D,family = binomial)
summary(modz)
anova(modz)
modz=update(modz,x=TRUE)
anova(modz)


mod0=flipscores(y~x+z,data=D,family = binomial,score_type = "effect")
anova(mod0)

modI=glm(y~1,data=D,family = binomial)
modI=flipscores(y~1,data=D,family = binomial,score_type = "effect")

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
mod=flipscores(y~x*z,data=D,family = poisson)
summary(mod)
anova(mod)

#### negative binomial
summary(MASS::glm.nb(y~x*z,data=D))
mod=flipscores(y~x*z,data=D,family="negbinom")
summary(mod)
# anova(mod)

################### correlation- multivariate
# D$y=rpois(40,exp(1+.25*x+latent))
# D$y2=rpois(40,exp(1+.25*x+latent))
# D$y=rpois(40,exp(latent))
# D$y2=rpois(40,exp(latent))
# 
# mod=flipscores(y~x*z,data=D,family = poisson,score_type = "o")
# mod2=flipscores(y2~x*z,data=D,family = poisson,score_type = "o")
# 
# cor(mod$scores[,2],mod2$scores[,2])
# plot(mod$scores[,2],mod2$scores[,2])
# 
# 
# ############
# n=200
# x=(rep(0:1,length.out=n))
# D=data.frame(y=rbinom(n,1,.25+x*.5),x=x,
#              z=rnorm(n),id=rep(1:(n/2),each=2))
# 
# res=c()
# for(i in 1:500)
#   {
#   # set.seed(i)
#   D$latent=rnorm(n)
#   D$y=rpois(n,exp(D$latent))
#   D$y2=rpois(n,exp(D$latent))
#   
#   scores=compute_scores(model0 = glm(y~x+z,data=D,family = poisson),
#                         model1 = glm(y~x*z,data=D,family = poisson),
#                         score_type = "o")
#   scores2=compute_scores(model0 = glm(y2~x+z,data=D,family = poisson),
#                          model1 = glm(y2~x*z,data=D,family = poisson),
#                          score_type = "o")
#   res=rbind(res,c(cor(scores,scores2),sum(scores),sum(scores2)))
# }  
# 
# # sim <- function(){
# #   # set.seed(i)
# #   D$latent=rnorm(n)
# #   D$y=rpois(n,exp(D$latent))
# #   D$y2=rpois(n,exp(D$latent))
# #   
# #   scores=compute_scores(model0 = glm(y~x+z,data=D,family = poisson),
# #                         model1 = glm(y~x*z,data=D,family = poisson),
# #                         score_type = "e")
# #   scores2=compute_scores(model0 = glm(y2~x+z,data=D,family = poisson),
# #                          model1 = glm(y2~x*z,data=D,family = poisson),
# #                          score_type = "e")
# #   c(cor(scores,scores2),sum(scores),sum(scores2))
# # }
# # res=plyr::laply(1:10,sim)
# # sim()
# 
# summary(res)
# hist(res[,1])
# cor(res[,-1])
# mean(res[,1])
# 
# 
# X=matrix(rnorm(12),4,3)
# sv=svd(X)
# sv$u%*%t(sv$u)
# t(sv$u)%*%sv$u

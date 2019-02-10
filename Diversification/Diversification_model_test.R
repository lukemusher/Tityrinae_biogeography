#Some of the scripts here were initially provided by Dr. Frank Burbrink during his course at AMNH, "Applied Phylogenetics"
#DDD scripts were provided by Dr. Fabien Condamine

setwd("")
library(RPANDA)
library(laser)
library(geiger)
library(DDD)
library(TESS)

#Read your tree
read.tree("Tree.tre")->tree

#Plot your lineages through time
ltt.plot(tree, log="y", col="blue", lwd=3, xlim=c(-40, 0))

#simulate some trees under constant b&d and plot their LTTs with your LTT
replicate(10000,rbdtree(birth=bd.km(tree), death=0.0001, Tmax =max(branching.times(tree)), BIRTH = NULL,
                      DEATH = NULL, eps = 1e-6), F)->simmed
class(simmed)<-"multiPhylo"

mltt.plot(simmed, log="y", legend=F)
ltt.lines(tree, col="darkred", lwd=3)

#one-tailed t-test gamma-test to test if div rates decline through time
gamStat(branching.times(tree))

##Let's look at the 9 Models of diversification from Morlon Fig.1; 
#http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000493

Model1 <- fit_coal_cst(tree, tau0=1.e-3, gamma=-1, cst.rate=T, N0=50)
Model2<-fit_coal_cst(tree, tau0=1.e-3, gamma=-1, cst.rate=FALSE, N0=50)
Model3<-fit_coal_var(tree, lamb0=0.34, alpha=-0.001, mu0=0.12, beta=0, N0=50, cst.lamb = T, cst.mu =T)
Model4a<-fit_coal_var(tree, lamb0=0.34, alpha=-0.001, mu0=0.12, beta=0, N0=50, cst.lamb = F, cst.mu = T)
Model4b<-fit_coal_var(tree, lamb0=0.34, alpha=-0.001, mu0=0.12, beta=0, N0=50, cst.lamb = T, cst.mu = F)
Model4c<-fit_coal_var(tree, lamb0=0.34, alpha=-0.001, mu0=0.12, beta=0, N0=50, cst.lamb = F, cst.mu = F, fix.eps=T)
Model4d<-fit_coal_var(tree, lamb0=0.34, alpha=-0.001, mu0=0.12, beta=0, N0=50, cst.lamb = F, cst.mu = F, fix.eps=F)
Model5<-fit_coal_var(tree, lamb0=0.34, alpha=-0.001, mu0=0.12, beta=0, N0=50, cst.lamb = T, mu.0=T)
Model6<-fit_coal_var(tree, lamb0=0.34, alpha=-0.001, mu0=0.12, beta=0, N0=50, cst.lamb = F, mu.0=T)

weights<-aicw(c(Model1$aicc, Model2$aicc, Model3$aicc, Model4a$aicc, Model4b$aicc, Model4c$aicc, Model4d$aicc, Model5$aicc, Model6$aicc))

#environmental model selection in RPANDA: 
#environmental model selection: first condition on crown
#first models without envi variable
?fit_bd

#pure birth constant rate
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y}
lamb_par<-c(0.3)
mu_par<-c(0)

BD_model1<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, fix.mu =T, cond = "crown")
BD_model1

#birth-death constant rates
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3)
mu_par<-c(0.1)

BD_model2<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, cond = "crown")
BD_model2

#pure birth exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0)

BD_model3<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, fix.mu = T, cond = "crown", dt = 1e-3)
BD_model3

#speciation exponential dependence on time with constant extinction
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1)

BD_model4<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = T, cond = "crown")
BD_model4

#birth-death exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]*exp(y[2]*t)}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1,0.01)

BD_model5<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, cond = "crown")
BD_model5

#pure birth exponential dependence on temperature
env_data<-("Zachos_climate.csv")
f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0)

Env_model1<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.mu = T,fix.mu = T, expo.lamb = T, expo.mu = T)
Env_model1
plot_fit_env(Env_model1, env, 22.90199861)

#birth exponential dependence on temperature constant extinction

f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1)

Env_model2<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.mu = T,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model2
#plot_fit_env(Env_model2, env, 22.90199861)

#exponential dependence of speciation and extinction on temperature

f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1,0.01)

Env_model3<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.lamb = F, cst.mu = F,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model3
plot_fit_env(Env_model2, env, 22.90199861)

AICs<-c(BD_model1$aicc, BD_model2$aicc, BD_model3$aicc, BD_model4$aicc, BD_model5$aicc, Env_model1$aicc, Env_model2$aicc, Env_model3$aicc)
aic.w(AICs)

#environmental model selection: first condition on stem
#first models without envi variable

#pure birth constant rate
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y}
lamb_par<-c(0.3)
mu_par<-c(0)

BD_model1<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, fix.mu =T, cond = "stem")
BD_model1

#birth-death constant rates
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3)
mu_par<-c(0.1)

BD_model2<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, cond = "stem")
BD_model2

#pure birth exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0)

BD_model3<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, fix.mu = T, cond = "stem", dt = 1e-3)
BD_model3

#speciation exponential dependence on time with constant extinction
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1)

BD_model4<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = T, cond = "stem")
BD_model4

#birth-death exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]*exp(y[2]*t)}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1,0.01)

BD_model5<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, cond = "stem")
BD_model5

#pure birth exponential dependence on temperature
env_data<-("Zachos_climate.csv")
f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0)

Env_model1<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3, cst.mu = T,fix.mu = T, expo.lamb = T, expo.mu = T)
Env_model1
plot_fit_env(Env_model1, env, 22.90199861)

#birth exponential dependence on temperature constant extinction

f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1)

Env_model2<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3, cst.mu = T,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model2
#plot_fit_env(Env_model2, env, 22.90199861)

#exponential dependence of speciation and extinction on temperature

f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.3, 0.01)
mu_par<-c(0.1,0.01)

Env_model3<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3, cst.lamb = F, cst.mu = F,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model3
plot_fit_env(Env_model2, env, 22.90199861)

AICs<-c(BD_model1$aicc, BD_model2$aicc, BD_model3$aicc, BD_model4$aicc, BD_model5$aicc, Env_model1$aicc, Env_model2$aicc, Env_model3$aicc)
aic.w(AICs)

#DDD analyses based on condamine et al., 2018, biol letters.

library(TreePar)
library(DDD)
library(phytools)

Tityrinae<-read.tree("Tree.tre")

phylo<-list(Tityrinae)

names<-c("Tityrinae")

nbclades<-length(phylo)
phyloi<-Tityrinae
brtsi<-getx(phyloi)
initi<-1+length(brtsi)+8
resi<-10*(1+length(brtsi)+8)

#use "empirical priors" for initial parameters. Using 0.4 which is roughly the mean speciation rate from TESS
#and use 0.3 which is roughly the extinction rate from TESS 
#Because this is an initial parameter used to help optimize the model, precise values are not necessary
print("1. Linear dependence of speciation rate without extinction")
DDD_1<-dd_ML(brtsi, ddmodel=1, initparsopt=c(0.4,initi), idparsopt=c(1,3), idparsfix=c(2), parsfix=c(0), res=resi, missnumspec=8, cond=0, btorph=1, soc=1, maxiter=2000)
print(DDD_1)

#model does not easily converge (different results for different ititial parameters) 
#we therefore run the model with our initial parameters, 
#and then again with our new parameter estimates for lambda and mu
#we also use a very high maximum number of iterations of ten million
print("2. Linear dependence of speciation rate with extinction")
DDD_2<-dd_ML(brtsi, ddmodel=1, initparsopt=c(0.4,0.3,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000)
print(DDD_2)

#likelihood of first iteration is -108.9423, lambda = 0.3400318, mu = 0.1185399, K = 64321

DDD_2<-dd_ML(brtsi, ddmodel=1, initparsopt=c(0.34,0.12,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000)
print(DDD_2)

#new likelihood is marginally improved -108.9421, lambda = 0.3394826, mu = 0.1180327, K = 18618

#We run another similar two-step optimization for model 3
print("3. Exponential dependence of speciation rate with extinction")
DDD_3<-dd_ML(brtsi, ddmodel=2, initparsopt=c(0.4,0.3,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000)
print(DDD_3)

#likelihood of first iteration is -109.0863, lambda = 0.4502944, mu = 0.193024, K = infinity

DDD_3<-dd_ML(brtsi, ddmodel=2, initparsopt=c(0.45,0.19,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000)
print(DDD_3)

#again marginally improved likelihood of -109.0132, lambda = 0.3900675, mu = 0.1492525, K = infinity

#model 4 does not improve after a single additional optimization, and we therefore also
#change the initial K value for the second optimization
print("4. Linear dependence of extinction rate")
DDD_4<-dd_ML(brtsi, ddmodel=3, initparsopt=c(0.4,0.3,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000)
print(DDD_4)

#likelihood of first iteration is -108.9422, lambda = 0.3410106, mu = 0.1199411, K = 356419.9

DDD_4<-dd_ML(brtsi, ddmodel=3, initparsopt=c(0.34,0.12,200000), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000)
print(DDD_4)

#results converge on stable optimization so original model is used

DDD_4<-dd_ML(brtsi, ddmodel=3, initparsopt=c(0.4,0.3,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000)
print(DDD_4)

print("5. exponential dependence of extinction rate")
DDD_5<-dd_ML(brtsi, ddmodel=4, initparsopt=c(0.4,0.3,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000, optimmethod = "simplex")
print(DDD_5)

#initial optimization cannot be optimized. Need different set of parapeters. Use

DDD_5<-dd_ML(brtsi, ddmodel=4, initparsopt=c(0.36,0.09,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000, optimmethod = "simplex")
print(DDD_5)

#initial optimization results in log likelihood of -109.3438, lambda = 0.2765708, mu = 4.2e-07, K = 8.50

DDD_5<-dd_ML(brtsi, ddmodel=4, initparsopt=c(0.36,0.09,2000), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=10000000, optimmethod = "simplex")
print(DDD_5)

#marginally improved, but one more optimization because this model is much more complex

bd_1<-bd_ML(brts = brtsi, initparsopt = c(0.4,0.3), missnumspec = 8, cond = 0, btorph = 1, soc = 1, maxiter = 10000000)
print(bd_1)

#results table with AICc and AICc weights:

results<-matrix(NA,6,8)


#colnames(results)<-c("Model","df","logL","AICc","Lambda","Mu","K","r")
#results[,1]<-c("DDL","DDL+E","DDX+E","DD+EL","DD+EX","DDL+EL")
colnames(results)<-c("Model","df","logL","AICc","AICw","Lambda","Mu","K")
results[,1]<-c("DDL","DDL+E","DDX+E","DD+EL","DD+EX", "BD")

#df
results[1,2]<-round(as.numeric(DDD_1[5]))
results[2,2]<-round(as.numeric(DDD_2[5]))
results[3,2]<-round(as.numeric(DDD_3[5]))
results[4,2]<-round(as.numeric(DDD_4[5]))
results[5,2]<-round(as.numeric(DDD_5[5]))
results[6,2]<-round(2)

#logL
results[1,3]<-round(as.numeric(DDD_1[4]),4)
results[2,3]<-round(as.numeric(DDD_2[4]),4)
results[3,3]<-round(as.numeric(DDD_3[4]),4)
results[4,3]<-round(as.numeric(DDD_4[4]),4)
results[5,3]<-round(as.numeric(DDD_5[4]),4)
results[6,3]<-round(as.numeric(bd_1[5]),4)

#AICc
results[1,4]<-round((2*(-round(as.numeric(DDD_1[4]),4))+2*round(as.numeric(DDD_1[5]))+(2*round(as.numeric(DDD_1[5]))*(round(as.numeric(DDD_1[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_1[5]))-1)),3)
results[2,4]<-round((2*(-round(as.numeric(DDD_2[4]),4))+2*round(as.numeric(DDD_2[5]))+(2*round(as.numeric(DDD_2[5]))*(round(as.numeric(DDD_2[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_2[5]))-1)),3)
results[3,4]<-round((2*(-round(as.numeric(DDD_3[4]),4))+2*round(as.numeric(DDD_3[5]))+(2*round(as.numeric(DDD_3[5]))*(round(as.numeric(DDD_3[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_3[5]))-1)),3)
results[4,4]<-round((2*(-round(as.numeric(DDD_4[4]),4))+2*round(as.numeric(DDD_4[5]))+(2*round(as.numeric(DDD_4[5]))*(round(as.numeric(DDD_4[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_4[5]))-1)),3)
results[5,4]<-round((2*(-round(as.numeric(DDD_5[4]),4))+2*round(as.numeric(DDD_5[5]))+(2*round(as.numeric(DDD_5[5]))*(round(as.numeric(DDD_5[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_5[5]))-1)),3)
results[6,4]<-round((2*(-round(as.numeric(bd_1[5]),4))+2*round(as.numeric(bd_1[6]))+(2*round(as.numeric(bd_1[6]))*(round(as.numeric(bd_1[6]))+1))/(Ntip(phyloi)-round(as.numeric(bd_1[6]))-1)),3)

#Lambda
results[1,6]<-round(as.numeric(DDD_1[1]),4)
results[2,6]<-round(as.numeric(DDD_2[1]),4)
results[3,6]<-round(as.numeric(DDD_3[1]),4)
results[4,6]<-round(as.numeric(DDD_4[1]),4)
results[5,6]<-round(as.numeric(DDD_5[1]),4)
results[6,6]<-round(as.numeric(bd_1[1]),4)

#Mu
results[2,7]<-round(as.numeric(DDD_2[2]),5)
results[3,7]<-round(as.numeric(DDD_3[2]),5)
results[4,7]<-round(as.numeric(DDD_4[2]),5)
results[5,7]<-round(as.numeric(DDD_5[2]),5)
results[6,7]<-round(as.numeric(bd_1[2]),5)

#K
results[1,8]<-round(as.numeric(DDD_1[3]),2)
results[2,8]<-round(as.numeric(DDD_2[3]),2)
results[3,8]<-round(as.numeric(DDD_3[3]),2)
results[4,8]<-round(as.numeric(DDD_4[3]),2)
results[5,8]<-round(as.numeric(DDD_5[3]),2)


weights<-aic.w(as.numeric(results[,4]))

results[,5]<-round(weights,3)

results


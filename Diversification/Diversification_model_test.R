#Some of the scripts here were initially provided by Dr. Frank Burbrink during his course at AMNH, "Applied Phylogenetics"
#DDD scripts were provided by Dr. Fabien Condamine

setwd("")
library(RPANDA)
library(laser)

#Read your tree
read.tree("Tree.tre")->tree

#one-tailed t-test gamma-test to test if div rates decline through time
gamStat(branching.times(tree))

#environmental model selection: first condition on stem
#pure birth constant rate
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y}
lamb_par<-c(0.3)
mu_par<-c(0)

BD_model1<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, fix.mu =T, cond = "stem")
BD_model1
#plot_fit_bd(BD_model1, tot_time = 22.90199861)
#plot_dtt(BD_model1, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#birth-death constant rates
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3)
mu_par<-c(0.1)

BD_model2<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, cond = "stem")
BD_model2
#plot_fit_bd(BD_model2, tot_time = 22.90199861)
#plot_dtt(BD_model2, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#pure birth exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0)

BD_model3<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, fix.mu = T, cond = "stem", dt = 1e-3)
BD_model3

plot_fit_bd(BD_model3, tot_time = 22.90199861)
plot_dtt(BD_model3, tot_time = 22.90199861, N0 = 56)
ltt.lines(tree, col="darkred", lwd=3)

#speciation exponential dependence on time with constant extinction
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1)

BD_model4<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = T, cond = "stem")
BD_model4
#plot_fit_bd(BD_model4, tot_time = 22.90199861)
#plot_dtt(BD_model4, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#exponential dependence of extinction on time with constant speciation
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y[1]*exp(y[2]*t)}
lamb_par<-c(0.3)
mu_par<-c(0.1, 0.01)

BD_model5<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = F, cond = "stem")
BD_model5
#plot_fit_bd(BD_model5, tot_time = 22.90199861)
#plot_dtt(BD_model5, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#birth-death exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]*exp(y[2]*t)}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1,0.01)

BD_model6<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, cond = "stem")
BD_model6
#plot_fit_bd(BD_model6, tot_time = 22.90199861)
#plot_dtt(BD_model6, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#pure birth exponential dependence on temperature
#env<-read.csv("Zachos_climate.csv")
data("InfTemp")
env<-InfTemp
f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0)

Env_model1<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3, cst.mu = T,fix.mu = T, expo.lamb = T, expo.mu = T)
Env_model1
#plot_fit_env(Env_model1, env, 22.90199861)

#birth exponential dependence on temperature constant extinction
f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1)

Env_model2<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3, cst.mu = T,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model2
#plot_fit_env(Env_model2, env, 22.90199861)

#exponential dependence of extinction on temperature constant speciation
f.lamb <-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.3)
mu_par<-c(0.1, 0.01)

Env_model3<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3, cst.mu = F, cst.lamb = T, fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model3
#plot_fit_env(Env_model3, env, 22.90199861)

#exponential dependence of speciation and extinction on temperature

f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1,0.01)

Env_model4<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3, cst.lamb = F, cst.mu = F,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model4
#plot_fit_env(Env_model4, env, 22.90199861)

AICs<-c(BD_model1$aicc, BD_model2$aicc, BD_model3$aicc, 
        BD_model4$aicc, BD_model5$aicc, BD_model6$aicc, Env_model1$aicc, 
        Env_model2$aicc, Env_model3$aicc,Env_model4$aicc)

weights2<-aic.w(AICs)
weights2

#environmental model selection: next condition on crown
#first models without envi variable

#pure birth constant rate
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y}
lamb_par<-c(0.3)
mu_par<-c(0)

BD_model1<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, fix.mu =T, cond = "crown")
BD_model1
#plot_fit_bd(BD_model1, tot_time = 22.90199861)
plot_dtt(BD_model1, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#birth-death constant rates
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3)
mu_par<-c(0.1)

BD_model2<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = T, cond = "crown")
BD_model2
#plot_fit_bd(BD_model2, tot_time = 22.90199861)
#plot_dtt(BD_model2, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#pure birth exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0)

BD_model3<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, fix.mu = T, cond = "crown", dt = 1e-3)
BD_model3

#plot_fit_bd(BD_model3, tot_time = 22.90199861)
#plot_dtt(BD_model3, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#speciation exponential dependence on time with constant extinction
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1)

BD_model4<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = T, cond = "crown")
BD_model4
#plot_fit_bd(BD_model4, tot_time = 22.90199861)
#plot_dtt(BD_model4, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#exponential dependence of extinction on time with constant speciation
f.lamb<-function(t,y){y[1]}
f.mu<-function(t,y){y[1]*exp(y[2]*t)}
lamb_par<-c(0.3)
mu_par<-c(0.1, 0.01)

BD_model5<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=T, cst.mu = F, cond = "crown")
BD_model5
#plot_fit_bd(BD_model5, tot_time = 22.90199861)
#plot_dtt(BD_model5, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#birth-death exponential dependence on time
f.lamb<-function(t,y){y[1]*exp(y[2]*t)}
f.mu<-function(t,y){y[1]*exp(y[2]*t)}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1,0.01)

BD_model6<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cst.lamb=F, cst.mu = F, cond = "crown")
BD_model6
#plot_fit_bd(BD_model6, tot_time = 22.90199861)
#plot_dtt(BD_model6, tot_time = 22.90199861, N0 = 56)
#ltt.lines(tree, col="darkred", lwd=3)

#pure birth exponential dependence on temperature
f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0)

Env_model1<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.mu = T,fix.mu = T, expo.lamb = T, expo.mu = T)
Env_model1
#plot_fit_env(Env_model1, env, 22.90199861)

#birth exponential dependence on temperature constant extinction
f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1)

Env_model2<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.mu = T,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model2

#exponential dependence of extinction on temperature constant speciation
f.lamb <-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.3)
mu_par<-c(0.1, 0.01)

Env_model3<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.mu = F, cst.lamb = T, fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model3

#exponential dependence of speciation and extinction on temperature

f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.3, -0.1)
mu_par<-c(0.1,0.01)

Env_model4<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.lamb = F, cst.mu = F,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model4

f.lamb <-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.35, -0.06)
mu_par<-c(0.0006,0.36)

Env_model4<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3, cst.lamb = F, cst.mu = F,fix.mu = F, expo.lamb = T, expo.mu = T)
Env_model4

#plot_fit_env(Env_model4, env, 22.90199861)

AICs<-c(BD_model1$aicc, BD_model2$aicc, BD_model3$aicc, 
        BD_model4$aicc, BD_model5$aicc, BD_model6$aicc, Env_model1$aicc, 
        Env_model2$aicc, Env_model3$aicc,Env_model4$aicc)
aic.w(AICs)


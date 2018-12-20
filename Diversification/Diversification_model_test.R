#Some of the scripts here were initially provided by Dr. Frank Burbrink during his course at AMNH, "Applied Phylogenetics"
#DDD scripts were provided by Dr. Fabien Condamine

setwd("")
library(RPANDA)
library(laser)
library(geiger)
library(DDD)
library(TESS)

read.tree("Tree.tre")->tree

ltt.plot(tree, log="y", col="blue", lwd=3, xlim=c(-40, 0))

replicate(10000,rbdtree(birth=bd.km(tree), death=0.0001, Tmax =max(branching.times(tree)), BIRTH = NULL,
                      DEATH = NULL, eps = 1e-6), F)->simmed
class(simmed)<-"multiPhylo"

mltt.plot(simmed, log="y", legend=F)
ltt.lines(tree, col="darkred", lwd=3)

gamStat(branching.times(tree))

##Let's look at the 9 Models of diversification from Morlon Fig.1; 
#http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000493

Model3<-fit_coal_var(tree, lamb0=0.15, alpha=-0.001, mu0=0.01, beta=0, N0=70, cst.lamb = T, cst.mu =T)
Model4a<-fit_coal_var(tree, lamb0=0.15, alpha=-0.001, mu0=0.01, beta=0, N0=70, cst.lamb = F, cst.mu = T)
Model4b<-fit_coal_var(tree, lamb0=0.15, alpha=-0.001, mu0=0.01, beta=0, N0=70, cst.lamb = T, cst.mu = F)
Model4c<-fit_coal_var(tree, lamb0=0.15, alpha=-0.001, mu0=0.01, beta=0, N0=70, cst.lamb = F, cst.mu = F, fix.eps=T)
Model4d<-fit_coal_var(tree, lamb0=0.15, alpha=-0.001, mu0=0.01, beta=0, N0=70, cst.lamb = F, cst.mu = F, fix.eps=F)
Model5<-fit_coal_var(tree, lamb0=0.15, alpha=-0.001, mu0=0.01, beta=0, N0=70, cst.lamb = T, mu.0=T)
Model6<-fit_coal_var(tree, lamb0=0.15, alpha=-0.001, mu0=0.01, beta=0, N0=70, cst.lamb = F, mu.0=T)

as.table(c(Model1$aicc, Model2$aicc, Model3$aicc, Model4a$aicc, Model4b$aicc, Model4c$aicc, Model4d$aicc, Model5$aicc, Model6$aicc))->table
names(table)<-c("Model1", "Model2", "Model3", "Model4a", "Model4b", "Model4c", "Model4d", "Model5", "Model6")
table
min(table)-table->tree_table

#compare model tables
tree_table


#environmental model selection: first condition on crown

f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par<-c(0.09)
mu_par<-c()

BD_model<-fit_bd(phylo = tree, tot_time = 22.90199861, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cond = "crown")

env<-read.csv("Zachos_climate.csv")
f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)} #create lambda parameter function where lambda is exponential function of envi variable
f.mu<-function(t,x,y){0}
lamb_par<-c(0.15, 0.01)
mu_par<-c()

Env_model<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 22.90199861, lamb_par = lamb_par, mu_par = mu_par, cond = "crown", dt = 1e-3)
plot_fit_env(Env_model, env, 22.90199861)

#condition on stem
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par<-c(0.09)
mu_par<-c()

BD_model_stem<-fit_bd(phylo = tree, tot_time = 25.34, f.lamb = f.lamb, f.mu = f.mu, lamb_par = lamb_par, mu_par = mu_par, f = 0.85, cond = "stem")

f.lamb <-function(t,x,y){y[1] * exp(y[2] * x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.15, 0.01)
mu_par<-c()

Env_model_stem<-fit_env(tree, env_data = env, f.lamb = f.lamb, f.mu = f.mu, f=0.85, tot_time = 25.34, lamb_par = lamb_par, mu_par = mu_par, cond = "stem", dt = 1e-3)

plot_fit_env(Env_model_stem, env, 25.34)

as.table(c(BD_model$aicc, BD_model_stem$aicc, Env_model$aicc, Env_model_stem$aicc))->table
names(table)<-c("BD_crown", "BD_stem", "Env_crown", "Env_stem")
table
min(table)-table->tree_table

tree_table

#DDD analyses

setwd("/Users/lmusher/AMNH/Pachyramphus/biogeography_2/DDD/")
library(TreePar)
library(DDD)

Tityrinae<-read.tree("Tree.tre")

phylo<-list(Tityrinae)

names<-c("Tityrinae")

nbclades<-length(phylo)

finalTityrinaeDDD<-list()
resTityrinaeDDD<-c()

phyloi<-Tityrinae
brtsi<-getx(phyloi)
initi<-1+length(brtsi)+8
resi<-10*(1+length(brtsi)+8)

print("1. Linear dependence of speciation rate without extinction")
DDD_1<-dd_ML(brtsi, ddmodel=1, initparsopt=c(0.04,initi), idparsopt=c(1,3), idparsfix=c(2), parsfix=c(0), res=resi, missnumspec=8, cond=0, btorph=1, soc=1, maxiter=500)
print(DDD_1)

print("2. Linear dependence of speciation rate with extinction")
DDD_2<-dd_ML(brtsi, ddmodel=1, initparsopt=c(0.04,0.03,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=500)
print(DDD_2)

print("3. Exponential dependence of speciation rate with extinction")
DDD_3<-dd_ML(brtsi, ddmodel=2, initparsopt=c(0.04,0.03,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=500)
print(DDD_3)

print("4. Linear dependence of extinction rate")
DDD_4<-dd_ML(brtsi, ddmodel=3, initparsopt=c(0.04,0.03,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=500)
print(DDD_4)

print("5. Exponential dependence of extinction rate")
DDD_5<-dd_ML(brtsi, ddmodel=4.2, initparsopt=c(0.04,0.03,initi), res=resi, missnumspec= 8, cond=0, btorph=1, soc=1, maxiter=500)
print(DDD_5)

bd_1<-bd_ML(brts = brtsi, initparsopt = c(0.04,0.03), missnumspec = 8, cond = 0, btorph = 1, soc = 1, maxiter = 500)
print(bd_1)

results<-matrix(NA,6,7)

#colnames(results)<-c("Model","df","logL","AICc","Lambda","Mu","K","r")
#results[,1]<-c("DDL","DDL+E","DDX+E","DD+EL","DD+EX","DDL+EL")
colnames(results)<-c("Model","df","logL","AICc","Lambda","Mu","K")
results[,1]<-c("DDL","DDL+E","DDX+E","DD+EL","DD+EX", "LDSE, BD")

#df
results[1,2]<-round(as.numeric(DDD_1[5]))
results[2,2]<-round(as.numeric(DDD_2[5]))
results[3,2]<-round(as.numeric(DDD_3[5]))
results[4,2]<-round(as.numeric(DDD_4[5]))
results[5,2]<-round(as.numeric(DDD_5[5]))
results[6,2]<-round(as.numeric(bd_1[5]))

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
results[1,5]<-round(as.numeric(DDD_1[1]),4)
results[2,5]<-round(as.numeric(DDD_2[1]),4)
results[3,5]<-round(as.numeric(DDD_3[1]),4)
results[4,5]<-round(as.numeric(DDD_4[1]),4)
results[5,5]<-round(as.numeric(DDD_5[1]),4)
results[6,5]<-round(as.numeric(bd_1[1]),4)

#Mu
results[2,6]<-round(as.numeric(DDD_2[2]),5)
results[3,6]<-round(as.numeric(DDD_3[2]),5)
results[4,6]<-round(as.numeric(DDD_4[2]),5)
results[5,6]<-round(as.numeric(DDD_5[2]),5)
results[6,6]<-round(as.numeric(bd_1[2]),5)

#K
results[1,7]<-round(as.numeric(DDD_1[3]),2)
results[2,7]<-round(as.numeric(DDD_2[3]),2)
results[3,7]<-round(as.numeric(DDD_3[3]),2)
results[4,7]<-round(as.numeric(DDD_4[3]),2)
results[5,7]<-round(as.numeric(DDD_5[3]),2)

results

save(results,file="finalTityrinae_stem_DDD.Rdata")
write.table(results,file="Results_Tityrinae_stem_DDD.txt",quote=FALSE,sep="\t",row.names=FALSE)

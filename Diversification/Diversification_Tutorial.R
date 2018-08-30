#scripts here provided by Frank Burbrink during his course, "Applied Phylogenetics"

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


#This script is written by Lukas J Musher. Please cite:
#Musher, L.J., M. Ferreira, A. Auerbach, J. Mckay, and J. Cracraft. 
#Why is Amazonia a source of biodiversity? Climate mediated dispersal and synchronous 
#speciation across the Andes in an avian group (Tityrinae) (In revision) Proceedings of the Royal Society B

require("ape")
library(gplots)
library(phytools)
library(ape)
library(phangorn)
library(ips)

#must have exactly same genes and trees--if there are 20 genes there must be 20 genes 
#and named appropriately so you know you have the right tree with each gene
#ideally you will have uces and gts named exactly the same e.g. uce-10.nexus uce-10.tre
#the easiest way to do this is make a list of the genes, and use a for loop to change
# the names of the RAxML outputs to just the locus name
#e.g. in linux/unix: for i in `cat names.genes`; do mv RAxML*$i* $i.tre; done
#runs slowly on datasets with lots of tips
#Start fresh R to maximize efficiency

setwd("example_data/")

#If you have your alignments in fasta format, you can quickly convert them in R with the following
#files1 <- list.files(path="/Users/lmusher/Desktop/mafft-nexus-clean-54taxa-95per/clock_test/", pattern="*.fasta", full.names=F, recursive=FALSE)
#convert_fasta_to_nexus<-function(files){
#  for (i in 1:length(files)){
#    x<-read.dna(files[i], format = "fasta")
#    write.nexus.data(x,file=paste(files[i],".nexus",sep=""))
#  }
#}
#convert_fasta_to_nexus(files1)

files <- list.files(path="/Users/lmusher/AMNH/Pachyramphus/biogeography_2/mafft-nexus-clean-54taxa-95per/", pattern="*.nexus", full.names=F, recursive=FALSE)
trees <- list.files(path="/Users/lmusher/AMNH/Pachyramphus/biogeography_2/mafft-nexus-clean-54taxa-95per/",pattern="*.tre",full.names=F,recursive=FALSE)

spec_tree<- read.tree("../RAxML_bestTree.best.tre")
plot(spec_tree)
outgroup<-c("acachl")
remove.taxa<-c("A","B") #any taxa in the tree that you don't want involved in the likelihood estimations, for example taxa that might overinfluence the likelihoodds due to missing data etc

pis<-function(x){
  len<-length(x[[1]])
  x<-read.nexus.data(x)
  site<-c()
  counts=0
  for(j in 1:length(x[[1]])){
    for(i in 1:length(x)){
      site[i]<-(x[[i]][j])
    }
    if(!("?" %in% unique(site))&!("-" %in% unique(site))){
      if (length(unique(site))>=2){
        counts=counts+1
      }
    }
    else if(!("?" %in% unique(site))){
      if (length(unique(site))>=3){
        counts=counts+1
      }
    }
    else if(!("-" %in% unique(site))){
      if (length(unique(site))>=3){
        counts=counts+1
      }
    }
    else if (length(unique(site))>4){
      #print(unique(site))
      counts=counts+1
    }
    else counts=counts
  }
  return(counts)
}

#how many bp sites in the alignment
site_count<-function(x){
  x<-read.nexus.data(x)
  len<-length(x[[1]])
  return(len)
}

#lambda is a smoothing parameter for estimating divergence times using penalized likelihood
#lambda is correlated with "clock-likeness" higher values indicate closer to clock i.e. even/lined up tips
find_lambda_best<-function(x,tree,remove.taxa){
  x<-read.nexus.data(x)
  x<-as.DNAbin(x)
  X<-as.phyDat(x)
  t<-c()
  counts=0
  for(i in 1:length(tree$tip.label)){
    if(!(tree$tip.label[i] %in% labels(X))){
      counts=counts+1
      t[counts]<-tree$tip.label[i]
    }
  }
  t<-t[!is.na(t)]
  tree1<-drop.tip(tree,tree$tip.label[na.omit(match(t,tree$tip.label))])
  remove.taxa_new<-c()
  for(i in 1:length(remove.taxa)){
    if(remove.taxa[i] %in% labels(X)){
      remove.taxa_new[i]<-remove.taxa[i]
    }
  }
  remove.taxa_new<-remove.taxa_new[!is.na(remove.taxa_new)]
  tree1<-drop.tip(tree1,tree$tip.label[na.omit(match(remove.taxa_new,tree$tip.label))])
  fit<-pml(tree1, X, model="GTRG")   #do likelihood estimation of the distance tree l<-10^(-1:6)
  l<-c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1.0,10,100)
  cv<-sapply(l, function(x) sum(attr(chronopl(midpoint(tree1), lambda=x, CV=T), "D2")))
  df<-as.data.frame(cbind(l,cv))
  l.min<-df$l[df$cv==min(df$cv)]
  returnValue(l.min)
}

#Use likelihood ratio tests to estimate likelihood ratio of clocklike to not clocklike model
is.clocklike<-function(dna,tree,lambda,outgroup,remove.taxa){
  x<-read.nexus.data(dna)
  x<-as.DNAbin(x)
  X<-as.phyDat(x)
  t<-c()
  counts=0
  for(i in 1:length(tree$tip.label)){
    if(!(tree$tip.label[i] %in% labels(X))){
      counts=counts+1
      t[counts]<-tree$tip.label[i]
    }
  }
  t<-t[!is.na(t)]
  tree1<-drop.tip(tree,tree$tip.label[na.omit(match(t,tree$tip.label))])
  #plot(tree1)
  outgroup_new<-c()
  for(i in 1:length(outgroup)){
    if(outgroup[i] %in% labels(X)){
      outgroup_new[i]<-outgroup[i]
    }
  }
  outgroup_new<-outgroup_new[!is.na(outgroup_new)]
  remove.taxa_new<-c()
  for(i in 1:length(remove.taxa)){
    if(remove.taxa[i] %in% labels(X)){
      remove.taxa_new[i]<-remove.taxa[i]
    }
  }
  remove.taxa_new<-remove.taxa_new[!is.na(remove.taxa_new)]
  tree1<-drop.tip(tree1,tree$tip.label[na.omit(match(remove.taxa_new,tree$tip.label))])
  plot(tree1)
  fit<-pml(tree1, X, model="GTRG")   #do likelihood estimation of tree 
  #par(mfrow=c(1,2))
  #plot(midpoint(fit$tree), main=dna)
  chr<-chronopl(phy=midpoint(tree1), lambda=lambda, age.min=1, CV=T) #create chronogram
  fit2<-pml(chr, X, model="GTRG")   #do likelihood estimation of chronogram 
  LR=(2*(fit$logLik-fit2$logLik))
  df=length(fit$tree$tip.label)-2
  p=pchisq(LR,df,lower.tail = F)
  if (p<0.05){
    plot(midpoint(chr),main="Reject clocklike model")
  }
  if (p>=0.05){
    plot(midpoint(chr), main="Accept clocklike model")
  }
  return(LR)
}

#create the dataframe
Make_df_lambda<-function(x,trees,outgroup,remove.taxa){
  loc<-c()
  lambda<-c() #create empty vector
  n.taxa<-c()
  log10lambda<-c()
  var.sites<-c()
  var.frac<-c()
  seq.len<-c()
  df<-c()
  LR<-c()
  p<-c()
  tst<-c()
  Tre<-c()
  for(i in 1:length(x)){ #for each locus in the directory
    loc[i]<-x[i]
    tree<-read.tree(trees[i])
    Tre[i]<-trees[i]
    remove.taxa_new<-c()
    for(j in 1:length(remove.taxa)){
      if(remove.taxa[j] %in% tree$tip.label){
        remove.taxa_new[j]<-remove.taxa[j]
      }
    }
    remove.taxa_new<-remove.taxa_new[!is.na(remove.taxa_new)]
    tree<-drop.tip(tree,tree$tip.label[na.omit(match(remove.taxa_new,tree$tip.label))])
    print(paste("Locus",i,": ",loc[i],sep=""))
    lambda[i]<-find_lambda_best(x[i],tree,remove.taxa) #the ith element is the result of finding the var.frac of the ith locus in the directory
    n.taxa[i]<-length(read.nexus.data(x[i]))
    log10lambda[i]<-log10(lambda[i])
    var.sites[i]<-pis(x[i])
    seq.len[i]<-site_count(x[i])
    var.frac[i]<-var.sites[i]/seq.len[i]
    cl<-is.clocklike(x[i],tree,lambda[i],outgroup,remove.taxa)
    deg.free<-length(read.nexus.data(x[i]))-2
    df[i]<-deg.free
    pval<-(pchisq(cl,df,lower.tail = FALSE))
    p[i]<-pval
    t<-c()
    counts<-0
    for(k in 1:length(spec_tree$tip.label)){
      if(!(spec_tree$tip.label[k] %in% tree$tip.label)){
        counts=counts+1
        t[counts]<-spec_tree$tip.label[k]
      }
    }
    t<-t[!is.na(t)]
    tree1<-drop.tip(spec_tree,spec_tree$tip.label[na.omit(match(t,spec_tree$tip.label))])
    if (pval<0.05){
      print(paste("Rejected m/c Clock for",loc[i]))
      tst[i]<-"Reject_clock"
    }
    if (pval>=0.05){
      print(paste("Accepted m/c Clock for",loc[i]))
      tst[i]<-"Accept_clock"
    }
    LR[i]<-cl
    print(paste("Smoothing par (lambda) =",lambda[i]))
    print(paste("log10 of Smoothing par (lambda) =",log10lambda[i]))
    print(paste("LR for clock model =",LR[i]))
    data.f<-as.data.frame(cbind(loc,Tre,n.taxa,seq.len,var.sites,var.frac,lambda,log10lambda,LR,p,df,tst))
    write.csv(data.f,file = "../timetree-4376-3May18.csv")
  }
  return(data.f)
}


uce_table<-Make_df_lambda(files,trees,outgroup,remove.taxa)

###############################
##########ANALYSES#############
###############################
library(BayesianFirstAid)
library(rjags)

par(mfrow=c(1,1))
tab<-read.csv("../timetree-4376-3May18.csv")
attach(tab)
names(tab)
hist(LR)
min(log10lambda)
#look at correlations. Might take a while under bayesian as below
outputCor1<-bayes.cor.test(var.sites, LR)
outputCor2<-bayes.cor.test(var.sites, log10lambda)
outputCor3<-bayes.cor.test(log10lambda,LR)
outputCor4<-bayes.cor.test(var.frac,LR)
plot(outputCor1)
plot(outputCor2)
plot(outputCor3)
plot(outputCor4)

range(LR)
mean(LR)
sd(LR)

mod<-glm(LR~var.sites+var.frac+log10lambda)
cooksd<-cooks.distance(mod)
summary(mod)
plot(cooksd, pch="+", cex=1, main="Influential Obs by Cooks distance")
abline(h=3*mean(cooksd, na.rm=T), col="red")
points(x=ifelse(cooksd>3*mean(cooksd, na.rm=T), names(cooksd),""), col="red")
#text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>3*mean(cooksd, na.rm=T), names(cooksd),""),col="red")

influential <- as.numeric(names(cooksd)[cooksd > 3*mean(cooksd, na.rm=T)])
tab[influential, ]
x<-c(1:1160)
tab2<-tab[! x %in% influential,]

write.csv(tab2,file="../tab2.csv")

tab3<-tab[influential,]
write.csv(tab3,file="../tab3.csv")

mod2<-glm(tab2$LR~tab2$var.sites+tab2$var.frac+tab2$log10lambda)
summary(mod2)
par(mfrow=c(5,5), mar=c(0.1,0.1,0.8,0.1))

###plot some outliers
for (i in tab3[15:39,]$Tre){
  tree<-read.tree(i)
  if("acachl" %in% tree$tip.label){
    tree<-midpoint(ladderize(tree), "acachl")
  }
  if(!("acachl" %in% tree$tip.label)){
    tree<-midpoint(ladderize(tree))
  }
  plot.phylo(tree, show.tip.label = F, main=paste("IS:",tab[tab$Tre==i,5], "LR:",round(tab[tab$Tre==i,9])))
}

#plot low and high lambda trees
par(mfrow=c(5,5), mar=c(0.1,0.1,0.8,0.1))
for (i in tab2[1:25]$Tre){
  tree<-read.tree(i)
  if(tab[tab$Tre==i,8]==-6){
    plot.phylo(ladderize(midpoint(tree)), show.tip.label = F, main=paste("LR=",round(tab[tab$Tre==i,9])))
  }
}

par(mfrow=c(10,3), mar=c(0.1,0.1,0.8,0.1))
for (i in tab$Tre){
  tree<-read.tree(i)
  if(tab[tab$Tre==i,8]==2){
    plot.phylo(ladderize(midpoint(tree)), show.tip.label = F, main=paste("LR=",round(tab[tab$Tre==i,9])))
  }
}

tab2<-read.csv("../tab2.csv")

#plot low LR trees
par(mfrow=c(5,5), mar=c(0.1,0.1,0.8,0.1))

for (i in tab2[1:25,]$Tre){
  tree<-read.tree(i)
  plot.phylo(ladderize(midpoint(tree)), show.tip.label = F, main=paste("IS:",tab[tab$Tre==i,5], "LR:",round(tab[tab$Tre==i,9])))
}

par(mfrow=c(10,3), mar=c(0.1,0.1,0.8,0.1))
for (i in tab$Tre){
  tree<-read.tree(i)
  if(tab[tab$Tre==i,9]>16000){
    plot.phylo(ladderize(midpoint(tree)), show.tip.label = F, main=paste("LR=",round(tab[tab$Tre==i,9])))
  }
}

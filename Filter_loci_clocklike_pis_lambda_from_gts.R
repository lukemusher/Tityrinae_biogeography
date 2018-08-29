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

setwd("/Users/lmusher/AMNH/Pachyramphus/biogeography_2/mafft-nexus-clean-54taxa-95per/")

#files1 <- list.files(path="/Users/lmusher/Desktop/mafft-nexus-clean-54taxa-95per/clock_test/", pattern="*.fasta", full.names=F, recursive=FALSE)
#convert_phylip_to_nexus<-function(files){
#  for (i in 1:length(files)){
    x<-read.dna(files[i], format = "fasta")
    write.nexus.data(x,file=paste(files[i],".nexus",sep=""))
  }
#}
#convert_phylip_to_nexus(files1)

files <- list.files(path="/Users/lmusher/AMNH/Pachyramphus/biogeography_2/mafft-nexus-clean-54taxa-95per/", pattern="*.nexus", full.names=F, recursive=FALSE)
trees <- list.files(path="/Users/lmusher/AMNH/Pachyramphus/biogeography_2/mafft-nexus-clean-54taxa-95per/",pattern="*.tre",full.names=F,recursive=FALSE)

spec_tree<- read.tree("../RAxML_bestTree.best.tre")
plot(spec_tree)
outgroup<-c("acachl")#"melund","nesnot"
toepads<-c("H")#HUMAN","ALLIG","ANOCA","CHEMY"
#is.monophyletic(read.tree("uce-10.tre"),outgroup)
#plot(read.tree("uce-10.tre"))
#t<-read.tree("1001.tre")
#t$tip.label
#spec_tree$tip.label
#setdiff(t$tip.label,spec_tree$tip.label)

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

#t<-find_lambda_best(files[1],read.tree(trees[1]))

#plot(log(l),cv)
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
  #tree1<-drop.tip(tree1,tree$tip.label[na.omit(match(outgroup_new,tree$tip.label))])
  remove.taxa_new<-c()
  for(i in 1:length(remove.taxa)){
    if(remove.taxa[i] %in% labels(X)){
      remove.taxa_new[i]<-remove.taxa[i]
    }
  }
  remove.taxa_new<-remove.taxa_new[!is.na(remove.taxa_new)]
  tree1<-drop.tip(tree1,tree$tip.label[na.omit(match(remove.taxa_new,tree$tip.label))])
  plot(tree1)
  fit<-pml(tree1, X, model="GTRG")   #do likelihood estimation of the distance tree 
  #fit<-optim.pml(fit, rearrangement = "NNI",optGamma=TRUE, optInv=TRUE,model="GTR") #with NNI rearrangement
  #fit<-optim.pml(fit, rearrangement = "ratchet",optGamma=TRUE, optInv=TRUE,model="GTR") #with NNI rearrangement
  #par(mfrow=c(1,2))
  plot(midpoint(fit$tree), main=dna)
  chr<-chronopl(phy=midpoint(tree1), lambda=lambda, age.min=1, CV=T)
  fit2<-pml(chr, X, model="GTRG")   #do likelihood estimation of the distance tree 
  LR=(2*(fit$logLik-fit2$logLik))
  df=length(fit$tree$tip.label)-2
  #print(paste("-lnL(nonClock)=",fit$logLik))
  #print(paste("-lnL(Clock)=",fit2$logLik))
  #print(paste("LR=",(2*(fit2$logLik-fit$logLik))))
  #print(paste("p =",pchisq(LR,df)))
  p=pchisq(LR,df,lower.tail = F)
  if (p<0.05){
    plot(midpoint(chr),main="Reject clocklike model")
  }
  if (p>=0.05){
    plot(midpoint(chr), main="Accept clocklike model")
  }
  return(LR)
}

Make_df_lambda<-function(x,trees,outgroup,remove.taxa,spec_tree){
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
  #RF<-c()
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
    #RF[i]<-RF.dist(tree1,midpoint(tree))
    if (pval<0.05){
      print(paste("Rejected m/c Clock for",loc[i]))
      tst[i]<-"Reject_clock"
    }
    if (pval>=0.05){
      print(paste("Accepted m/c Clock for",loc[i]))
      tst[i]<-"Accept_clock"
    }
    LR[i]<-cl
    #print(paste("The RF Dist of GT to ST =",RF[i]))
    print(paste("Smoothing par (lambda) =",lambda[i]))
    print(paste("log10 of Smoothing par (lambda) =",log10lambda[i]))
    print(paste("LR for clock model =",LR[i]))
    #data.f<-as.data.frame(cbind(loc,Tre,n.taxa,seq.len,var.sites,var.frac,RF,lambda,log10lambda,LR,p,df,tst))
    data.f<-as.data.frame(cbind(loc,Tre,n.taxa,seq.len,var.sites,var.frac,lambda,log10lambda,LR,p,df,tst))
    write.csv(data.f,file = "../timetree-4376-3May18.csv")
  }
  return(data.f)
}


uce_table<-Make_df_lambda(files,trees,outgroup,toepads,spec_tree)

uce_table


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


par(mfrow=c(1,3))

#plot(var.sites,LR)
plot(tab2$var.sites,tab2$LR, xlim=c(0,250), ylim=c(0,26000), pch=19, col="black", xlab="Total informative sites", ylab="LR")
points(tab[influential,]$var.sites,tab[influential,]$LR, pch=19, col="gray")
abline(glm(LR~var.sites), col="gray",lwd=3)
abline(glm(tab2$LR~tab2$var.sites), col="red", lwd=3)

#plot(var.frac,LR)
plot(tab2$var.frac,tab2$LR, xlim=c(0,0.5), ylim=c(0,26000), pch=19, col="black", xlab="Proportion informative sites", ylab="")
points(tab[influential,]$var.frac,tab[influential,]$LR, pch=19, col="gray")
abline(glm(LR~var.frac),col="gray", lwd=3)
abline(glm(tab2$LR~tab2$var.frac), col="red", lwd=3)
#hist(tab[influential,]$var.frac)


#plot(log10lambda,LR)
plot(tab2$log10lambda,tab2$LR, xlim=c(-6,2), ylim=c(0,26000), pch=19, col="black",xlab="log10(lambda)", ylab="")
points(tab[influential,]$log10lambda,tab[influential,]$LR, pch=19, col="gray")
abline(lm(LR~log10lambda),col="gray", lwd=3)
abline(lm(tab2$LR~tab2$log10lambda), col="red", lwd=3)

#plot(var.frac,log10lambda)
plot(tab2$var.frac,tab2$log10lambda, ylim=c(-6,2), xlim=c(0,0.5), pch=1, col="darkgreen")
points(tab[influential,]$var.frac,tab[influential,]$log10lambda, pch="+", col="darkred")

summary(lm(LR~var.frac))
summary(lm(tab2$LR~tab2$var.frac))

summary(lm(LR~log10lambda))
summary(lm(tab2$LR~tab2$log10lambda))

hist(LR,xlim = c(0,25000), breaks = 100, col="darkred")
hist(tab2$LR, xlim = c(0,25000), breaks = 75, col="darkgreen", add=T)
hist(var.sites,xlim = c(0,250), ylim=c(0,80), breaks = 50, col="darkred")
hist(tab2$var.sites, xlim = c(0,250), ylim=c(0,80), breaks = 50, col="darkgreen", add=T)
hist(log10lambda,xlim = c(-6,2), ylim=c(0,400), breaks = 25, col="darkred")
hist(tab2$log10lambda, xlim = c(-6,2), ylim=c(0,400), breaks = 25, col="darkgreen", add=T)


tab0<-read.csv("../48-birds-uce-alignments-clean-95/48-birds-95-uce-metrics.csv")
tab1<-read.csv("../random_genes1/ranom1-cds-metrics.csv")
tab2<-read.csv("../random_genes2/ranom2-cds-metrics.csv")
tab3<-read.csv("../random_genes3/ranom3-cds-metrics.csv")
par(mfrow=c(1,1))
attach(tab1)
names(tab1)

dat.all<-read.csv("All_data.csv")
attach(dat.all)
names(dat.all)
#length(dat.all[[1]])
#LR<-as.numeric(LR)
class(LR)

hist(LR)
#RF<-as.numeric(RF)
#var.sites<-as.numeric(var.sites)
#pis.frac<-as.numeric(pis.frac)
#log10lambda<-as.numeric(log10lambda)
#seq.len<-as.numeric(seq.len)
par(mfrow=c(1,1))
plot(dataset,LR, ylab="Likelihood Ratio")
plot(dataset,scaled.rf, ylab="RF Distance")
plot(dataset,n.taxa)
plot(dataset,log10lambda, "")
plot(dataset,var.sites)
plot(dataset,pis.frac)
plot(dataset,seq.len)

plot(LR,RF)
abline(lm(RF~LR))
summary(lm(RF~LR))

plot(LR,log10lambda)
abline(lm(log10lambda~LR))
summary(lm(log10lambda~LR))

plot(pis.frac,LR)
abline(lm(LR~pis.frac))
summary(lm(LR~pis.frac))

plot(dat.all[7:12])

par(mfrow=c(2,1))

plot(var.sites,LR, main="Pre-filtering",xlab="Number of Variable Sites", ylab="Likelihood Ratio")
abline(lm(LR~var.sites))
summary(lm(LR~var.sites))
text(5000,400000,"Adj. R-squared = 0.4134")

plot(var.sites[9692:9776],LR[9692:9776],main="Post-filtering", xlab="Number of Variable Sites", ylab="Likelihood Ratio")
abline(lm(LR[9692:9776]~var.sites[9692:9776]))
summary(lm(LR[9692:9776]~var.sites[9692:9776]))
text(200,4000,"Adj. R-squared = 0.06043")


plot(var.sites,RF, xlab="Number of Variable Sites", ylab="RF genetree to species tree")
abline(lm(RF~var.sites))
summary(lm(RF~var.sites))
text(5200,75,"Adj. R-squared = 0.1106")

min(LR)

setwd("/Users/lmusher/AMNH/Fasta_Aves_genomes/48-birds-uces/")
dat.all<-read.csv("All_data.csv")

outputCor1<-bayes.cor.test(var.sites[dataset=="random1"],LR[dataset=="random1"], main="Random 1")
outputCor1
pdf(file="BayesCor-pisXLR-Rand1.pdf")
plot(outputCor1)
dev.off()

outputCor2<-bayes.cor.test(var.sites[dataset=="random2"],LR[dataset=="random2"], main="Random 2")
outputCor2
pdf(file="BayesCor-pisXLR-Rand2.pdf")
plot(outputCor2)
dev.off()

outputCor3<-bayes.cor.test(var.sites[dataset=="random3"],LR[dataset=="random3"], main="Random 3")
outputCor3
pdf(file="BayesCor-pisXLR-Rand3.pdf")
plot(outputCor3)
dev.off()

outputCor4<-bayes.cor.test(var.sites[dataset=="uce"],LR[dataset=="uce"], main="Ultraconserved Elements")
outputCor4
pdf(file="BayesCor-pisXLR-uce.pdf")
plot(outputCor4)
dev.off()

outputCor4_filt1<-bayes.cor.test(var.sites[dataset=="uce-95CI-pis-frac"],LR[dataset=="uce-95CI-pis-frac"], main="Ultraconserved Elements")
outputCor4_filt1
pdf(file="BayesCor-pisXLR-uce-filt-pis.pdf")
plot(outputCor4_filt1)
dev.off()

outputCor4_filt2<-bayes.cor.test(var.sites[dataset=="uce-95CI-pis-frac-5per-LR"],LR[dataset=="uce-95CI-pis-frac-5per-LR"], main="Ultraconserved Elements")
outputCor4_filt2
pdf(file="BayesCor-pisXLR-uce-filt-LR_pis.pdf")
plot(outputCor4_filt2)
dev.off()

outputCor5<-bayes.cor.test(var.sites[dataset=="random1"],scaled.rf[dataset=="random1"], main="Random 1")
outputCor5
pdf(file="BayesCor-pisXRF-rand1.pdf")
plot(outputCor5)
dev.off()

outputCor6<-bayes.cor.test(var.sites[dataset=="random2"],scaled.rf[dataset=="random2"], main="Random 2")
outputCor6
pdf(file="BayesCor-pisXRF-rand2.pdf")
plot(outputCor6)
dev.off()

outputCor7<-bayes.cor.test(var.sites[dataset=="random3"],scaled.rf[dataset=="random3"], main="Random 3")
outputCor7
pdf(file="BayesCor-pisXRF-rand3.pdf")
plot(outputCor7)
dev.off()

outputCor8<-bayes.cor.test(var.sites[dataset=="uce"],scaled.rf[dataset=="uce"], main="Ultraconserved Elements")
outputCor8
pdf(file="BayesCor-pisXRF-uce.pdf")
plot(outputCor8)
dev.off()

outputCor8_filt1<-bayes.cor.test(var.sites[dataset=="uce-95CI-pis-frac"],scaled.rf[dataset=="uce-95CI-pis-frac"], main="Ultraconserved Elements")
outputCor8_filt1
pdf(file="BayesCor-pisXRF-uce-filt-pis.pdf")
plot(outputCor8_filt1)
dev.off()

outputCor8_filt2<-bayes.cor.test(var.sites[dataset=="uce-95CI-pis-frac-5per-LR"],scaled.rf[dataset=="uce-95CI-pis-frac-5per-LR"], main="Ultraconserved Elements")
outputCor8_filt2
pdf(file="BayesCor-pisXRF-uce-filt-pis-LR.pdf")
plot(outputCor8_filt2)
dev.off()

outputCor9<-bayes.cor.test(LR[dataset=="random1"],scaled.rf[dataset=="random1"])
outputCor9
pdf(file="BayesCor-LRXRF-rand1.pdf")
plot(outputCor9)
dev.off()

outputCor10<-bayes.cor.test(LR[dataset=="random2"],scaled.rf[dataset=="random2"])
outputCor10
pdf(file="BayesCor-LRXRF-rand2.pdf")
plot(outputCor10)
dev.off()

outputCor11<-bayes.cor.test(LR[dataset=="random3"],scaled.rf[dataset=="random3"])
outputCor11
pdf(file="BayesCor-LRXRF-rand3.pdf")
plot(outputCor11)
dev.off()

outputCor12<-bayes.cor.test(LR[dataset=="uce"],scaled.rf[dataset=="uce"])
outputCor12
pdf(file="BayesCor-LRXRF-uce.pdf")
plot(outputCor12)
dev.off()

outputCor12_filt1<-bayes.cor.test(LR[dataset=="uce-95CI-pis-frac"],scaled.rf[dataset=="uce-95CI-pis-frac"])
outputCor12_filt1
pdf(file="BayesCor-LRXRF-uce-filt-pis.pdf")
plot(outputCor12_filt1)
dev.off()

outputCor12_filt2<-bayes.cor.test(LR[dataset=="uce-95CI-pis-frac-5per-LR"],scaled.rf[dataset=="uce-95CI-pis-frac-5per-LR"])
outputCor12_filt2
pdf(file="BayesCor-LRXRF-uce-filt-pis-LR.pdf")
plot(outputCor12_filt2)
dev.off()

outputCor13<-bayes.cor.test(LR[dataset=="random1"],log10lambda[dataset=="random1"])
outputCor13
pdf(file="BayesCor-LRXLambda-rand1.pdf")
plot(outputCor13)
dev.off()

outputCor14<-bayes.cor.test(LR[dataset=="uce-95CI-pis-frac-5per-LR"],log10lambda[dataset=="uce-95CI-pis-frac-5per-LR"])
outputCor14
pdf(file="BayesCor-LRXLambda-uce-filt-pis-LR.pdf")
plot(outputCor14)
dev.off()

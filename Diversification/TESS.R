library(ape)
library(phytools)
library(phangorn)
library(devtools)
library(TESS)
#We can run this with two expected changes, who knows? We come up with some expected survival prob after an extinction. 
#Using the expected survival probability, we compute the α and β parameters of
#the beta distribution. We set the value of β to be large, which focuses the prior
#density more tightly around the expected survival probability. Then, we compute
#α based on the expected survival probability and the specified β value

setwd("/Users/lmusher/AMNH/Pachyramphus/biogeography_2/TESS/")
read.tree("Tree.tre")->tyran

plot(tyran)

expectedSurvivalProbability <- 0.05

numExpectedRateChanges<-3
numExpectedMassExtinctions<-0

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}
tyran<-force.ultrametric(tyran)

tess.analysis(tyran,
              empiricalHyperPriors = TRUE,
              samplingProbability = 0.85,
              estimateNumberMassExtinctions = FALSE,
              MAX_ITERATIONS = 100000,
              dir = "comet_no_mass_extinctions_tyran_CC15")

output <- tess.process.output("comet_no_mass_extinctions_tyran_CC15",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)


par(mfrow=c(2,2))
tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times"),las=2)
par(mfrow=c(2,1))
tess.plot.output(output,
                 fig.types = c("net-diversification rates"
                               ),las=2, plot.tree=T, yaxt = )



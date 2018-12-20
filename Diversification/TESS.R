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

#evaluate convergence
tess.plot.singlechain.diagnostics(output = output, 
                                  parameters = 
                                    c("speciation rates","extinction rates"),
                                  las=2
)

#now let's divide the rates into 20 intevervals and correlate div rates with paleoclimate

#First create a new ouput where numIntevals = 20
output <- tess.process.output("comet_no_mass_extinctions_tyran_CC15",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions,
                              numIntervals = 20
                              )
#evaluate convergence
tess.plot.singlechain.diagnostics(output = output, 
                                  parameters = 
                                    c("speciation rates","extinction rates"),
                                  las=2
)

#write your the dates of your 20 intervals to a file
write.csv(output$intervals, file = "intervals.csv")

#right now you have posterior output for 20 intervals, but we want to summarize these as means 
#so we create list of the means of posterior net-div-rates of the 20 intevals

net.div.rates<c()

for (i in 1:20){
  net.div.rates[i]<-mean(output$`net-diversification rates`[,i])
}

#write the mean net-div-rates for each interval to a file
write.csv(net.div.rates, file = "net.div.rates.csv")

#I then take means of all the d18O values for each time interval and combine all the data into one "intervals.csv" by hand, 
#but this could easily be scripted if you have more than 20 intervals

tab<-read.csv("intervals.csv")

attach(tab)

names(tab)

plot(d18O,net.div.rates)
abline(lm(net.div.rates~d18O))

summary(lm(net.div.rates~d18O))

library(ape)
library(phytools)
library(phangorn)
library(TESS)

#We can run this with two expected changes, who knows? We come up with some expected survival prob after an extinction. 
#Using the expected survival probability, we compute the α and β parameters of
#the beta distribution. We set the value of β to be large, which focuses the prior
#density more tightly around the expected survival probability. Then, we compute
#α based on the expected survival probability and the specified β value

setwd("")
read.tree("Tree.tre")->tyran

plot(tyran, show.tip.label = F)

expectedSurvivalProbability <- 0.05

numExpectedRateChanges<-2
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
              estimateNumberMassExtinctions = F,
              BURNIN = 10000,
              MAX_ITERATIONS = 100000,
              verbose = T,
              dir = "comet_no_mass_extinctions_tyran_CC15")

output <- tess.process.output("comet_no_mass_extinctions_tyran_CC15",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions,
                              numIntervals = 50
                              )

#get 50 interval times
write.csv(output$intervals, file = "intervals50.csv")

#get div rates for each interval
net.div.rates<c()

for (i in 1:50){
  net.div.rates[i]<-mean(output$`net-diversification rates`[,i])
}


write.csv(net.div.rates, file = "net.div.rates50.csv")

par(mfrow=c(2,2))
tess.plot.singlechain.diagnostics(output = output, 
                                  parameters = 
                                    c("speciation rates","extinction rates"),
                                  las=2
)

#look at data
par(mfrow=c(2,3))
tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "speciation Bayes factors",
                               "extinction rates",
                               "extinction shift times",
                               "extinction Bayes factors"),las=2)
par(mfrow=c(1,1))
tess.plot.output(output,
                 fig.types = c("net-diversification rates"
                               ),las=2, plot.tree=T, yaxt = )


#glms for net-div and temperature

clim<-InfTemp #RPANDA Temp Data
attach(clim)
names(clim)

temp_means<-c()

#Get mean temp values for each TESS interval
for (i in 1:50){
  print(output$intervals[i])
  temp_means[i]<-mean(clim[Age<=output$intervals[i] 
                     & Age>=output$intervals[i+1],"Temperature"])
}
write.csv(temp_means, file = "temp2.csv")

#stitch together three csvs that have times, div data, and temp means and save as net.div.rates.csv (could do this in R)

tab<-read.csv("net.div.rates.csv")

attach(tab)

intervals_neg<-Interval*(-1) #make time negative
names(tab)

par(mfrow=c(1,1))
plot(d18O,net.div.rates)
abline(lm(net.div.rates~d18O))

#plot Figure S5
plot(Interval,net.div.rates, xlim = c(25,0), type="l", xlab="Time (mya)", ylab = "Net-diversification Rates")
#abline(lm(net.div.rates~Interval))
par(new=T)
plot(Interval,Temperature_means, xlim = c(25,0),type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", ylim = c(0,10), col="darkred")
axis(side=4, labels = T, col = "darkred")

#Get stats:
summary(lm(net.div.rates~d18O))

lm1<-glm(net.div.rates~intervals_neg)
lm2<-glm(net.div.rates~Temperature_means)
lm3<-glm(net.div.rates~Temperature_means*intervals_neg)
sum1<-summary(lm1)
sum2<-summary(lm2)
sum3<-summary(lm3)

Rsq1<-summary(lm(net.div.rates~intervals_neg))
Rsq2<-summary(lm(net.div.rates~d18O))
Rsq3<-summary(lm(net.div.rates~d18O*intervals_neg))

#summarize results in table S4

Models<-c("M1: Time", "M2: Temp", "M3: Time and Temp")
Rsqs<-c(Rsq1$adj.r.squared,Rsq2$adj.r.squared,Rsq3$adj.r.squared)
AICs<-c(sum1$aic,sum2$aic,sum3$aic)
weights<-aic.w(c(sum1$aic,sum2$aic,sum3$aic))
weights

results<-matrix(NA,3,4)

results[,1]<-Models
results[,2]<-Rsqs
results[,3]<-AICs
results[,4]<-weights

write.csv(results,file="glms.csv")


#Also to create figure3 for main text
output <- tess.process.output("comet_no_mass_extinctions_tyran_CC15",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions
                              
)

par(mfrow=c(1,1))
tess.plot.output(output,
                 fig.types = c("net-diversification rates"
                 ),las=2, plot.tree=T, yaxt = )
par(new=T)
plot(Age,Temperature, xlim = c(25,0),type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", ylim = c(-2,13), col="darkred")
axis(side=4, labels = T, col = "darkred")

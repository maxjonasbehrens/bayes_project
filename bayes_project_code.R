
# Load data ---------------------------------------------------------------

library(boot)
data(coal)
x <- 1851:1962
y <- tabulate(floor(coal[[1]]))
y <- y[x]

# Visualise the original data

data0 <- data.frame(c(1851:1962),cumsum(y))
names(data0) <- c("year","y")

require(ggplot2)
ggplot(data0,aes(year,y))+
  geom_line()+
  ylab("Cumulative number of disasters")+
  xlab("Year")+
  theme_bw(base_size = 16)+
  xlim(c(1850,1962))

ggsave("zero_cp_line.png", height = 6, width = 8)

# One Changepoint ---------------------------------------------------------

# Possible T values

T <- seq(2,112)

# Marginal Posterior T Distribution
pT <- rep(0,112)
pT[1] <- log(0)
for (t in T) {
  pT[t] <- log(114-t)*(-2-sum(y[t:112]))+ lgamma(sum(y[t:112])+2)+
          log(t)*(-2-sum(y[1:(t-1)])) + lgamma(sum(y[1:(t-1)])+2) # Formula of p(T|y), integrate everything else out
}

pT <- pT - max(pT)
pT <- exp(pT)
pT <- pT/sum(pT) # Standardise probabilities
sum(pT)

# MAP - 1 Changepoint
map1 <- which(pT==max(pT))+1850 # 1892

# Sample from T 
sampleT <- sample(x = 1:112, replace = TRUE, size = 10000, prob = pT)

# Sample from Lambda 1
sampleL1 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL1[sample] <- rgamma(1,2+sum(y[1:(sampleT[sample]-1)]),sampleT[sample])
}

# Sample from Lambda 2
sampleL2 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL2[sample] <- rgamma(1,2+sum(y[sampleT[sample]:112]),114-sampleT[sample])
}

# Concat the samples
postSam <- matrix(c(sampleT,sampleL1,sampleL2),nrow = 10000,ncol = 3)

# Visualise marginal likelihood functions
require(RColorBrewer)
require(MASS)

layout(matrix(c(1,2,3,4,5,6,7,8,9),nrow=3,byrow=T))
par(mar=c(3,3,1,1))
mycol = c("#FFFFFF",colorRampPalette(brewer.pal(3, "Reds"))(250))

d <- density(postSam[,2], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 1", line=2, cex.lab=1)

d <- kde2d(postSam[,3],postSam[,2],lims=c(0,5,0,5),n=256)
image(d$x,d$y,d$z,xlab="",ylab="",col=mycol,useRaster=T)
title(ylab="Lambda 1", line=2, cex.lab=1)
title(xlab="Lambda 2", line=2, cex.lab=1)

d <- kde2d(postSam[,1],postSam[,2],lims=c(0,112,0,5),n=256)
image(d$x,d$y,d$z,xlab="",ylab="",col=mycol,useRaster=T)
title(ylab="Lambda 1", line=2, cex.lab=1)
title(xlab="T", line=2, cex.lab=1)

plot(NA,frame=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1))

d <- density(postSam[,3], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "",ylab = "")
title(xlab = "Lambda 2")

d <- kde2d(postSam[,1],postSam[,3],lims=c(0,112,0,5),n=256)
image(d$x,d$y,d$z,xlab="",ylab="",col=mycol,useRaster=T)
title(ylab="Lambda 2", line=2, cex.lab=1)
title(xlab="T", line=2, cex.lab=1)

plot(NA,frame=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1))
plot(NA,frame=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1))

hist(postSam[,1], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
title(xlab="T", line=2, cex.lab=1)

# Draw predictions from from predictive posterior distributions
pred <- matrix(0,nrow = 10000, ncol = 112)
for (sample in 1:10000) {
  for (year in 1:112) {
    if(year < postSam[sample,1]){
      pred[sample,year] <- rpois(1,postSam[sample,2])
    } else {
      pred[sample,year] <- rpois(1,postSam[sample,3])
    }
  }
}

# Cummulative Sum of predictions
cp <- t(apply(pred,1,cumsum))

# Derive the quantiles of sampled predictions
median <- rep(0,112)
low <- rep(0,112)
upper <- rep(0,112)

for (year in 1:112) {
  median[year] <- quantile(cp[,year], probs = c(0.5))
  low[year] <- quantile(cp[,year], probs = c(0.025))
  upper[year] <- quantile(cp[,year], probs = c(0.975))
}

# Visualise the predictions
par(mfrow=c(1,1))

data1 <- data.frame(c(1851:1962),cumsum(y), median, low, upper)
names(data1) <- c("year","y","median","low","upper")

require(ggplot2)
ggplot(data1)+
  geom_ribbon(aes(year,ymin=low,ymax=upper,colour="95% Prediction Interval"),alpha = 0.3, fill = 'green')+
  geom_line(aes(year,y,colour="Observations"))+
  geom_line(aes(year,median,colour="Prediction Median"))+
  geom_vline(xintercept = map1, linetype='dotted')+
  annotate("text",label="MAP(T1) = 1892", x=1905,y=25)+
  scale_color_manual("",breaks=c("95% Prediction Interval","Observations","Prediction Median"),values = c("green","black","red"))+
  xlab("Years")+
  ylab("Cummulative number of disasters")+
  theme_bw(base_size = 16)+
  theme(legend.position = "top")
  
ggsave("one_cp_predictions.png", height = 6, width = 8)

# Two Changepoints --------------------------------------------------------

# Possible values for T1 and T2
T <- seq(2,112)
Tgrid <- expand.grid(T,T) 
names(Tgrid) <- c("t1","t2")
Tgrid <- Tgrid[Tgrid$t1 < Tgrid$t2,]

# Marginal Posterior T Distribution
pT <- rep(0,nrow(Tgrid)+1)
pT[1] <- log(0)
for (i in 1:nrow(Tgrid)) {
  pT[i+1] <- log(114-Tgrid[i,2])*(-2-sum(y[Tgrid[i,2]:112]))+ lgamma(sum(y[Tgrid[i,2]:112])+2)+
            log(Tgrid[i,2]-Tgrid[i,1])*(-2-sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)]))+ lgamma(sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)])+2)+
            log(Tgrid[i,1])*(-2-sum(y[1:(Tgrid[i,1]-1)])) + lgamma(sum(y[1:(Tgrid[i,1]-1)])+2) # Formula of p(T|y), integrate everything else out
}

pT <- pT - max(pT)
pT <- exp(pT)
pT <- pT/sum(pT) # Standardise probabilities
sum(pT) # Check if probabilities sum up to 1

# MAP of two changepoints
map_idx <- which(pT==max(pT))
map2 <- Tgrid[map_idx,]+1850

# Concat possible T values
postT <- matrix(c(c(0,Tgrid$t1),c(0,Tgrid$t2),pT),nrow = length(pT),ncol = 3)

# Sample from T1 and T2
sampleT <- sample(x = 1:nrow(postT), replace = TRUE, size = 10000, prob = postT[,3])
sampleT1 <- rep(0,10000)
sampleT2 <- rep(0,10000)
for (i in 1:length(sampleT)) {
  sampleT1[i] <- postT[sampleT[i],1]
  sampleT2[i] <- postT[sampleT[i],2]
}

# Sample from Lambda 1
sampleL1 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL1[sample] <- rgamma(1,2+sum(y[1:(sampleT1[sample]-1)]),sampleT1[sample])
}

# Sample from Lambda 2
sampleL2 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL2[sample] <- rgamma(1,2+sum(y[sampleT1[sample]:(sampleT2[sample]-1)]),sampleT2[sample]-sampleT1[sample])
}

# Sample from Lambda 3
sampleL3 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL3[sample] <- rgamma(1,2+sum(y[sampleT2[sample]:112]),114-sampleT2[sample])
}

# Concat the samples
postSam <- matrix(c(sampleT1,sampleT2,sampleL1,sampleL2,sampleL3),nrow = 10000,ncol = 5)

# Visualise distributions
layout(matrix(c(1,2,3,4,5,6),nrow=2,byrow=T))
par(mar=c(3,3,1,1))

d <- density(postSam[,3], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 1", line=2, cex.lab=1.2)
d <- density(postSam[,4], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 2", line=2, cex.lab=1.2)
d <- density(postSam[,5], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 3", line=2, cex.lab=1.2)
hist(postSam[,1], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
title(xlab="T1", line=2, cex.lab=1.2)
hist(postSam[,2], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
title(xlab="T2", line=2, cex.lab=1.2)

# Draw predictions from the sample
pred <- matrix(0,nrow = 10000, ncol = 112)
for (sample in 1:10000) {
  for (year in 1:112) {
    if(year < postSam[sample,1]) {
      pred[sample,year] <- rpois(1,postSam[sample,3])
    } else if (year >= postSam[sample,1] & year < postSam[sample,2]) {
      pred[sample,year] <- rpois(1,postSam[sample,4])
    } else {
      pred[sample,year] <- rpois(1,postSam[sample,5])
    }
  }
}

# Cummulative Predictions
cp <- t(apply(pred,1,cumsum))

# Quantiles of predictions
median <- rep(0,112)
low <- rep(0,112)
upper <- rep(0,112)

for (year in 1:112) {
  median[year] <- quantile(cp[,year], probs = c(0.5))
  low[year] <- quantile(cp[,year], probs = c(0.025))
  upper[year] <- quantile(cp[,year], probs = c(0.975))
}

# Visualise the predictions
par(mfrow=c(1,1))

data2 <- data.frame(c(1851:1962),cumsum(y), median, low, upper)
names(data2) <- c("year","y","median","low","upper")

require(ggplot2)
ggplot(data2)+
  geom_ribbon(aes(year,ymin=low,ymax=upper,colour="95% Prediction Interval"),alpha = 0.3, fill = 'green')+
  geom_line(aes(year,y,colour="Observations"))+
  geom_line(aes(year,median,colour="Prediction Median"))+
  geom_vline(xintercept = unlist(map2[1]), linetype='dotted')+
  geom_vline(xintercept = unlist(map2[2]), linetype='dotted')+
  annotate("text",label="MAP(T1) = 1892", x=1905,y=25)+
  annotate("text",label="MAP(T2) = 1892", x=1955,y=25)+
  scale_color_manual("",breaks=c("95% Prediction Interval","Observations","Prediction Median"),values = c("green","black","red"))+
  xlab("Years")+
  ylab("Cummulative number of disasters")+
  theme_bw(base_size = 16)+
  theme(legend.position = "top")

ggsave("two_cp_predictions.png", height = 6, width = 8)

# 3 Changepoints ----------------------------------------------------------

# All possible values for T1, T2 and T3
T <- seq(2,112)
Tgrid <- expand.grid(T,T,T) 
names(Tgrid) <- c("t1","t2","t3")
Tgrid <- Tgrid[Tgrid$t1 < Tgrid$t2 & Tgrid$t2 < Tgrid$t3-2,]

# Marginal Posterior T Distribution
pT <- rep(0,nrow(Tgrid)+1)
pT[1] <- log(0)
for (i in 1:nrow(Tgrid)) {
  pT[i+1] <- log(114-Tgrid[i,3])*(-2-sum(y[Tgrid[i,3]:112]))+ lgamma(sum(y[Tgrid[i,3]:112])+2)+
    log(Tgrid[i,3]-Tgrid[i,2])*(-2-sum(y[Tgrid[i,2]:(Tgrid[i,3]-1)]))+ lgamma(sum(y[Tgrid[i,2]:(Tgrid[i,3]-1)])+2)+
    log(Tgrid[i,2]-Tgrid[i,1])*(-2-sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)]))+ lgamma(sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)])+2)+
    log(Tgrid[i,1])*(-2-sum(y[1:(Tgrid[i,1]-1)])) + lgamma(sum(y[1:(Tgrid[i,1]-1)])+2) 
}

pT <- pT - max(pT)
pT <- exp(pT)
pT <- pT/sum(pT) # Standardise probabilities
sum(pT)

# MAP - 3 Changepoints
map_idx <- which(pT==max(pT))
map3 <- Tgrid[map_idx,]+1850

# Concat possible T values
postT <- matrix(c(c(0,Tgrid$t1),c(0,Tgrid$t2),c(0,Tgrid$t3),pT),nrow = length(pT),ncol = 4)

# Sample from T1, T2 and T3
sampleT <- sample(x = 1:nrow(postT), replace = TRUE, size = 10000, prob = postT[,4])
sampleT1 <- rep(0,10000)
sampleT2 <- rep(0,10000)
sampleT3 <- rep(0,10000)
for (i in 1:length(sampleT)) {
  sampleT1[i] <- postT[sampleT[i],1]
  sampleT2[i] <- postT[sampleT[i],2]
  sampleT3[i] <- postT[sampleT[i],3]
}

# Sample from Lambda 1
sampleL1 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL1[sample] <- rgamma(1,2+sum(y[1:(sampleT1[sample]-1)]),sampleT1[sample])
}

# Sample from Lambda 2
sampleL2 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL2[sample] <- rgamma(1,2+sum(y[sampleT1[sample]:(sampleT2[sample]-1)]),sampleT2[sample]-sampleT1[sample])
}

# Sample from Lambda 3
sampleL3 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL3[sample] <- rgamma(1,2+sum(y[sampleT2[sample]:(sampleT3[sample]-1)]),sampleT3[sample]-sampleT2[sample])
}

# Sample from Lambda 4
sampleL4 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL4[sample] <- rgamma(1,2+sum(y[sampleT3[sample]:112]),114-sampleT3[sample])
}

# Concat all distributions
postSam <- matrix(c(sampleT1,sampleT2,sampleT3,sampleL1,sampleL2,sampleL3,sampleL4),nrow = 10000,ncol = 7)

# Visualise distributions
layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,byrow=T))
par(mar=c(3,3,1,1))

d <- density(postSam[,4], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 1", line=2, cex.lab=1.2)
d <- density(postSam[,5], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 2", line=2, cex.lab=1.2)
d <- density(postSam[,6], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 3", line=2, cex.lab=1.2)
d <- density(postSam[,7], from = 0, to = 5, n = 1000)
plot(d$x,d$y,type = 'l',xlab = "Lambda 1",ylab = "")
title(xlab="Lambda 4", line=2, cex.lab=1.2)
hist(postSam[,1], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
title(xlab="T1", line=2, cex.lab=1.2)
hist(postSam[,2], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
title(xlab="T2", line=2, cex.lab=1.2)
hist(postSam[,3], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
title(xlab="T3", line=2, cex.lab=1.2)

# Make predictions from drawn posterior samples
pred <- matrix(0,nrow = 10000, ncol = 112)
for (sample in 1:10000) {
  for (year in 1:112) {
    if(year < postSam[sample,1]) {
      pred[sample,year] <- rpois(1,postSam[sample,4])
    } else if (year >= postSam[sample,1] & year < postSam[sample,2]) {
      pred[sample,year] <- rpois(1,postSam[sample,5])
    } else if (year >= postSam[sample,2] & year < postSam[sample,3]) {
      pred[sample,year] <- rpois(1,postSam[sample,6])
    } else {
      pred[sample,year] <- rpois(1,postSam[sample,7])
    }
  }
}

# Calculate the cummulative predictions
cp <- t(apply(pred,1,cumsum))

# Derive the quantiles
median <- rep(0,112)
low <- rep(0,112)
upper <- rep(0,112)

for (year in 1:112) {
  median[year] <- quantile(cp[,year], probs = c(0.5))
  low[year] <- quantile(cp[,year], probs = c(0.025))
  upper[year] <- quantile(cp[,year], probs = c(0.975))
}

# Visualise the predictions
par(mfrow=c(1,1))

data3 <- data.frame(c(1851:1962),cumsum(y), median, low, upper)
names(data3) <- c("year","y","median","low","upper")

require(ggplot2)
ggplot(data3)+
  geom_ribbon(aes(year,ymin=low,ymax=upper,colour="95% Prediction Interval"),alpha = 0.3, fill = 'green')+
  geom_line(aes(year,y,colour="Observations"))+
  geom_line(aes(year,median,colour="Prediction Median"))+
  geom_vline(xintercept = unlist(map3[1]), linetype='dotted')+
  geom_vline(xintercept = unlist(map3[2]), linetype='dotted')+
  geom_vline(xintercept = unlist(map3[3]), linetype='dotted')+
  annotate("text",label="MAP(T1) = 1892", x=1905,y=50)+
  annotate("text",label="MAP(T2) = 1930", x=1935,y=25)+
  annotate("text",label="MAP(T2) = 1948", x=1955,y=50)+
  scale_color_manual("",breaks=c("95% Prediction Interval","Observations","Prediction Median"),values = c("green","black","red"))+
  xlab("Years")+
  ylab("Cummulative number of disasters")+
  theme_bw(base_size = 16)+
  theme(legend.position = "top")

ggsave("three_cp_predictions.png", height = 6, width = 8)

# 4 Changepoints - Test ----------------------------------------------------------

# All possible values for T1, T2, T3 and T4
T <- seq(2,111)
Tgrid <- expand.grid(T,T,T,T) 
names(Tgrid) <- c("t1","t2","t3","t4")
Tgrid <- Tgrid[Tgrid$t1 < Tgrid$t2 & Tgrid$t2 < Tgrid$t3-2 & Tgrid$t3+2 < Tgrid$t4,]

# Marginal Posterior T Distribution
pT <- rep(0,nrow(Tgrid)+1)
pT[1] <- log(0)
for (i in 1:nrow(Tgrid)) {
  pT[i+1] <- log(114-Tgrid[i,4])*(-2-sum(y[Tgrid[i,3]:112]))+ lgamma(sum(y[Tgrid[i,4]:112])+2)+
    log(Tgrid[i,4]-Tgrid[i,3])*(-2-sum(y[Tgrid[i,3]:(Tgrid[i,4]-1)]))+ lgamma(sum(y[Tgrid[i,3]:(Tgrid[i,4]-1)])+2)+
    log(Tgrid[i,3]-Tgrid[i,2])*(-2-sum(y[Tgrid[i,2]:(Tgrid[i,3]-1)]))+ lgamma(sum(y[Tgrid[i,2]:(Tgrid[i,3]-1)])+2)+
    log(Tgrid[i,2]-Tgrid[i,1])*(-2-sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)]))+ lgamma(sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)])+2)+
    log(Tgrid[i,1])*(-2-sum(y[1:(Tgrid[i,1]-1)])) + lgamma(sum(y[1:(Tgrid[i,1]-1)])+2) 
}

pT <- pT - max(pT)
pT <- exp(pT)
pT <- pT/sum(pT) # Standardise probabilities
sum(pT)

# MAP - 4 Changepoints
map_idx <- which(pT==max(pT))
map4 <- Tgrid[map_idx,]+1850

# Concat possible T values
postT <- matrix(c(c(0,Tgrid$t1),c(0,Tgrid$t2),c(0,Tgrid$t3),c(0,Tgrid$t4),pT),nrow = length(pT),ncol = 5)

# Sample from T1, T2 and T3
sampleT <- sample(x = 1:nrow(postT), replace = TRUE, size = 10000, prob = postT[,5])
sampleT1 <- rep(0,10000)
sampleT2 <- rep(0,10000)
sampleT3 <- rep(0,10000)
sampleT4 <- rep(0,10000)
for (i in 1:length(sampleT)) {
  sampleT1[i] <- postT[sampleT[i],1]
  sampleT2[i] <- postT[sampleT[i],2]
  sampleT3[i] <- postT[sampleT[i],3]
  sampleT4[i] <- postT[sampleT[i],4]
}

# Sample from Lambda 1
sampleL1 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL1[sample] <- rgamma(1,2+sum(y[1:(sampleT1[sample]-1)]),sampleT1[sample])
}

# Sample from Lambda 2
sampleL2 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL2[sample] <- rgamma(1,2+sum(y[sampleT1[sample]:(sampleT2[sample]-1)]),sampleT2[sample]-sampleT1[sample])
}

# Sample from Lambda 3
sampleL3 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL3[sample] <- rgamma(1,2+sum(y[sampleT2[sample]:(sampleT3[sample]-1)]),sampleT3[sample]-sampleT2[sample])
}

# Sample from Lambda 4
sampleL4 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL4[sample] <- rgamma(1,2+sum(y[sampleT3[sample]:(sampleT4[sample]-1)]),sampleT4[sample]-sampleT3[sample])
}

# Sample from Lambda 5
sampleL5 <- rep(0,10000)
for (sample in 1:10000) {
  sampleL5[sample] <- rgamma(1,2+sum(y[sampleT4[sample]:112]),114-sampleT4[sample])
}

# Concat all distributions
postSam <- matrix(c(sampleT1,sampleT2,sampleT3,sampleT4,sampleL1,sampleL2,sampleL3,sampleL4,sampleL5),nrow = 10000,ncol = 9)

# Visualise
hist(postSam[,1], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
hist(postSam[,2], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
hist(postSam[,3], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)
hist(postSam[,4], freq = FALSE, breaks = 112, xlim = c(0,112), main = NULL)

# Make predictions from drawn posterior samples
pred <- matrix(0,nrow = 10000, ncol = 112)
for (sample in 1:10000) {
  for (year in 1:112) {
    if(year < postSam[sample,1]) {
      pred[sample,year] <- rpois(1,postSam[sample,5])
    } else if (year >= postSam[sample,1] & year < postSam[sample,2]) {
      pred[sample,year] <- rpois(1,postSam[sample,6])
    } else if (year >= postSam[sample,2] & year < postSam[sample,3]) {
      pred[sample,year] <- rpois(1,postSam[sample,7])
    } else if (year >= postSam[sample,3] & year < postSam[sample,4]) {
      pred[sample,year] <- rpois(1,postSam[sample,8])
    } else {
      pred[sample,year] <- rpois(1,postSam[sample,9])
    }
  }
}

# Calculate the cummulative predictions
cp <- t(apply(pred,1,cumsum))

# Derive the quantiles
median <- rep(0,112)
low <- rep(0,112)
upper <- rep(0,112)

for (year in 1:112) {
  median[year] <- quantile(cp[,year], probs = c(0.5))
  low[year] <- quantile(cp[,year], probs = c(0.025))
  upper[year] <- quantile(cp[,year], probs = c(0.975))
}

par(mfrow=c(1,1))

data4 <- data.frame(c(1851:1962),cumsum(y), median, low, upper)
names(data4) <- c("year","y","median","low","upper")

require(ggplot2)
ggplot(data3)+
  geom_ribbon(aes(year,ymin=low,ymax=upper,colour="95% Prediction Interval"),alpha = 0.3, fill = 'green')+
  geom_line(aes(year,y,colour="Observations"))+
  geom_line(aes(year,median,colour="Prediction Median"))+
  geom_vline(xintercept = unlist(map4[1]), linetype='dotted')+
  geom_vline(xintercept = unlist(map4[2]), linetype='dotted')+
  geom_vline(xintercept = unlist(map4[3]), linetype='dotted')+
  geom_vline(xintercept = unlist(map4[4]), linetype='dotted')+
  annotate("text",label="MAP(T1) = 1892", x=1905,y=50)+
  annotate("text",label="MAP(T2) = 1930", x=1935,y=25)+
  annotate("text",label="MAP(T2) = 1948", x=1955,y=75)+
  annotate("text",label="MAP(T1) = 1892", x=1905,y=50)+
  scale_color_manual("",breaks=c("95% Prediction Interval","Observations","Prediction Median"),values = c("green","black","red"))+
  xlab("Years")+
  ylab("Cummulative number of disasters")+
  theme_bw()+
  theme(legend.position = "top")

ggsave("four_cp_predictions.png", height = 6, width = 8)

# Task 3 - Marginal Likelihood Evaluation on no. of changepoints ----------

## 1 Changepoint
pT <- rep(0,length(2:111))
idx <- 1
for (t in 2:111) {
  pT[idx] <-  log(1/choose(111,1))+(
              (log(112-t)*(-2-sum(y[t:112])) + lgamma(sum(y[t:112])+2)+
              log(t)*(-2-sum(y[1:(t-1)])) + lgamma(sum(y[1:(t-1)])+2))-
              (sum(log(factorial(y))))
              )
  idx <- idx+1
}
marginal1 <- sum(exp(pT))

## 2 Changepoints
T <- seq(2,111)
Tgrid <- expand.grid(T,T) 
names(Tgrid) <- c("t1","t2")
Tgrid <- Tgrid[Tgrid$t1 < Tgrid$t2,]

pT <- nrow(Tgrid)
idx <- 1
for (i in 1:nrow(Tgrid)) {
    pT[idx] <-  log(1/choose(111,2))+
                log(114-Tgrid[i,2])*(-2-sum(y[Tgrid[i,2]:112]))+ lgamma(sum(y[Tgrid[i,2]:112])+2)+
                log(Tgrid[i,2]-Tgrid[i,1])*(-2-sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)]))+ lgamma(sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)])+2)+
                log(Tgrid[i,1])*(-2-sum(y[1:(Tgrid[i,1]-1)])) + lgamma(sum(y[1:(Tgrid[i,1]-1)])+2)-
                (sum(log(factorial(y))))
    
    idx <- idx+1
}
marginal2 <- sum(exp(pT))

## 3 Changepoints
T <- seq(2,111)
Tgrid <- expand.grid(T,T,T) 
names(Tgrid) <- c("t1","t2","t3")
Tgrid <- Tgrid[Tgrid$t1 < Tgrid$t2 & Tgrid$t2 < Tgrid$t3,]

# Marginal T Distribution - p(T|y)
pT <- nrow(Tgrid)
idx <- 1
for (i in 1:nrow(Tgrid)) {
  pT[idx] <-  log(1/choose(111,3))+
              log(114-Tgrid[i,3])*(-2-sum(y[Tgrid[i,3]:112]))+ lgamma(sum(y[Tgrid[i,3]:112])+2)+
              log(Tgrid[i,3]-Tgrid[i,2])*(-2-sum(y[Tgrid[i,2]:(Tgrid[i,3]-1)]))+ lgamma(sum(y[Tgrid[i,2]:(Tgrid[i,3]-1)])+2)+
              log(Tgrid[i,2]-Tgrid[i,1])*(-2-sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)]))+ lgamma(sum(y[Tgrid[i,1]:(Tgrid[i,2]-1)])+2)+
              log(Tgrid[i,1])*(-2-sum(y[1:(Tgrid[i,1]-1)])) + lgamma(sum(y[1:(Tgrid[i,1]-1)])+2)-
              (sum(log(factorial(y))))
  
  idx <- idx+1
}
marginal3 <- sum(exp(pT))

## Visualise
barplot(c(marginal1,marginal2,marginal3), names.arg = c("1","2","3"), xlab = "Number of changepoints considered", ylab = "Marginal likelihood for a model", cex.names = 1.5, cex.lab = 1.5, cex.axis = 1.2)

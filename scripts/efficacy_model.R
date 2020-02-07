## Install packages

install.packages('MCMCglmm')
install.packages("PerformanceAnalytics")
install.packages("rgl")
## Load packages
library("MCMCglmm")
library("lattice")
library(corrplot)
library(parallel)
library(coda)
library(rgl)
library(ggplot2)
library(plyr)
## set working directory
wd <- getwd()
setwd("C:/R-3.4.2/bin/efficacy_analysis")
survival<-read.csv('~/Thesis references/Data/Data bases/CSV files/Survival_and_response_rates_Experiments_Clinical.csv',
		header = TRUE, sep = ';', quote = "\"", dec = ',')
head(survival,10)

survival$MOS_COR<-as.numeric(sub(',', '.',levels(survival$MOS_COR)))[survival$MOS_COR]
survival$RR<-as.numeric(sub(',', '.',levels(survival$RR)))[survival$RR]
survival<-within(survival, MAB<-relevel(MAB, 'Control'))
survival<-within(survival, STRAIN<-relevel(STRAIN,'C57BL/6'))
survival<-within(survival, CELL<-relevel(CELL,'B16F10'))

therapy_indx <- survival$MAB=='antiPDL1'|survival$MAB=='antiPD1'|survival$MAB=='antiCTLA4'|survival$MAB=='antiCTLA4+antiPD1'|survival$MAB=='antiCTLA4+antiPDL1'|survival$MAB=='Control'
therapy_indx2 <- survival$MAB=='antiPDL1'|survival$MAB=='antiPD1'|survival$MAB=='antiCTLA4'|survival$MAB=='antiCTLA4+antiPD1'|survival$MAB=='antiCTLA4+antiPDL1'|survival$MAB=='Control'|survival$MAB=='iIDO'|survival$MAB=='antiOX40'|survival$MAB=='antiGITR'|survival$MAB=='antiCTLA4+iIDO'|survival$MAB=='antiCTLA4+antiGITR'

survival_subset <-subset(survival, !is.nan(survival$MOS_COR)&therapy_indx, select = c(STUDY_ID,MOS_COR, CELL, STRAIN, MAB, LAB,N))
survival_subset2 <-subset(survival, !is.nan(survival$MOS_COR)&therapy_indx2, select = c(STUDY_ID,MOS_COR, CELL, STRAIN, MAB, LAB,N))

survival_subset$MAB<-droplevels(survival_subset$MAB)
survival_subset$CELL<-droplevels(survival_subset$CELL)
survival_subset$LAB<-droplevels(survival_subset$LAB)

survival_subset2$MAB<-droplevels(survival_subset2$MAB)
survival_subset2$CELL<-droplevels(survival_subset2$CELL)
survival_subset2$LAB<-droplevels(survival_subset2$LAB)


levels(survival_subset$MAB)
levels(survival_subset$LAB)
levels(survival_subset$CELL)
summary(survival_subset$CELL)
summary(survival_subset$MAB)
beta_survival<-aggregate(log(survival_subset$MOS_COR),list(survival_subset$MAB), mean)

exp(beta_survival[1,2])/(exp(beta_survival[2,2]))
exp(beta_survival[1,2])/(exp(beta_survival[3,2]))
exp(beta_survival[1,2])/(exp(beta_survival[4,2]))
exp(beta_survival[1,2])/(exp(beta_survival[5,2]))
exp(beta_survival[1,2])/(exp(beta_survival[6,2]))

nrow(survival_subset2)
# Using the subset of data
# Model 0: Fixed effect:MAB
prior0 = list(R = list(V = 1, nu = .002))

model0_glm <- MCMCglmm(log(MOS_COR)~MAB, family = "gaussian", data = survival_subset,
			prior=prior0,pr=TRUE,nitt=1e6, burnin=1e4, thin=100, singular.ok = TRUE)
summary(model0_glm)
plot(model0_glm)

% Model 1: Fixed effect:MAB+Random effect:CELL
prior1 = list(R = list(V = 1, nu = .002), G = list(G1 = list(V = 1,
     		nu =.02)))
model1_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL, family = "gaussian", data = survival_subset,
			prior=prior1,pr=TRUE,nitt=1e6, burnin=1e4, thin=100, singular.ok = TRUE)

%% Model 2: Fixed effect:MAB+Random effect:LAB
prior2 = list(R = list(V = 1, nu = .002), G = list(G1 = list(V = 1,
     		nu =.02)))

model2_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~LAB, family = "gaussian", data = survival_subset,
			prior=prior2,pr=TRUE,nitt=1e6, burnin=1e4, thin=100, singular.ok = TRUE)


%% Model 3
prior3 = list(R = list(V = 1, nu = .002), G = list(G1 = list(V = 1,
     		nu =.02), G2 = list(V = 1, nu = .02)))
model3_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL+LAB, data=survival_subset,
				prior = prior3,nitt=1e6,burnin=1e4, thin=1e2)


# Model 4: Fixed effect:MAB+RE:CELL+LAB+(MAB\CELL)
prior4 = list(R = list(V = diag(1), nu = 0.2),
     		G = list(G1 = list(V = diag(6), nu = 0.2),
         	G2 = list(V = 1, nu = 0.2)))
model4_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~idh(MAB):CELL+LAB, data=survival_subset,
				prior = prior4,nitt=1e6,burnin=1e4, thin=1e2)



# Model 5: Fixed effect:MAB+RE:STUDY+CELL
prior5 = list(R = list(V = 1, nu = .002), G = list(G1 = list(V = 1,
     		nu =.02), G2 = list(V = 1, nu = .02)))

model5_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL+STUDY_ID, data=survival_subset,
				prior = prior5,nitt=1e6,burnin=1e4, thin=1e2)

# Model 6: Fixed effect:MAB+RE:STUDY+CELL+LAB
prior6 = list(R = list(V = 1, nu = .002), G = list(G1 = list(V = 1,
     		nu =.02), G2 = list(V = 1, nu = .02),G3 = list(V = 1, nu = .02)))

model6_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL+STUDY_ID+LAB, data=survival_subset,
				prior = prior6,nitt=1e6,burnin=1e4, thin=1e2)


# Model 7: Fixed effect:MAB+RE:(MAB\CELL))
prior7 = list(R = list(V = diag(1), nu = 0.2),
     		G = list(G1 = list(V = diag(6), nu = 0.2)))
model7_glm <- MCMCglmm(log(MOS_COR)~STUDY_ID, random=~idh(MAB):CELL, data=survival_subset,
				prior = prior7,nitt=1e6,burnin=1e4, thin=1e2)



%% comparing all 4 models
summary(model0_glm)
summary(model1_glm)
summary(model2_glm)
summary(model3_glm)
DIC <-c(summary(model0_glm)$DIC,summary(model1_glm)$DIC,summary(model2_glm)$DIC,
		summary(model3_glm)$DIC,summary(model4_glm)$DIC,summary(model5_glm)$DIC,
		summary(model6_glm)$DIC,summary(model7_glm)$DIC)
data.frame(DIC,row.names=c('Model0', 'Model1', 'Model2', 'Model3', 'Model4', 'Model5', 'Model6', 'Model7'  ))

model1_mc <- mclapply(1:4, function(i) {
  MCMCglmm(log(MOS_COR)~MAB,random=~LAB+CELL, 
		data=survival_subset,
		family='gaussian', 
		prior=prior,thin=100, 
		burnin=1e3, nitt=1e5)
}, mc.cores=1)

model1_mc<-lapply(model1_mc,function(m) m$Sol)
model1_mc<-do.call(mcmc.list,model1_mc)
par(mfrow=c(3,2)mar=c(2,2,1,2))
gelman.plot(model1_mc,auto.layout=F)
gelman.diag(model1_mc)


par(mfrow=c(6,2), mar=c(2, 1, 1, 1))
plot(model1_mc, ask=F, auto.layout=F)



# Identifiability of fixed effects
fixed.effects.cormat <-cor(model3_glm$Sol[,])
random.effects.cormat <-cor(model3_glm$VCV)
corrplot(fixed.effects.cormat[1:6,1:6],method='color')
corrplot(random.effects.cormat,method='color')


% Autocorrelation of variance components samples
diag(autocorr(model3_glm$VCV)[2,,])
effectiveSize(model3_glm$VCV)

% Proportion of variance explained by each term
HPDinterval(model3_glm$VCV[,1:3]/rowSums(model3_glm$VCV),prob=.9)

% Posterior prediction intervals of the proportional increases in median survival

HR_HPD_3<-HPDinterval(exp(-(model3_glm$Sol[,2:6])))
HR_modes_3<-posterior.mode(exp(-(model3_glm$Sol[,2:6])))
HR_CI_3<-cbind(HR_modes_3,HR_HPD_3)
HR_HPD_0<-HPDinterval(exp(-(model0_glm$Sol[,2:6])))
HR_modes_0<-posterior.mode(exp(-(model0_glm$Sol[,2:6])))
HR_CI_0<-cbind(HR_modes_0,HR_HPD_0)
HR_HPD_1<-HPDinterval(exp(-(model1_glm$Sol[,2:6])))
HR_modes_1<-posterior.mode(exp(-(model1_glm$Sol[,2:6])))
HR_CI_1<-cbind(HR_modes_1,HR_HPD_1)

HR_HPD_2<-HPDinterval(exp(-(model2_glm$Sol[,2:6])))
HR_modes_2<-posterior.mode(exp(-(model2_glm$Sol[,2:6])))
HR_CI_2<-cbind(HR_modes_2,HR_HPD_2)



HR_CI_0
HR_CI_1
HR_CI_2
HR_CI_3

# Predicitons

survival_subset$predict<-predict(model3_glm,type ="response",marginal=~LAB+CELL, posterior='mode')
survival_subset$residuals<- survival_subset$predict-log(survival_subset$MOS_COR)
ggplot(survival_subset, aes(x=log(MOS_COR),y=residuals,colour=CELL))+
	geom_point()+
	facet_wrap(~LAB)

plot(log(survival_subset$MOS_COR),(survival_subset$predict), xlim=c(0, 5), ylim=c(0, 5))
abline(0,1)

beta_MAB = c(posterior.mode(model1_glm$Sol[,1]), posterior.mode(model1_glm$Sol[,2]),
	posterior.mode(model1_glm$Sol[,3]),posterior.mode(model1_glm$Sol[,4]),
	posterior.mode(model1_glm$Sol[,5]),posterior.mode(model1_glm$Sol[,6]))
plot(predict~MAB,data=survival_subset)

ggplot(survival_subset, aes(x=log(MOS_COR),y=predict,colour=CELL))+
	geom_point()+
	facet_wrap(~CELL)

% Posterior means and 95% C.I
plot.estimates <- function(x) {
  if (class(x) != "summary.mcmc")
    x <- summary(x)
  n <- dim(x$statistics)[1]
  par(mar=c(2, 7, 4, 1))
  plot(x$statistics[,1], n:1,
       yaxt="n", ylab="",
       xlim=range(x$quantiles)*1.1,
       pch=15,
       main="Posterior means and 95% credible intervals")
  grid()
  axis(2, at=n:1, rownames(x$statistics), las=2)
  arrows(x$quantiles[,1], n:1, x$quantiles[,5], n:1, code=0)
  abline(v=0, lty=2)
}
plot.estimates(model1_mc)

Vcell <- diag(colMeans(model2_glm$VCV))
summary(model2_glm)
plot(model2_glm)

# Identifiability of fixed effects
fixed.effects.cormat2 <-cor(model3_glm$Sol[,])
random.effects.cormat2 <-cor(model3_glm$VCV)

corrplot(fixed.effects.cormat2[1:6,1:6],method='color')
corrplot(random.effects.cormat2,method='color')


# Proportion of variance explained by each term
HPDinterval(model0_glm$VCV[,1]/rowSums(model0_glm$VCV))
HPDinterval(model1_glm$VCV[,1:2]/rowSums(model1_glm$VCV))
HPDinterval(model2_glm$VCV[,1:2]/rowSums(model2_glm$VCV))
HPDinterval(model3_glm$VCV[,1:3]/rowSums(model3_glm$VCV))
HPDinterval(model4_glm$VCV[,1:8]/rowSums(model4_glm$VCV))
HPDinterval(model5_glm$VCV[,1:3]/rowSums(model5_glm$VCV))

HPDinterval(model6_glm$VCV[,1:4]/rowSums(model6_glm$VCV))

# Counting factors
count(survival_subset, c('LAB', 'STUDY_ID'))
count(survival_subset, c('LAB', 'CELL'))
count(survival_subset, c('CELL', 'LAB'))
count(survival_subset, c('MAB', 'CELL'))

# Define extensions of .csv files containing the survival curves from preclinical expts.
library(plyr)
library(ggplot2)
library(survival)

# Define functions to use
processSurvivalData <- function(df,N, Time){
  ncol = dim(df)[2]
  nrow = dim(df)[1]
  indx = 1
  experiment_levels = factor(names(df)[seq(1,ncol,2)], ordered = F)
  experiment_list = strsplit(names(df)[seq(1,ncol,2)], "_")
  
  cell_line_levels = unlist(lapply(experiment_list,function(x)x[length(x)]))
  treatment_levels = (lapply(experiment_list,function(x)x[1:(length(x)-1)]))
  treatment_levels_1 = unlist(lapply(treatment_levels, function(x)x[1]))
  treatment_levels_2 = unlist(lapply(treatment_levels, function(x)x[2]))
  treatment_levels_3 = unlist(lapply(treatment_levels, function(x)x[3]))
  
  control_indx =  treatment_levels_1=='Control'
  experiment_levels_ordered = experiment_levels
  experiment_levels_ordered[control_indx] <-experiment_levels[control_indx]
  experiment_levels_ordered[!control_indx] <- experiment_levels[!control_indx]
  experiment_levels_ordered <- factor(experiment_levels_ordered,experiment_levels_ordered,experiment_levels_ordered, ordered = TRUE)
  
  survival <-data.frame(Time = rep(0,length(experiment_levels)),
                        Survival = rep(0,length(experiment_levels)),
                        AtRisk = rep(0, length(experiment_levels)),
                        Deaths = rep(0, length(experiment_levels)),
                        Censored = rep(0, length(experiment_levels)),
                        N_mice = rep(0, length(experiment_levels)),
                        Experiment = experiment_levels_ordered,
                        Cell = cell_line_levels,
                        Treatment_1 = treatment_levels_1,
                        Treatment_2 = treatment_levels_2,
                        Treatment_3 = treatment_levels_3)
  for (i in 1:(ncol/2)){
    time_i = as.numeric(levels(df[2:nrow,(i-1)*2 + 1]))
    nanIndx = !is.nan(time_i)
    time_i = sort(round(time_i[nanIndx], digits = 0))
    S_t_i = as.numeric(levels(df[2:nrow, (i-1)*2 + 2]))[df[2:nrow, (i-1)*2 + 2]]
    S_t_i = sort(round(S_t_i[nanIndx]*100)/100, decreasing = TRUE)
    N_i = N[i]
    
    survivors = round(S_t_i*N_i)
    deaths = c(0 ,diff(survivors)*(-1))
    atRisk = survivors+deaths
    censored_i = time_i >= Time[i]
    survival[indx:(indx+length(time_i)-1),1] = time_i
    survival[indx:(indx+length(time_i)-1),2] = S_t_i
    survival[indx:(indx+length(time_i)-1),3] = atRisk
    survival[indx:(indx+length(time_i)-1),4] = deaths
    survival[indx:(indx+length(time_i)-1),5] = censored_i
    survival[indx:(indx+length(time_i)-1),6] = N_i
    
    survival[indx:(indx+length(time_i)-1),7] = experiment_levels[i]
    survival[indx:(indx+length(time_i)-1),8] = cell_line_levels[i]
    survival[indx:(indx+length(time_i)-1),9] = treatment_levels_1[i]
    survival[indx:(indx+length(time_i)-1),10] = treatment_levels_2[i]
    survival[indx:(indx+length(time_i)-1),11] = treatment_levels_3[i]
    
    
    indx = indx+length(time_i)
  }
  survival$Experiment <- ordered(survival$Experiment, experiment_levels_ordered,experiment_levels_ordered)
  return(survival)
}
plotSurvival <- function(df){
  cell<-ggplot(df,aes(y=(Survival),x=Time,colour=Experiment)) +
    geom_point() +
    geom_step(data=df,aes(y=(Survival),x=Time))+ 
    theme(axis.text.x = element_text(angle = 45))
  
  cell+labs(color='Treatment ',  y='Survival',
            title='Survival in preclinical experiments',
            subtitle='Effect of treatment')
  
}
expandSurvivalFunction <- function(surv.df){
  Time <- rep(surv.df$Time, surv.df$Deaths+surv.df$Censored)
  Event <- rep(!surv.df$Censored, surv.df$Deaths+surv.df$Censored)
  Experiment <- rep(surv.df$Experiment, surv.df$Deaths+surv.df$Censored)
  Cell <-rep(surv.df$Cell, surv.df$Deaths+surv.df$Censored)
  Treatment_1 <- rep(surv.df$Treatment_1, surv.df$Deaths+surv.df$Censored)
  Treatment_2 <- rep(surv.df$Treatment_2, surv.df$Deaths+surv.df$Censored)
  Treatment_3 <- rep(surv.df$Treatment_3, surv.df$Deaths+surv.df$Censored)
  N_mice <- rep(surv.df$N_mice, surv.df$Deaths+surv.df$Censored)
  survival <- data.frame(Time, Event,N_mice, Experiment,Cell, Treatment_1, Treatment_2, Treatment_3)
  return(survival)
}
SurvMedians <- function(df){
  med_hat <- log(2)*sum(df$Time)/sum(df$Event)
  var_Bartholomew <- med_hat^2/sum(1-exp(-df$Time*log(2)/med_hat))
  CI_LB_Barth <- med_hat-1.96*sqrt(var_Bartholomew)
  CI_UB_Barth <-med_hat + 1.96*sqrt(var_Bartholomew)
  d_Tot <- sum(df$Event)
  lambda_hat <- sum(df$Event)/sum(df$Time)
  CI_LB_VarEst <- log(2)/exp(log(lambda_hat)+d_Tot^(-1/2)*1.96)
  CI_UB_VarEst <- log(2)/exp(log(lambda_hat)-d_Tot^(-1/2)*1.96)
  medians.df <-data.frame(Median = med_hat,
                          CI_LB_Barth, CI_UB_Barth, CI_LB_VarEst, CI_UB_VarEst)
  return(medians.df)
}
createDF <- function(medians.list,df, publication_name){
  indx <- match( names(medians.list), df$Experiment)
  data.frame(Study=rep(publication_name, length(names(medians.list))),
             Treatment = names(medians.list),
             Median=unlist(lapply(medians.list, function(x) x[c(1)])),
             LB=unlist(lapply(medians.list, function(x) x[c(4)])),
             UB=unlist(lapply(medians.list, function(x) x[c(5)])), row.names = c(),
             Cell=df$Cell[indx],
             N_mice = df$N_mice[indx],
             Treatment_1 = df$Treatment_1[indx],
             Treatment_2 = df$Treatment_2[indx],
             Treatment_3 = df$Treatment_3[indx]
  )
}
phi.u <- function(m,n,X,Y,psi){
  phi_tilde<-m*n/(m+n)^2*(psi*(n-X)*Y/n*(1+(psi-1)*Y/m)+X*(m-Y)/m*(psi-(psi-1)*X/n))
  return(phi_tilde)
} 
MH.Estimate <- function(df){
  treatments <- df$Treatment_1
  control_indx <-which(treatments=="Control")
  treatments_indx <- which(treatments!="Control")
  treatments <- unique(df$Experiment[treatments_indx])
  betaMH <- data.frame(beta = rep(0, length(treatments)),
                       sigma = rep(0, length(treatments)),
                       CI_LB = rep(0, length(treatments)),
                       CI_UB = rep(0, length(treatments)), row.names = treatments)
  
  for (i in 1:length(treatments)){
    treatment_indx = which(df$Experiment==treatments[i])
    times = sort(unique(c(df$Time[control_indx], df$Time[treatment_indx])))
    times_indx_Control = match( df$Time[control_indx],times)
    times_indx_Treatment = match(df$Time[treatment_indx],times)
    beta.df <- data.frame(Time = times,
                          AtRisk_Control = rep(0, length(times)),
                          Deaths_Control = rep(0, length(times)),
                          AtRisk_Treatment = rep(0, length(times)),
                          Deaths_Treatment = rep(0, length(times)))
    beta.df$Deaths_Control[times_indx_Control] = df$Deaths[control_indx]
    beta.df$Deaths_Treatment[times_indx_Treatment] = df$Deaths[treatment_indx]
    beta.df$AtRisk_Control = df$AtRisk[control_indx[1]]-cumsum(beta.df$Deaths_Control)+beta.df$Deaths_Control
    beta.df$AtRisk_Treatment = df$AtRisk[treatment_indx[1]]-cumsum(beta.df$Deaths_Treatment)+beta.df$Deaths_Treatment
    
    beta.df$AtRisk_Control[2:length(times)]<-beta.df$AtRisk_Control[1:length(times)-1]- beta.df$Deaths_Control[1:length(times)-1]
    beta.df$AtRisk_Treatment[2:length(times)]<-beta.df$AtRisk_Treatment[1:length(times)-1] - beta.df$Deaths_Treatment[1:length(times)-1]
    
    
    beta.df$S_k <- (beta.df$Deaths_Control*(beta.df$AtRisk_Treatment-beta.df$Deaths_Treatment)/(beta.df$AtRisk_Control+beta.df$AtRisk_Treatment))
    beta.df$R_k <-(beta.df$Deaths_Treatment*(beta.df$AtRisk_Control - beta.df$Deaths_Control)/(beta.df$AtRisk_Control+beta.df$AtRisk_Treatment))
    betaOR <- sum(beta.df$R_k)/sum(beta.df$S_k)
    beta.df$phi_k <- phi.u(beta.df$AtRisk_Control,
                           beta.df$AtRisk_Treatment,
                           beta.df$Deaths_Treatment, 
                           beta.df$Deaths_Control, 
                           betaOR)
    betaVar <-sum(beta.df$phi_k/length(times),na.rm = T)/sum(beta.df$S_k/length(times),na.rm=T)^2/length(times)
    ln_betaVar <- betaVar/betaOR^2
    betaLB <- exp(log(betaOR)-sqrt(ln_betaVar)*1.96)
    betaUB <- exp(log(betaOR) + sqrt(ln_betaVar)*1.96)
    
    betaMH$beta[i]<-betaOR
    betaMH$simga[i]<- sqrt(betaVar)
    betaMH$CI_LB[i] <- betaLB
    betaMH$CI_UB[i] <- betaUB
  }
  return(betaMH)
}
CoxModel <-function(df){
  treatments <- df$Treatment_1
  control_indx <-which(treatments=="Control")
  treatments_indx <- which(treatments!="Control")
  treatments <- unique(df$Experiment[treatments_indx])
  HR_Cox <- data.frame(Treatment = treatments,
                       HR = rep(0, length(treatments)),
                       CI_LB = rep(0, length(treatments)),
                       CI_UB = rep(0, length(treatments)),
                       logrank = rep(0, length(treatments)),
                       p_value = rep(0, length(treatments)))
  control_subset <-  df[df$Treatment_1==c("Control"),]
  for(i in 1:length(treatments)){
    subset<- rbind(control_subset, df[df$Experiment==as.character(treatments[i]),])
    subset$Experiment <- droplevels(subset$Experiment)
    subset$Experiment <- factor(subset$Experiment, ordered = F)
    cox.model <- coxph(Surv(Time, Event) ~ Experiment, data = subset)
    summ.cox.model <-summary(cox.model)
    HR <- summ.cox.model$conf.int[1]
    HR_LB <- summ.cox.model$conf.int[3]
    HR_UB <- summ.cox.model$conf.int[4]
    logrank <- summ.cox.model$sctest[1]
    pval <- summ.cox.model$sctest[3]
    HR_Cox[i,2:6] = c(HR, HR_LB, HR_UB, logrank, pval)
  }
  return(HR_Cox)
}

# Loading data
setwd("~/Documents/Data/CSV/Survival/Datasets")
study <- sub('./','', list.dirs()[-1])
temp <- (lapply(list.dirs()[-1], function(x) list.files(path = x, pattern = '*.csv')))

figures <- sub('(.*)_Fig', '', unlist(temp))
figures <- sub('(.*)_SF', 'SF',figures) 
figures <- sub('.csv', '', figures)
studies <- sub('_Fig(.*)', '', unlist(temp))
studies <- sub('_SF(.*)', '', studies)

fileExt <- paste(getwd(), studies, unlist(temp), sep = '/')
datasets <- lapply(fileExt, read.csv)

experiments <- lapply(datasets, function(x) names(x)[seq(1, length(names(x)), 2)])

n_studies <- unlist(lapply(experiments, length))
studies.factor <- rep(studies, n_studies)
figures.factor <- rep(figures, n_studies)
experiments.factor <- data.frame((studies.factor), figures.factor, (unlist(experiments)))
names(experiments.factor) <- c('STUDY_ID', 'FIG', 'EXPERIMENT')
head(experiments.factor)

Survival_Data <- read.csv2("~/Documents/Data/CSV/Survival/Survival_Data.csv")
head(Survival_Data)
EXPERIMENT <- paste(Survival_Data$TREATMENT_1, 
                                  Survival_Data$TREATMENT_2, 
                                  Survival_Data$TREATMENT_3,
                                  Survival_Data$CELL,sep='_')
EXPERIMENT <- sub('__', '_', EXPERIMENT) # apply twice
EXPERIMENT <- sub('__', '_', EXPERIMENT) # apply twice

Survival_Data$EXPERIMENT <- EXPERIMENT
head(Survival_Data)
# Match datasets to meta-data
survival.subset <-match_df(Survival_Data, experiments.factor, on = c('STUDY_ID', 'FIG', 'EXPERIMENT'))
differences <- data.frame(setdiff(experiments.factor$EXPERIMENT, survival.subset$EXPERIMENT))

# Check that all records are available
setdiff(experiments.factor$STUDY_ID, survival.subset$STUDY_ID)
setdiff(survival.subset$EXPERIMENT, experiments.factor$EXPERIMENT)

N_mice <-list()
finalT <-list()
for(i in 1:length(datasets)){
  indx <-  grepl(studies[i],survival.subset$STUDY_ID) & (survival.subset$EXPERIMENT %in%  experiments[[i]]) & survival.subset$FIG %in% figures[i]
  N_mice[[i]] <- survival.subset$N_mice[indx]
  finalT[[i]] <- survival.subset$T_F[indx]
  
}

dataSurvival <- list()
for(i in 1:length(datasets)){
  dataSurvival[[i]] = processSurvivalData(datasets[[i]], N_mice[[i]], finalT[[i]])
}
cell<-list()
cell<-lapply(dataSurvival, plotSurvival)

processSurvival.df <-lapply(dataSurvival,expandSurvivalFunction)
survival.df <- processSurvival.df[[1]]
survival.df$Study <- rep(studies[1], dim(survival.df)[1])
for(i in 2:length(studies)){
  subset_i <- processSurvival.df[[i]]
  subset_i$Study <- rep(studies[i], dim(subset_i)[1])
  survival.df <- rbind(survival.df, subset_i)
}

survivalMedian.df <-lapply(processSurvival.df,function(x)by(x, x$Experiment, SurvMedians))

for(i in 1:length(studies)){
  survivalMedian.df[[i]] <- createDF(survivalMedian.df[[i]], processSurvival.df[[i]], studies[i])
}

survMedian.df <- survivalMedian.df[[1]]
for(i in 2:length(studies)){
  survMedian.df <- rbind(survMedian.df, survivalMedian.df[[i]])
}

OS_PFS<-ggplot(survMedian.df,aes(y=Median,x=Cell,colour=paste(survMedian.df$Treatment_1,survMedian.df$Treatment_2,survMedian.df$Treatment_3,sep = '+')))+
  geom_point()+
  geom_errorbar(aes(ymin=LB, ymax=UB))

OS_PFS+labs(color='Treatment',  y='Time [days]',x="Treatment",
            title='Median survival in preclinical experiments',
            subtitle='Progression-free vs overall survival')
  
MS_PFS<-ggplot(survMedian.df,aes(y=Median,x=Treatment_1,colour=Cell,sep = '+'))+
  geom_point()+
  geom_errorbar(aes(ymin=LB, ymax=UB))

MS_PFS+labs(color='Treatment',  y='Time [days]',x="Treatment",
            title='Median survival in preclinical experiments',
            subtitle='Survival by cell line and first therapy')



HR_Estimates <- by(survival.df, survival.df$Study, CoxModel)
HR <- HR_Estimates$ALLARD_2013
HR$Study <- rep(names(HR_Estimates)[1], dim(HR)[1])
for(i in 2:length(HR_Estimates)){
  subset <- HR_Estimates[[i]]
  subset$Study <- rep(names(HR_Estimates)[i], dim(subset)[1])
  HR <- rbind(HR, subset)
}
HR_PFS<-ggplot(HR,aes(y=HR,x=Treatment,color = Treatment))+
  geom_point()
  

HR_PFS+labs(color='Treatment',  y='Time [days]',x="Treatment",
            title='Hazard ratio [Preclinical]',
            subtitle='Treatment vs control')


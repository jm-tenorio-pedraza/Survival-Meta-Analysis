processSurvivalData <- function(df,N, Time,exps, figure,study){
  ncol = dim(df)[2]
  nrow = dim(df)[1]
  indx = 1
  experiment_levels = factor(exps, ordered = F)
  experiment_list = strsplit(exps, "_")
  
  cell_line_levels = unlist(lapply(experiment_list,function(x)x[length(x)]))
  treatment_levels = (lapply(experiment_list,function(x)x[1:(length(x)-1)]))
  treatment_levels_1 = unlist(lapply(treatment_levels, function(x)x[1]))
  treatment_levels_2 = unlist(lapply(treatment_levels, function(x)x[2]))
  treatment_levels_3 = unlist(lapply(treatment_levels, function(x)x[3]))
  treatment_levels_4 = unlist(lapply(treatment_levels, function(x)x[4]))
  
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
                        Treatment_3 = treatment_levels_3,
                        Treatment_4 = treatment_levels_4)
  for (i in 1:(ncol/2)){
    timecolIndx <- (i-1)*2 +1;
    survcolIndx <- timecolIndx+1;
    time_i = (df[1:nrow,timecolIndx])
    nanIndx = !(is.na(time_i) | is.nan(time_i))
    time_i = sort(round(time_i[nanIndx], digits = 0))
    if(any(time_i<0)) time_i[time_i<0]<-0 # Correct for values that might be slightly below 0
     
    S_t_i <- df[1:nrow, (i-1)*2 + 2]
    S_t_i <- sort(round(S_t_i[nanIndx]*100)/100, decreasing = TRUE)
    if(all(time_i!=0) & S_t_i[1]<1){time_i<-c(0,time_i)
    S_t_i<-c(1,S_t_i)}
    if(any(S_t_i<0))S_t_i[S_t_i<0]<-0
    if(any(S_t_i>1))S_t_i[S_t_i>1]<-1
    
    N_i <- N[i]
    censTime_i <-time_i
    # Check that the censoring time is higher than or equal to the last observed time point 
    if(censTime_i[length(time_i)]>Time[i]) Time[i]<-censTime_i[length(time_i)]
    
    # add extra row if there is censoring and the row is not present
    if(ceiling(S_t_i[length(S_t_i)]*100)>=floor(1/N_i*100)){
      if (ceiling(S_t_i[length(S_t_i)-1]*100)-ceiling(S_t_i[length(S_t_i)]*100)>=floor(1/N_i*100)){
        time_i<-c(time_i, Time[i])    # Set the last observed survival to the final experiment time point if there is at least 1 survivor

        censTime_i<-c(censTime_i,Time[i])
        S_t_i<-c(S_t_i, S_t_i[length(S_t_i)])
        
      }
        else{censTime_i[length(time_i)] <-Time[i]}
      
    }  
    if(time_i[length(time_i)]-time_i[length(time_i)-1]<1){
      time_i<-time_i[1:(length(time_i)-1)]
      censTime_i<-censTime_i[1:(length(censTime_i)-1)]
      S_t_i<-S_t_i[1:(length(S_t_i)-1)]
    }

    survivors = round(S_t_i*N_i)
    deaths = c(0 ,diff(survivors)*(-1))
    atRisk = survivors+deaths
    censored_i = (censTime_i >= Time[i]) & survivors>0 
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
    survival[indx:(indx+length(time_i)-1),12] = treatment_levels_4[i]
    
    
    indx = indx+length(time_i)
  }
  survival$Experiment <- ordered(survival$Experiment, experiment_levels_ordered,experiment_levels_ordered)
  survival$Figure <- figure
  survival$Study <- study
  return(survival)
}
plotSurvival <- function(df, study, figure){
  cell<-ggplot(df,aes(y=(Survival),x=Time,colour=Experiment)) +
    geom_point() +
    geom_step(data=df,aes(y=(Survival),x=Time))+ 
    theme(axis.text.x = element_text(angle = 45))
  titleName <- paste("Survival in", study, sep = " ")
  subtitleName <- paste("Figure", figure, sep = " ")
  cell+labs(color='Treatment ',  y='Survival',
            title=titleName,
            subtitle= subtitleName)
  
}
expandSurvivalFunction <- function(surv.df){
  Time <- rep(surv.df$Time, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Event <- rep(!surv.df$Censored, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Experiment <- rep(surv.df$Experiment, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Cell <-rep(surv.df$Cell, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Treatment_1 <- rep(surv.df$Treatment_1, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Treatment_2 <- rep(surv.df$Treatment_2, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Treatment_3 <- rep(surv.df$Treatment_3, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Treatment_4 <- rep(surv.df$Treatment_4, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  N_mice <- rep(surv.df$N_mice, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Figure <- rep(surv.df$Figure, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  Experiment<- as.factor(paste(Experiment, Figure, sep='_'))
  Study <- rep(surv.df$Study, surv.df$Deaths+surv.df$Censored*(surv.df$AtRisk-surv.df$Deaths))
  
  
  survival <- data.frame(Time, Event,N_mice, Experiment,Cell, Treatment_1, Treatment_2, Treatment_3, Treatment_4,Figure,Study)
  return(survival)
}
SurvMedians <- function(df){
  NR<-df$N_mice[1]*.5>sum(df$Event)
  surv.df<-Surv(df$Time, df$Event)
  median.surv.est <-median(surv.df)
  med_hat <- log(2)*sum(df$Time)/sum(df$Event)
  med_samp <- median(df$Time)
  mean_samp <-mean(df$Time)
  var_Bartholomew <- med_hat^2/sum(1-exp(-df$Time*log(2)/med_hat))
  CI_LB_Barth <- med_hat-1.96*sqrt(var_Bartholomew)
  CI_UB_Barth <-med_hat + 1.96*sqrt(var_Bartholomew)
  d_Tot <- sum(df$Event)
  lambda_hat <- sum(df$Event)/sum(df$Time)
  CI_LB_VarEst <- log(2)/exp(log(lambda_hat)+d_Tot^(-1/2)*1.96)
  CI_UB_VarEst <- log(2)/exp(log(lambda_hat)-d_Tot^(-1/2)*1.96)
  medians.df <-data.frame(Median_Exp = med_hat,
                          Median_Surv = median.surv.est$quantile[1],
                          Median_Surv_LB = median.surv.est$lower[1],
                          Median_Surv_UB = median.surv.est$upper[1],
                          CI_LB_Barth, CI_UB_Barth, CI_LB_VarEst, CI_UB_VarEst,var_Bartholomew,
                          NR, 
                          MedS = med_samp, 
                          MeanS = mean_samp)
  return(medians.df)
}
createDF <- function(medians.list, df, publication_name){
  indx <- match( names(medians.list), df$Experiment)
  data.frame(Study=rep(publication_name, length(names(medians.list))),
             Experiment = names(medians.list),
             Median_Surv = unlist(lapply(medians.list, function(x) x['Median_Surv'])),
             Median_Surv_LB = unlist(lapply(medians.list, function(x) x['Median_Surv_LB'])),
             Median_Surv_UB = unlist(lapply(medians.list, function(x) x['Median_Surv_UB'])),
             Median_Exponential=unlist(lapply(medians.list, function(x) x['Median_Exp'])),
             LB=unlist(lapply(medians.list, function(x) x['CI_LB_VarEst'])),
             UB=unlist(lapply(medians.list, function(x) x['CI_UB_VarEst'])),
             Var=unlist(lapply(medians.list, function(x) x['var_Bartholomew'])),
             NR=unlist(lapply(medians.list, function(x) x['NR'])), row.names = c(),
             SampleMedian=unlist(lapply(medians.list, function(x) x['MedS'])),
             SampleMean=unlist(lapply(medians.list, function(x) x['MeanS'])),
             Cell=df$Cell[indx],
             N_mice = df$N_mice[indx],
             Treatment_1 = df$Treatment_1[indx],
             Treatment_2 = df$Treatment_2[indx],
             Treatment_3 = df$Treatment_3[indx],
             Treatment_4 = df$Treatment_4[indx]
             
  )
}
phi.u <- function(m,n,X,Y,psi){
  phi_tilde<-m*n/(m+n)^2*(psi*(n-X)*Y/n*(1+(psi-1)*Y/m)+X*(m-Y)/m*(psi-(psi-1)*X/n))
  return(phi_tilde)
} 
controlAssignment <- function(df){
  df$Experiment <- droplevels(df$Experiment)
  medianSurv <- by(df$Time, df$Experiment, median, simplify = T)
  medianSurv <- medianSurv[1:dim(medianSurv)]
  minMedian <- min(medianSurv)
  control <- names(medianSurv)[minMedian==medianSurv]
  control <- control[1]
  control_indx <- df$Experiment == control
  return(control_indx)
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
    betaLB <- exp(log(betaOR) - sqrt(ln_betaVar)*1.96)
    betaUB <- exp(log(betaOR) + sqrt(ln_betaVar)*1.96)
    
    betaMH$beta[i]<-betaOR
    betaMH$simga[i]<- sqrt(betaVar)
    betaMH$CI_LB[i] <- betaLB
    betaMH$CI_UB[i] <- betaUB
  }
  return(betaMH)
}
CoxModel <-function(df){
  # Identify number of figures
  figs<-df$Figure
  uniqueFigs<-unique(figs)
  # Identify control mice
  treatment_1 <- df$Treatment_1
  control_indx <-treatment_1=="Control"
  monotherapy_indx<-is.na(df$Treatment_2) & !control_indx
  is.element(df$Treatment_1,c('antiCTLA4','antiPD1','antiPDL1'))
  if(!any(control_indx)) control_indx <- controlAssignment(df)
  controlgroup<-unique(treatment_1[control_indx])
  # Create HR output structure
  HR_Cox <- data.frame(Experiment = NA,
                         Study = NA,
                         Figure = NA,
                         Cell = NA,
                         Treat_1 = NA,
                         Treat_2 = NA,
                         Treat_3 = NA,
                         Treat_4 = NA,
                         RefGroup=NA,
                         Median_Survival_logRatio = NA,
                         N=NA,
                         Treatment_Median_Sample=NA,
                         Control_Median_Sample=NA,
                         NR=NA,
                         Concordance=NA,
                         HR = NA,
                         SE_coef=NA,
                         HR_LB = NA,
                         HR_UB = NA,
                         logrank = NA,
                         p_value = NA)
  HR_indx<-1
  # For each figure:
  for(j in 1:length(uniqueFigs)){
    control_j_indx<-treatment_1=='Control' & figs==uniqueFigs[j]
    if(all(!control_j_indx)) control_j_indx <-control_indx
    # Identify treated mice in Figure j
    treatments_indx <- !control_indx & figs==uniqueFigs[j]
    # Create temporary df with control and mice from each treatment
    control_subset <-  df[control_j_indx,]
    control_subset$Experiment <-factor(paste('Control', control_subset$Cell, sep='_'))
    treatments_subset <- df[treatments_indx,]
    # Experiment names
    exps<- (treatments_subset$Experiment)
    # Unique treatments
    uniqueExp <- unique(exps)
    uniquetreatsIndx <- match(uniqueExp, exps)
    uniqueTreats <- treatments_subset[uniquetreatsIndx, c('Cell', 'Treatment_1', 'Treatment_2', 'Treatment_3','Treatment_4', 'Figure')]
    # Calculate median survival in control group
    control_median.df<-SurvMedians(control_subset)
    if(!any(treatments_indx)){
      HR_Cox[HR_indx,c('Experiment','Study','Figure','Cell','Treat_1','RefGroup')] <- c(as.factor(paste(control_subset$Experiment[1], control_subset$Figure[1],sep='_')),
                           control_subset$Study[1],
                           control_subset$Figure[1],
                           control_subset$Cell[1],
                           control_subset$Treatment_1[1],
                           controlgroup)
      HR_indx<-HR_indx+1
      
    }
    else{
      # df strcture for output 
      for(i in 1:length(uniqueExp)){
        HR_Cox[HR_indx,c('Experiment','Study','Figure',
                         'Cell','Treat_1','Treat_2',
                         'Treat_3','Treat_4','RefGroup')] <- c(as.character(uniqueExp[i]),(df$Study[1]),
                                                               uniqueTreats$Figure[i],uniqueTreats$Cell[i],
                                                               uniqueTreats$Treatment_1[i],
                                                               uniqueTreats$Treatment_2[i],uniqueTreats$Treatment_3[i],
                                                               uniqueTreats$Treatment_4[i],controlgroup)
        treatIndx_i <- exps==uniqueExp[i]
        subset <- rbind(control_subset, treatments_subset[treatIndx_i,])
        subset$Experiment <- droplevels((subset$Experiment))
        subset$Experiment <- factor(subset$Experiment, ordered = F)
        cox.model <- coxph(Surv(Time, Event) ~ Experiment, data = subset,ties='efron')
        treat_median.df<-SurvMedians(treatments_subset[treatIndx_i,])
        summ.cox.model <-summary(cox.model)
        HR <- summ.cox.model$conf.int[1]
        HR_LB <- summ.cox.model$conf.int[3]
        HR_UB <- summ.cox.model$conf.int[4]
        logrank <- summ.cox.model$sctest[1]
        pval <- summ.cox.model$sctest[3]
        conc<-summ.cox.model$concordance[1]
        se<-summ.cox.model$coefficients[1,'se(coef)']
        
        median_surv_logratio <-log(treat_median.df$Median_Surv/control_median.df$Median_Surv)
        med_sample <-treat_median.df$MedS
        NR<-treat_median.df$NR
        N<- dim(subset)[1]
        HR_Cox[HR_indx,c("Median_Survival_logRatio",'N','Treatment_Median_Sample','Control_Median_Sample','NR',"Concordance","HR",'SE_coef', "HR_LB", "HR_UB", "logrank", "p_value")]<-
          c(median_surv_logratio,N,med_sample,control_median.df$Median_Surv,NR,conc,HR,se, HR_LB, HR_UB, logrank, pval)
        HR_indx<-HR_indx+1
      }
    }
  }
  # HR_Cox$Experiment<-as.factor(HR_Cox$Experiment)
  return(HR_Cox)
  
}
CoxModel.modified <-function(df){
  # Identify number of figures
  figs<-df$Figure
  uniqueFigs<-unique(figs)
  # Identify control mice
  treatment_1 <- df$Treatment_1
  control_indx <-treatment_1=="Control"
  monotherapy<-is.na(df$Treatment_2) & !control_indx & is.element(df$Treatment_1,c('antiCTLA4','antiPD1','antiPDL1'))
  # Choose a different monotherapy as the reference when none of the anti-CTLA-4, anti-PD-1 or anti-PD-L1 therapies are present
  if(all(!monotherapy)) monotherapy<-is.na(df$Treatment_2) & !control_indx
  # Choose the control as the monotherapy when there is no other comparator group
  if(all(!monotherapy)) monotherapy<- control_indx
  # Assign the monotherapy referece based on the least effective
  monotherapy_indx<-controlAssignment(df[monotherapy,])
  if(!any(control_indx)) control_indx <- controlAssignment(df)
  controlgroup<-unique(treatment_1[control_indx])
  monotherapygroup<-unique(treatment_1[monotherapy][monotherapy_indx])
  # Create HR output structure
  HR_Cox <- data.frame(Experiment = NA,
                       Study = NA,
                       Figure = NA,
                       Cell = NA,
                       Treat_1 = NA,
                       Treat_2 = NA,
                       Treat_3 = NA,
                       Treat_4 = NA,
                       RefGroup=NA,
                       Median_Survival_logRatio = 0,
                       N=0,
                       Treatment_Median_Sample=0,
                       Control_Median_Sample=0,
                       NR=FALSE,
                       Concordance=0,
                       HR = 0,
                       SE_coef=0,
                       HR_LB = 0,
                       HR_UB = 0,
                       logrank = 0,
                       p_value = 0)
  HR_indx<-1
  # For each figure:
  for(j in 1:length(uniqueFigs)){
    control_j_indx<-treatment_1=='Control' & figs==uniqueFigs[j]
    # Assign study-level control if no figure-specific control is available
    if(all(!control_j_indx)) control_j_indx <-control_indx
    # Identify treated mice in Figure j
    treatments_indx <- !control_indx & figs==uniqueFigs[j]
    treatments_subset <- df[treatments_indx,]
    # Experiment names
    exps<- (treatments_subset$Experiment)
    # Unique treatments
    uniqueExp <- unique(exps)
    uniquetreatsIndx <- match(uniqueExp, exps)
    uniqueTreats <- treatments_subset[uniquetreatsIndx, c('Cell', 'Treatment_1', 'Treatment_2', 'Treatment_3','Treatment_4', 'Figure')]
    
    if(!any(treatments_indx)){
      control_subset <-  df[control_j_indx,]
      control_subset$Experiment <-factor(paste('Control', control_subset$Cell, sep='_'))
      HR_Cox[HR_indx,c('Experiment','Study','Figure','Cell','Treat_1','RefGroup')] <- c(as.factor(paste(control_subset$Experiment[1], control_subset$Figure[1],sep='_')),
                                                                                        control_subset$Study[1],
                                                                                        control_subset$Figure[1],
                                                                                        control_subset$Cell[1],
                                                                                        control_subset$Treatment_1[1],
                                                                                        controlgroup)
      HR_indx<-HR_indx+1
      
    }
    else{
      # df strcture for output 
      for(i in 1:length(uniqueExp)){
        HR_Cox[HR_indx,c('Experiment','Study','Figure',
                         'Cell','Treat_1','Treat_2',
                         'Treat_3','Treat_4')] <- c(as.character(uniqueExp[i]),(df$Study[1]),
                                                               uniqueTreats$Figure[i],uniqueTreats$Cell[i],
                                                               uniqueTreats$Treatment_1[i],
                                                               uniqueTreats$Treatment_2[i],uniqueTreats$Treatment_3[i],
                                                               uniqueTreats$Treatment_4[i])
        treatIndx_i <- exps==uniqueExp[i]
        # Assign monotherapy control in case of combination therapy
        if(!is.na(treatments_subset[treatIndx_i,'Treatment_2'][1])){
          control_subset<-df[monotherapy,][monotherapy_indx,];refgroup_i<-monotherapygroup
          } else {
            control_subset <-  df[control_j_indx,];control_subset$Experiment <-factor(paste('Control', control_subset$Cell, sep='_'));refgroup_i<-controlgroup}
        # Calculate median survival in control group
        control_median.df<-SurvMedians(control_subset)
        subset <- rbind(control_subset, treatments_subset[treatIndx_i,])
        subset$Experiment <- droplevels((subset$Experiment))
        subset$Experiment <- factor(subset$Experiment,c(as.character(control_subset$Experiment[1]),as.character(uniqueExp[i])), ordered = T)
        cox.model <- coxph(Surv(Time, Event) ~ Experiment, data = subset,ties='efron')
        treat_median.df<-SurvMedians(treatments_subset[treatIndx_i,])
        summ.cox.model <-summary(cox.model)
        HR <- summ.cox.model$conf.int[1]
        HR_LB <- summ.cox.model$conf.int[3]
        HR_UB <- summ.cox.model$conf.int[4]
        logrank <- summ.cox.model$sctest[1]
        pval <- summ.cox.model$sctest[3]
        conc<-summ.cox.model$concordance[1]
        se<-summ.cox.model$coefficients[1,'se(coef)']
        
        median_surv_logratio <-log(treat_median.df$Median_Surv/control_median.df$Median_Surv)
        med_sample <-treat_median.df$MedS
        NR<-treat_median.df$NR
        N<- dim(subset)[1]
        HR_Cox$RefGroup[HR_indx]<-refgroup_i
        HR_Cox[HR_indx,c("Median_Survival_logRatio",'N','Treatment_Median_Sample','Control_Median_Sample','NR',"Concordance","HR",'SE_coef', "HR_LB", "HR_UB", "logrank", "p_value")]<-
          c(median_surv_logratio,N,med_sample,control_median.df$Median_Surv,NR,conc,HR,se, HR_LB, HR_UB, logrank, pval)
        HR_indx<-HR_indx+1
      }
    }
  }
  # HR_Cox$Experiment<-as.factor(HR_Cox$Experiment)
  return(HR_Cox)
  
}

doseExpansion<-function(df){
  id<-1:dim(df)[1]
  df$id<-(id)
  dosingTimes.ls<-list()
  df$dosingTimes<-df$BASELINE_days
 
  for(j in 1:4){
    freq_j<-paste('FREQ', j, sep='_')
    dose_j<-paste('DOSE',j,sep='_')
    cummdose_j<-paste('CUMMULATIVE_DOSE', j,sep='_')
    dosingTimes_j<-paste('dosingTimes',j,sep='')
    treatStart_j<-paste('TREAT_START',j,sep='_')
    ncycles_j<-paste('ncycles', 1, sep='')
    nlag_j<-paste('schedule1_nlag', j,sep='')
    ndoses_j<-paste('schedule1_ndoses', j,sep = '')
    sch2_nlag_j<-paste('schedule2_nlag', j,sep='')
    sch2_ndoses_j<-paste('schedule2_ndoses', j,sep='')
    
    # Rows of subjects that received no dose of treatment_j
    na_indx <- is.na(df[,freq_j]) | is.element(df[,freq_j],'Control')
    
    
    # Subset of individuals that received more that 1 dose of treatment_j
    treat.df<-subset(df,!(na_indx))
    if (dim(treat.df)[1]!=0){
    treat.df<-droplevels(treat.df)
    
    sd_indx <- is.element(treat.df[,freq_j],'SD')
    
    # treat1.df<-treat.df[!sd_indx,]
    # treat1.df<-droplevels(treat1.df)
    # Identify individuals that had two different treatment schedules
    cycleINDX<-grepl('+',treat.df[!sd_indx,freq_j],fixed=T)

    # Obtain the first treatment schedules
    schedule1<-gsub('\\+.*','', treat.df[!sd_indx,freq_j])
    schedule2<-gsub('.*\\+','',treat.df[!sd_indx,freq_j])
    schedule2[!cycleINDX]<-NA
    
    # Obtain the number of cycles of the first treatment schedule
    cycles1 <-sub('x*Q.{1,2}Dx.+','',schedule1)
    ncycles1 <- as.numeric((sub('Q.*W','',cycles1)))
    ncycles1[is.na(ncycles1)]<-1
    
    # Obtain the lag between cycles of the first treatment schedule
    cyclelag1 <- as.numeric(sub('.*Q','',sub(c('W.*'),'', cycles1)))
    cyclelag1[is.na(cyclelag1)]<-0
    
    # Obtain the lag between doses of the first treatment schedule
    doseCycle1<-sub('.+Wx','',schedule1)
    nlag1<- as.numeric(sub('D.+','',sub('.*Q','',doseCycle1)))
    nlag1[is.na(nlag1)]<-0
    
    # Obtain the number of doses of the first treatment schedule
    ndoses1<- as.numeric(sub('.*Q.+D.','',doseCycle1))
    ndoses1[is.na(ndoses1)]<-1
    
    # Time lag between doses of the second treatment schedule
    nlag2<- as.numeric(sub('Dx.+','',sub(c('.*Q'),'', schedule2)))
    nlag2[is.na(nlag2)]<-0
    # Number of doses of the second treatment schedule
    ndoses2<-as.numeric(sub('.*Q.+D.','',schedule2))
    ndoses2[is.na(ndoses2)]<-0
    
    # All all to the dataframe
    treat.df[!sd_indx,ncycles_j]<-ncycles1
    treat.df[!sd_indx,nlag_j] <-nlag1
    treat.df[!sd_indx,ndoses_j]<-ndoses1
    treat.df[!sd_indx,sch2_nlag_j]<-nlag2
    treat.df[!sd_indx,sch2_ndoses_j]<-ndoses2
    
    treat.df[sd_indx,ncycles_j]<-1
    treat.df[sd_indx,nlag_j] <-0
    treat.df[sd_indx,ndoses_j]<-1
    treat.df[sd_indx,sch2_nlag_j]<-0
    treat.df[sd_indx,sch2_ndoses_j]<-0
    
    # Obtain sequence of dosing times for each individual according to the number of schedules, cycles, time lags and number of doses
    dosingTimes<-as.list(by(treat.df,treat.df$id,function(x)dosingTimes.fx(x,nlag_j,ndoses_j,sch2_nlag_j,sch2_ndoses_j,treatStart_j),simplify = F))
    rowNames<-names(dosingTimes)
    # Transform to matrix with one row for each individual 
    dosingTimes<-matrix((dosingTimes),byrow=T,nrow=length(dosingTimes))
    # Assing col and row names
    colnames(dosingTimes)<-dosingTimes_j
    rownames(dosingTimes)<-rowNames
    # Obtain the length of each row element in the matrix
    lengths<-sapply(dosingTimes[,1],length)
    # Create data frame by expanding the id factor by lengths and unlisting the elemnts of the matrix
    dosingTimes.df<-data.frame(id=rep(rownames(dosingTimes),lengths),
                               apply(dosingTimes,2,unlist),row.names=NULL)
    # Match the dosing times with the event/censoring time and dose given by the id
    treat_j.df<-join(treat.df[,c('id', dose_j,'Time','BASELINE_days')],dosingTimes.df,by=c('id'))
    # Remove rows where doses are given after the event
    treat_j.df<-treat_j.df[treat_j.df$Time>=treat_j.df[,dosingTimes_j],]
    # Set to 0 the doses that occur when the event happens
    treat_j.df[,dose_j][treat_j.df$Time==treat_j.df[,dosingTimes_j]]<-0 
    # Set to 0 the doses that occur when dosing times ==0 and BASELINE_days!=0
    treat_j.df[,dose_j][treat_j.df$dosingTimes==0 & treat_j.df$BASELINE_days!=0 ]<-0 
    treat_j.df<-subset(treat_j.df,select=-c(BASELINE_days))
    # Add a factor to identify which treatment it refers to
    treat_j.df$TREAT<-factor(j,levels = 1:3)
    }
    else{ treat_j.df<-data.frame('id'=NA,'DOSE'=NA,'Time'=NA,'dosingTimes'=NA,'TREAT'=NA)}

    names(treat_j.df)<-c('id', 'DOSE', 'Time', 'dosingTimes','TREAT')
    dosingTimes.ls[[j]]<-treat_j.df
    
  }
  # Row bind all the rows of the three dfs that contain the dosing times of individuals receiving 1,2, or 3 treatments
  dosingTimes.df<-rbind(dosingTimes.ls[[1]],dosingTimes.ls[[2]],dosingTimes.ls[[3]],dosingTimes.ls[[4]])
  id<-(as.numeric(dosingTimes.df$id))
  dosingTimes.df <-dosingTimes.df[order(dosingTimes.df$id,dosingTimes.df$dosingTimes),]
  
  # Identify the unique dosing time points for each individual over all treatments
  uniqueTimes.ls<-by(dosingTimes.df,dosingTimes.df$id,function(x)sort(unique(x$dosingTimes)))
  # Extracting id's from names of elements
  rowNames<-names(uniqueTimes.ls)
  # Convert to matrix with one row for each id
  uniqueDosingTimes<-matrix((uniqueTimes.ls),byrow=T,nrow=length(uniqueTimes.ls))
  colnames(uniqueDosingTimes)<-'dosingTimes'
  rownames(uniqueDosingTimes)<-rowNames
  # Extract length of each individual unique dosing times
  lengths<-sapply(uniqueDosingTimes[,1],length)
  # 
  uniqueDosingTimes.df1<-data.frame(id=rep(rowNames,lengths),
                             apply(uniqueDosingTimes,2,unlist),row.names=NULL)
  # Change names to allow for joining of datasets
  names(dosingTimes.ls[[1]])<-c('id', 'DOSE_1', 'Time', 'dosingTimes', 'TREAT1')
  names(dosingTimes.ls[[2]])<-c('id', 'DOSE_2', 'Time', 'dosingTimes', 'TREAT2')
  names(dosingTimes.ls[[3]])<-c('id', 'DOSE_3', 'Time', 'dosingTimes', 'TREAT3')
  names(dosingTimes.ls[[4]])<-c('id', 'DOSE_4', 'Time', 'dosingTimes', 'TREAT4')
  
  # Join datasets by id and unique dosing times
  uniqueDosingTimes.df1<-join(uniqueDosingTimes.df1,dosingTimes.ls[[1]][,c('id','Time','dosingTimes','DOSE_1')], by=c('id','dosingTimes'), type='left')
  uniqueDosingTimes.df1<-join(uniqueDosingTimes.df1,dosingTimes.ls[[2]][,c('id','Time','dosingTimes', 'DOSE_2')], by=c('id','dosingTimes'), type='left')
  uniqueDosingTimes.df1<-join(uniqueDosingTimes.df1,dosingTimes.ls[[3]][,c('id','Time','dosingTimes', 'DOSE_3')], by=c('id','dosingTimes'), type='left')
  uniqueDosingTimes.df1<-join(uniqueDosingTimes.df1,dosingTimes.ls[[4]][,c('id','Time','dosingTimes', 'DOSE_4')], by=c('id','dosingTimes'), type='left')
  
    # # Join the datasets from the treatments that received only one dose
  # uniqueDosingTimes.df1<-join(uniqueDosingTimes.df1,df[,c('id','dosingTimes','Time', 'DOSE_1')], by=c('id','dosingTimes'), type='full')
  # uniqueDosingTimes.df1<-join(uniqueDosingTimes.df1,df[,c('id','dosingTimes','Time', 'DOSE_2')], by=c('id','dosingTimes'), type='full')
  # uniqueDosingTimes.df1<-join(uniqueDosingTimes.df1,df[,c('id','dosingTimes','Time', 'DOSE_3')], by=c('id','dosingTimes'), type='full')
  # 
  names(uniqueDosingTimes.df1)<-c('id','dosingTimes','Time1', 'DOSE_1','Time2','DOSE_2', 'Time3','DOSE_3','Time4','DOSE_4')
  
  # Create a unique Time vector of event/censoring times 
  uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$Time1),'Time1']<-uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$Time1),'Time2']
  uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$Time1),'Time1']<-uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$Time1),'Time3']
  uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$Time1),'Time1']<-uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$Time1),'Time4']
  
  # Set dose to 0 wherever there is a NA
  uniqueDosingTimes.df1<-uniqueDosingTimes.df1[,c('id','dosingTimes', 'Time1','DOSE_1','DOSE_2','DOSE_3')]
  uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$DOSE_1),'DOSE_1']<-0
  uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$DOSE_2),'DOSE_2']<-0
  uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$DOSE_3),'DOSE_3']<-0
  uniqueDosingTimes.df1[is.na(uniqueDosingTimes.df1$DOSE_4),'DOSE_4']<-0
  
  # Calculate the cumulative dose per individual
  id<-(as.numeric(uniqueDosingTimes.df1$id))
  uniqueDosingTimes.df1$id<-id
  uniqueDosingTimes.df1 <-uniqueDosingTimes.df1[order(uniqueDosingTimes.df1$id,uniqueDosingTimes.df1$dosingTimes),]
  cumdose1 <-by(uniqueDosingTimes.df1,uniqueDosingTimes.df1$id,function(x)cumsum(x[,'DOSE_1']))
  cumdose1 <-matrix(unlist(cumdose1))
  cumdose2 <-by(uniqueDosingTimes.df1,uniqueDosingTimes.df1$id,function(x)cumsum(x[,'DOSE_2']))
  cumdose2 <-matrix(unlist(cumdose2))
  cumdose3 <-by(uniqueDosingTimes.df1,uniqueDosingTimes.df1$id,function(x)cumsum(x[,'DOSE_3']))
  cumdose3 <-matrix(unlist(cumdose3))
  cumdose4 <-by(uniqueDosingTimes.df1,uniqueDosingTimes.df1$id,function(x)cumsum(x[,'DOSE_4']))
  cumdose4 <-matrix(unlist(cumdose4))
  
  uniqueDosingTimes.df1[,'CUMULATIVE_DOSE1']<-cumdose1  
  uniqueDosingTimes.df1[,'CUMULATIVE_DOSE2']<-cumdose2  
  uniqueDosingTimes.df1[,'CUMULATIVE_DOSE3']<-cumdose3  
  uniqueDosingTimes.df1[,'CUMULATIVE_DOSE4']<-cumdose4  
  
  varNames<-names(uniqueDosingTimes.df1)
  
  df1<-join(uniqueDosingTimes.df1[,!is.element(varNames,'Time1')],df,by=c('id'),type='full')
  # df1$Event<-df1$dosingTimes==df1$Time
  df1 <-df1[order(df1$id,df1$dosingTimes),]
  df1$CUMULATIVE_DOSE1[df1$Treatment_1=='Control']<-0
  df1$CUMULATIVE_DOSE2[is.na(df1$Treatment_2)]<-0
  df1$CUMULATIVE_DOSE3[is.na(df1$Treatment_3)]<-0
  df1$CUMULATIVE_DOSE4[is.na(df1$Treatment_4)]<-0
  
  df1$dosingTimes[df1$Treatment_1=='Control']<-0
  
  df1$Event<-(df1$Time==df1$dosingTimes & df1$Event | df1$Treatment_1=='Control')
  
  return(df1) 
  
}
dosingTimes.fx<-function(df1,nlag_j,ndoses_j,sch2_nlag_j,sch2_ndoses_j,treatStart_j){
    
    singlecycle1 <-seq(0,(df1[1,ndoses_j]-1)*df1[1,nlag_j], df1[1,nlag_j])

    cycles1<-rep(singlecycle1,df1[1,'ncycles1'])
    cyclestart1<-rep(seq(0,7*(df1[1,'ncycles1']-1),7),each=df1[1,ndoses_j])

    dosingTimes1<-cycles1+cyclestart1
    
    singlecycle2<-seq(0,(df1[1,sch2_ndoses_j]-1)*df1[1,sch2_nlag_j],df1[1,sch2_nlag_j])
    cyclestart2<-dosingTimes1[length(dosingTimes1)]+14
    if(df1[1,sch2_ndoses_j]>0) dosingTimes2<-singlecycle2+cyclestart2 else dosingTimes2<-singlecycle2
    
    dosingTimes=sort(unique(c(0,dosingTimes1+sum(df1[1,c('BASELINE_days',treatStart_j)]),dosingTimes2+sum(df1[1,c('BASELINE_days',treatStart_j)]),df1[1,'Time'])))
    return(dosingTimes)
    
}
plotForest<-function(df,treat){
  re.est<-df$beta[1]
  re.lb<-df$ci.lb
  re.ub<-df$ci.ub
  y<-df$yi.f[1:length(df$yi.f)]
  v<-df$vi.f[1:length(df$vi.f)]
  y.names<-attr(df$yi.f,'slab')
  
  y.ub<-y+sqrt(v)*1.96
  y.lb <-y-sqrt(v)*1.96
  redf<-data.frame('Name'='Overall Estimate',"RE.est"=re.est,
                   'LB'=c(re.lb),'UB'=c(re.ub),'Class'=1,'order'='Overall Estimate')
  plotdf<-data.frame('Name'=c(y.names),"RE.est"=c(y),
                     'LB'=c(y.lb),'UB'=c(y.ub),'Class'=0)
  plotdf$order<-with(plotdf,reorder(Name,RE.est,median))
  plotdf<-rbind(plotdf,redf)
  
  p<-ggplot(plotdf,aes(y=order,x=exp(RE.est),xmin=exp(LB),xmax=exp(UB)))+
    geom_point()+
    geom_errorbarh(height=.1)+
    labs(title=("Individual effect sizes vs overall efect"),subtitle = treat)+
    xlab('log-MR estimates') +
    ylab('')+
    theme_minimal()+
    theme(text=element_text(family="Helvetica",size=8, color="black"),
          title=element_text(family='Helvetica',size=12,color='black',face='bold'))+
    theme(panel.spacing = unit(1, "lines"))+
    geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5)+
    geom_vline(xintercept=exp(re.est), color="red", linetype="dashed", alpha=1)+
    facet_grid(Class~., scales= "free", space="free")
  
  return(p)
}

modelComparison<-function(models.list,models.form){
  
  AIC.models<-sapply(models.list,function(x)x$fit.stats['AIC','REML'])
  I2.models <-sapply(models.list,function(x)x$I2)
  QEp.models <-sapply(models.list,function(x)x$QEp)
  QMp.models <-sapply(models.list,function(x)x$QMp)
  R2<-(sapply(models.list,function(x)x$R2))
  tau <-sapply(models.list,function(x)x$tau2)
  intercept <-sapply(models.list,function(x)x$beta[1])
  pvalue <-sapply(models.list,function(x)x$pval[1])
  if(is.null(R2)){R2<-rep(0,length(AIC.models))}
  m.formulas<-sapply(models.form,function(x)x)
  HR.M.summ<-data.frame('Formula'=c(as.character(m.formulas)),
                         'AIC'=round(AIC.models,4),'R^2'=round(R2,2),'I^2_perc'=round(I2.models,2),'tau'=round(tau,4),
                        'intercept'=round(intercept,4),
                         'Heterogeneity test p-value'=round(QEp.models,4),
                        'Modifiers test p-value'=round(QMp.models,4),
                        'intercept p-value'=round(pvalue,4))
  return(HR.M.summ)
  
}

plotTF<-function(tf.df){
  m0tf.plot<-ggplot(tf.df,aes(x=points,y=1/(se)))+
    geom_point(aes(color=Fill.in))+
    geom_vline(aes(xintercept =Estimate,linetype='Estimate'), data=tf.df)+
    geom_vline(aes(xintercept =Adjusted.Estimate,linetype='Adjusted Estimate'), data=tf.df,color='red')+
    scale_linetype_manual(name = 'Estimates',
                          values = c('Estimate' = 2,
                                     'Adjusted Estimate' = 2),
                          guide = guide_legend(override.aes = list(colour = c('red',
                                                                              'black'))))+
    labs(title='Trim-and-fill analysis for log-Hazard ratios ',subtitle = ' All ICB mAb Treatments')+
    xlab('log(HR)')+ ylab(" Precision (1/SE)") +
    theme_minimal()+
    theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
          axis.text.x=element_text(face="plain", size=8,family="Helvetica"),
          axis.title=element_text(size=12,face="bold",family="Helvetica"),
          strip.text.x=element_text(size=6))+ 
    scale_color_grey()+
    facet_wrap(~Treatment)
  return(m0tf.plot)
  
  
}
plotTFbyCell<-function(tf.df,factor.var){
  m0tf.plot<-ggplot(tf.df,aes(x=points,y=1/(se)))+
    geom_point(aes(color=Fill.in))+
    geom_vline(aes(xintercept =Estimate,linetype='Estimate'), data=tf.df)+
    geom_vline(aes(xintercept =Adjusted.Estimate,linetype='Adjusted Estimate'), data=tf.df,color='red')+
    scale_linetype_manual(name = 'Estimates',
                          values = c('Estimate' = 2,
                                     'Adjusted Estimate' = 2),
                          guide = guide_legend(override.aes = list(colour = c('red',
                                                                              'black'))))+
    labs(title='Trim-and-fill analysis for log-Hazard ratios ',subtitle = tf.df$Treatment[1])+
    xlab('log(HR)')+ ylab(" Precision (1/SE)") +
    theme_minimal()+
    theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
          axis.text.x=element_text(face="plain", size=8,family="Helvetica"),
          axis.title=element_text(size=12,face="bold",family="Helvetica"))+ 
    scale_color_grey()+
    facet_wrap(as.formula(paste('~', factor.var)))
  return(m0tf.plot)
}
plotTFbyLab<-function(tf.df){
  m0tf.plot<-ggplot(tf.df,aes(x=points,y=1/(se)))+
    geom_point(aes(color=Fill.in))+
    geom_vline(aes(xintercept =Estimate,linetype='Estimate'), data=tf.df)+
    geom_vline(aes(xintercept =Adjusted.Estimate,linetype='Adjusted Estimate'), data=tf.df,color='red')+
    scale_linetype_manual(name = 'Estimates',
                          values = c('Estimate' = 2,
                                     'Adjusted Estimate' = 2),
                          guide = guide_legend(override.aes = list(colour = c('red',
                                                                              'black'))))+
    labs(title='Trim-and-fill analysis for log-Hazard ratios ',subtitle = tf.df$Treatment[1])+
    xlab('log(HR)')+ ylab(" Precision (1/SE)") +
    theme_minimal()+
    theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
          axis.text.x=element_text(face="plain", size=8,family="Helvetica"),
          axis.title=element_text(size=12,face="bold",family="Helvetica"))+ 
    scale_color_grey()+
    facet_wrap(~Lab)
  return(m0tf.plot)
  
  
}

getTFdf<-function(m0,m0.tf){
  m0.points<-(m0.tf$yi.f[1:length(m0.tf$yi.f)])
  n_fill<-sum(m0.tf$fill)
  n_obs<-sum(!m0.tf$fill)
  m0.fillpoints<-m0.points[(n_obs+1):n_fill]
  m0.se<-sqrt(m0.tf$vi[1:length(m0.tf$vi)])
  m0.fill <-rep(c('Observed Data'),each=c(n_obs))
  m0.fill1<-rep(c('Fill-in data'),n_fill)
  
  fill.character<-c(m0.fill,m0.fill1)
  m0tf.df<-data.frame('Estimate'=m0$beta,'Adjusted Estimate'=m0.tf$beta,'points'=m0.points,'se'=m0.se,'Fill-in'=fill.character)
  
}

estimateComparison <-function(mcomp,mcomp.form,treats){
  HRest.df<-data.frame('Treatment'=NA,'Estimate'=0,'LB'=0,'UB'= 0, 'Model'=NA)
  hrindx<-1
  for(j in 1:length(mcomp)){
  mi.rownames<-rownames(mcomp[[j]]$beta)
    for(i in 1:length(treats)){
      namesIndx<-is.element(mi.rownames,paste('Treatments',treats[i],sep=''))
      if(all(!namesIndx)){
        namesIndx<-is.element(mi.rownames,'intrcpt')
        intrcptOnly<-T
      }
      else{intrcptOnly<-F}
      HRintercept<-mcomp[[j]]$beta['intrcpt',1]
      HRintercept.se<-(mcomp[[j]]$vb['intrcpt','intrcpt'])
      
      HRtreat<-mcomp[[j]]$beta[namesIndx,1]
      HRtreat.se<-(mcomp[[j]]$vb[namesIndx,namesIndx])
      if(intrcptOnly){
        HRestimate<-HRintercept
        HR.se<-HRintercept.se
      }
      else
      HRestimate <- HRintercept+HRtreat
      HR.se<-sqrt(HRintercept.se + HRtreat.se)
      znorm.95<-qnorm(.975)
      HRest.df[hrindx,c('Estimate','LB','UB')]<-c(HRestimate, HRestimate-znorm.95*HR.se,HRestimate+znorm.95*HR.se)
      HRest.df[hrindx,c('Treatment','Model')]<-c(treats[i],mcomp.form[j])
      hrindx<-hrindx+1
    }
  
  }
  return(HRest.df)
}
plotEstimateComparison<-function(HRest.df){
  # Plot estimates comparisons
  est.p<-ggplot(HRest.df,aes(x=Model,y=exp(Estimate),ymin=exp(LB),ymax=exp(UB)))+
    geom_point()+
    labs(title='Estimates of treatment effects',subtitle = 'Model comparison')+
    xlab('')+ ylab("Hazard Ratio (95% Confidence Interval)")+
    geom_errorbar(aes(ymin=exp(LB), ymax=exp(UB),col=Model),width=0.2,cex=1)+ 
    geom_hline(yintercept = 1,color='red',linetype=2)+
    facet_wrap(~Treatment) +
    theme_minimal()+
    theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
          axis.title=element_text(size=12,face="bold",family="Helvetica"),
          axis.text.x=element_blank(),
          axis.ticks=element_blank())+
    scale_color_grey()+
    theme(panel.spacing = unit(.5, "lines"))
    
  return(est.p)
}

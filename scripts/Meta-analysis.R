## Load packages
library(MCMCglmm)
library(lattice)
library(ggplot2)
library(plyr)
## set working directory
wd <- getwd()

# Load dataset
survival<-read.csv('~/Documents/Data/CSV/Survival/Survival_Preclinical.csv',
                   header = TRUE, sep = ';', quote = "\"", dec = ',')
survival_clinical<-read.csv('~/Documents/Data/CSV/Survival/Survival_Clinical.csv',
                            header = TRUE, sep = ';', quote = "\"", dec = ',')

## Set default levels in each of the categories
survival<-within(survival, MAB<-relevel(MAB, 'Control'))
survival<-within(survival, STRAIN<-relevel(STRAIN,'C57BL/6'))
survival<-within(survival, CELL<-relevel(CELL,'B16F10'))

## Subset relevant treatments
therapy_indx <- survival$MAB=='antiPDL1'|survival$MAB=='antiPD1'|survival$MAB=='antiCTLA4'|
  survival$MAB=='antiCTLA4+antiPD1'|survival$MAB=='antiCTLA4+antiPDL1'|survival$MAB=='Control'|
  survival$MAB=='iIDO'|survival$MAB=='antiCTLA4+iIDO'|survival$MAB=='antiPDL1+iIDO'|survival$MAB=='antiCTLA4+antiPDL1+iIDO'

survival_subset <-subset(survival, is.finite(survival$MOS_COR)&therapy_indx, 
                         select = c(STUDY_ID,MOS_COR, CELL, CANCER, STRAIN, MAB, LAB,N, EUTHANIZED,
                                    DOSE_mg.kg, DOSE_2, DOSE_3,T_0_cells))

## Drop irrelevant levels from each category
survival_subset$MAB<-droplevels(survival_subset$MAB)
survival_subset$CELL<-droplevels(survival_subset$CELL)
survival_subset$LAB<-droplevels(survival_subset$LAB)

## Check levels
levels(survival_subset$MAB)
levels(survival_subset$LAB)
levels(survival_subset$CELL)
levels(survival_subset$CANCER)

## Model proposal and evaluation
# Model 0: Fixed effect:MAB
prior0 = list(R = list(V = 1, nu = .002))

model0_glm <- MCMCglmm(log(MOS_COR)~MAB, family = "gaussian", data = survival_subset,
			prior=prior0,pr=TRUE,nitt=1e6, burnin=1e4, thin=100, singular.ok = TRUE)


# Model 1: Fixed effect:MAB+Random effect:CELL
prior1 = list(R = list(V = 1, nu = .002), G = list(G1 = list(V = 1,
     		nu =.02)))
model1_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL, family = "gaussian", data = survival_subset,
			prior=prior1,pr=TRUE,nitt=1e6, burnin=1e4, thin=100, singular.ok = TRUE)

# Model 2: Fixed effect:MAB+Random effect:LAB
prior2 = list(R = list(V = 1, nu = .002), G = list(G1 = list(V = 1,
     		nu =.02)))

model2_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~LAB, family = "gaussian", data = survival_subset,
			prior=prior2,pr=TRUE,nitt=1e6, burnin=1e4, thin=100, singular.ok = TRUE)


# Model 3: Fixed effect:MAB+Random effect:CELL+LAB
prior3 = list(R = list(V = 1, nu = 2), G = list(G1 = list(V = 1,
     		nu =2), G2 = list(V = 1, nu = 2)))
model3_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL+LAB, data=survival_subset,
				prior = prior3,pr=TRUE,nitt=1e6,burnin=1e4, thin=1e2)
model3_glm_red <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL+LAB, data=survival_subset,
                       prior = prior3,pr=FALSE,nitt=1e4,burnin=1e3, thin=10)
sum3 <- summary(model3_glm_red)

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

# Model 8: Fixed effect:MAB+DOSE+RE:(1\CELL))
prior8 = list(R = list(V = diag(1), nu = 0.2),
              G = list(G1 = list(V = diag(6), nu = 0.2)))
model8_glm <- MCMCglmm(log(MOS_COR)~MAB, random=~CELL, data=survival_subset,
                      nitt=1e6,burnin=1e4, thin=1e2)

## Comparing all models
DIC <-c(summary(model0_glm)$DIC,summary(model1_glm)$DIC,summary(model2_glm)$DIC,
		summary(model3_glm)$DIC,summary(model4_glm)$DIC,summary(model5_glm)$DIC,
		summary(model6_glm)$DIC,summary(model7_glm)$DIC)
data.frame(DIC,row.names=c('Model0', 'Model1', 'Model2', 'Model3', 'Model4', 'Model5', 'Model6', 'Model7'  ))

## Model selection
model_final<-model3_glm

## Model identifiability
par(mfrow=c(10,2), mar=c(2,2,1,0))
plot(model_final$Sol, auto.layout=F)

## Model estimates of HR
# Posterior prediction intervals of the proportional increases in median survival
mono_indx<-c(2, 7, 8, 10) # Column indexes where monotherapies are located
combo_indx<-c(3, 4, 5,6, 9)
HR_HPD_3_Monotherapy<-HPDinterval(exp(-(model3_glm$Sol[,mono_indx])))
HR_modes_3_Monotherapy<-posterior.mode(exp(-(model3_glm$Sol[,mono_indx])))

HR_HPD_3_Combination<-HPDinterval(exp(-(model3_glm$Sol[,combo_indx]-model3_glm$Sol[,mono_indx[1]]))) # set HR relative to antiCTLA4
HR_HPD_3_Combination[c(5,10)]<-HPDinterval(exp(-(model3_glm$Sol[,9]-model3_glm$Sol[,8])))# Set antiPDL1+iIDO relative to antiPDL1

HR_modes_3_Combination<-posterior.mode(exp(-(model3_glm$Sol[,combo_indx]-model3_glm$Sol[,mono_indx[1]])))
HR_modes_3_Combination[5]<-posterior.mode(exp(-(model3_glm$Sol[,9]-model3_glm$Sol[,8])))

HR_CI_3<-cbind(HR_modes_3_Monotherapy,HR_HPD_3_Monotherapy)
HR_CI_3<-rbind(HR_CI_3,cbind(HR_modes_3_Combination,HR_HPD_3_Combination))

# Proportion of variance explained by each term
fixed_var<-as.mcmc(diag((var(t(as.matrix(model_final$Sol[,2:10])))))) # sigma_f variation in fixed effects responses

var_comp<-as.mcmc(cbind(fixed_var,model_final$VCV[,1:3])) # MCMC samples of variance for fixed, random and error variances
posterior.mode(var_comp/rowSums(var_comp))
HPDinterval(var_comp/rowSums(var_comp),prob=.95) # Credible intervals

## Model predictions
survival_subset$fit<-predict(model_final,type ="response",posterior='all')
survival_subset$residuals<- survival_subset$predict-log(survival_subset$MOS_COR)
survival_subset$prediction<-predict(model_final,marginal=~LAB+CELL,type ="response",interval = 'prediction')

# Data vs Residuals 
data_res<-ggplot(survival_subset, aes(x=log(MOS_COR),y=residuals,colour=CELL))+
  geom_point()+
  facet_wrap(~LAB)
data_res
# Data vs fit
data_pred<-ggplot(survival_subset, aes(x=log(MOS_COR),y=fit,colour=CELL))+
  geom_point()+
  facet_wrap(~LAB)
data_pred+ geom_abline(intercept = 0, slope = 1)

# Data vs predictions
data_pred<-ggplot(survival_subset, aes(x=log(MOS_COR),y=prediction))+
  geom_point()+
  facet_wrap(~LAB)
data_pred+ geom_abline(intercept = 0, slope = 1)

# Plot predictions against clinical data
survival_clinical_subset=subset(survival_clinical,is.finite(OS_HR),
                                select = c(STUDY_ID,OS_HR, OS_HR_LB,
                                           OS_HR_UB, PFS_HR,PFS_HR_LB,
                                           PFS_HR_UB,MAB,CANCER,AGENT,DOSE,N_mAb))
# HR of clinical studies
OS_HR<-ggplot(survival_clinical_subset,aes(AGENT,OS_HR,colour=CANCER))+
  geom_errorbar(aes(ymin=OS_HR_LB, ymax=OS_HR_UB))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 45))

OS_HR+labs(color='Indication',  y='Hazard Ratio (Clinical)',x="Study",
        title='Survival in clinical studies',
        subtitle='Overall survival')

PFS_HR<-ggplot(survival_clinical_subset,aes(AGENT,PFS_HR,colour=CANCER))+ 
  geom_errorbar(aes(ymin=PFS_HR_LB, ymax=PFS_HR_UB))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 45))

PFS_HR+labs(color='Indication',  y='Hazard Ratio (Clinical)',x="Study",
         title='Survival in clinical studies',
         subtitle='Progression-free survival')

# Add predictions to dataset
HR_MAB<-sub("MAB","",rownames(HR_CI_3))
preclin.pred<- data.frame(HR_CI_3,row.names = HR_MAB)
for(i in 1:dim(survival_clinical_subset)[1]) {
  # Compare clinical  with preclinical treatment
  
  for (j in 1:dim(HR_CI_3)[1]){
    if(survival_clinical_subset$MAB[i]==row.names(preclin.pred)[j]){
      survival_clinical_subset$OS_HR_Pred[i]=preclin.pred[j,1]
    survival_clinical_subset$OS_HR_Pred_LB[i]=preclin.pred[j,2]
    survival_clinical_subset$OS_HR_Pred_UB[i]=preclin.pred[j,3]
    survival_clinical_subset$Preclinical_MAB[i]=HR_MAB[j]
    }
  }
  if(survival_clinical_subset$MAB[i]=="antiPD1+iIDO"){
    survival_clinical_subset$OS_HR_Pred[i]=preclin.pred[HR_MAB=="antiPDL1+iIDO",1]
    survival_clinical_subset$OS_HR_Pred_LB[i]=preclin.pred[HR_MAB=="antiPDL1+iIDO",2]
    survival_clinical_subset$OS_HR_Pred_UB[i]=preclin.pred[HR_MAB=="antiPDL1+iIDO",3]
    survival_clinical_subset$Preclinical_MAB[i]="antiPDL1+iIDO"
    
  }
}
# Plot OS vs PFS

OS_PFS<-ggplot(survival_clinical_subset,aes(PFS_HR,OS_HR,colour=CANCER))+
  geom_point()+
  geom_errorbar(aes(ymin=OS_HR_LB, ymax=OS_HR_UB)) +
  geom_errorbarh(aes(xmin=PFS_HR_LB,xmax=PFS_HR_UB))+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept=1)+
  geom_abline(linetype=2)
OS_PFS+labs(color='Indication',  y='OS Hazard Ratio (Clinical)',x="PFS Hazard Ratio (Clinical)",
        title='Survival in clinical trials of ICB',
        subtitle='Progression-free vs overall survival')

# Plot OS HR 
HR<-ggplot(survival_clinical_subset,aes(OS_HR_Pred,OS_HR,colour=CANCER))+ 
  geom_point()+
  geom_errorbar(aes(ymin=OS_HR_LB, ymax=OS_HR_UB)) +
  geom_errorbarh(aes(xmin=OS_HR_Pred_LB,xmax=OS_HR_Pred_UB))+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept=1)+
  geom_abline(linetype=2)+
  coord_trans(x="log")
HR+labs(color='Indication',  y='Hazard Ratio (Clinical)',x="Hazard Ratio (Preclinical)",
    title='Survival in preclinical vs clinical settings',
    subtitle='Overall survival')

# Plot PFS HR 
PFS_HR<-ggplot(survival_clinical_subset,aes(OS_HR_Pred,PFS_HR,colour=AGENT))+ 
  geom_point()+
  geom_errorbar(aes(ymin=PFS_HR_LB, ymax=PFS_HR_UB)) +
  geom_errorbarh(aes(xmin=OS_HR_Pred_LB,xmax=OS_HR_Pred_UB))+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept=1)+
  geom_abline(linetype=2)+
  coord_trans(x="log")

PFS_HR+labs(color='Treatment',  y='Hazard Ratio (Clinical)',x="Hazard Ratio (Preclinical)",
        title='Survival in preclinical vs clinical settings',
        subtitle='Progression-free survival')

# Regression of clinical hazard ratios to preclininal ones
trans_model<-lm(OS_HR~OS_HR_Pred,data=survival_clinical_subset)

# Correlation coefficient
cor.test(survival_clinical_subset$PFS_HR[!is.na(survival_clinical_subset$PFS_HR)],survival_clinical_subset$OS_HR_Pred[!is.na(survival_clinical_subset$PFS_HR)])

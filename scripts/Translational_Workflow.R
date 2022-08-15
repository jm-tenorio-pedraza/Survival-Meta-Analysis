# Translational workflow
rm(list=ls())
library(plyr)
library(metafor)
library(ggplot2)
library(latex2exp)
library(openxlsx)
outputFilePath<-'~/Documents/Thesis/Results/'
# Load models based on which vars to study:
clin.var<-'OS_HR'
pre.var<-'MR'
if(clin.var=='OS_HR'){
  load(file=paste(outputFilePath,'OS_HR_clinical.m0.RData',sep = ''))
  load(file=paste(outputFilePath,'OS_HR_clinical.m1.RData',sep = ''))
  
  y.label<-'OS hazard ratios'
  yaxis.label<-'Clincical OS hazard ratios'
} else {load(file=paste(outputFilePath,'PFS_HR_clinical.m0.RData',sep=''))
  load(file=paste(outputFilePath,'PFS_HR_clinical.m1.RData',sep=''))
  y.label<- 'PFS hazard ratios'
  yaxis.label<-'Clinical PFS hazard ratios'
}
if(pre.var=='HR'){
  load(file=paste(outputFilePath,'HR_preclinical.m0.RData',sep=''))
  load(file=paste(outputFilePath,'HR_preclinical.m1.RData',sep=''))
  
  x.label<-'Hazard ratios'
  xaxis.label<-'Preclinical OS hazard ratios'
  
} else {load(file=paste(outputFilePath,'MR_preclinical.m0.RData',sep=''))
  load(file=paste(outputFilePath,'MR_preclinical.m1.RData',sep=''))
  x.label<-'Median survival ratios'
  xaxis.label<-'Preclinical median survival ratios'
}
## Load table with cancer types in rows and cell lines in columns
cancer.mat<-read.table(file=paste(outputFilePath,'Preclinical_',pre.var,'CELL_predMat.csv',sep=''))
# Load table with leave-one-out cv prediction
pred.mat<-read.table(file=paste(outputFilePath,clin.var,'_','Clinicical_LOO_Predictions.m1.csv',sep=''))
# How many cell lines represent each type of cancer
nCellLines<-apply(cancer.mat,1,(function(x)sum(x>0)))
# Load clinical data
clinical.mat<-read.csv('~/Documents/GitHub/Survival-Meta-Analysis/output/clinical.red.df.csv',sep=',')

## Result 1: Preclinical predictions based on treatment and cancer combinations compared to actual data
# Match treatment effect estimates between preclinical and clinical models
preclinical.treatments<-sub('Treatments','',row.names(preclinical.m1$beta))
clinical.treatments<-sub('Treatments','',row.names(clinical.m1$beta))
# Identify which preclin treatments are present in the clinical treats
transIndx<-preclinical.treatments%in%clinical.treatments
treatments<-preclinical.treatments[transIndx]
# Create prediction matrix based on the identified treatments
newmods.mat<- diag(1,nrow=length(transIndx))*as.numeric(transIndx)
rowIndx<-apply(newmods.mat,1,function(x)all(x==0))
newmods.mat<-newmods.mat[!rowIndx,]
# Add the cell-related columns
newmods.mat<-kronecker(newmods.mat,rep(1,dim(cancer.mat)[1])) # Expand matrix of treatments 
newmods.mat[,grep('CELL_FAMILY',preclinical.treatments)]<-kronecker((rep(1,length(treatments))),as.matrix(cancer.mat))
# Generate predictions for each combination of cancer and treatment
preclinical.m1.df<-(predict(preclinical.m1,newmods =newmods.mat,transf = exp))
preclinical.m1.df<-data.frame('Treatments'=rep(treatments,each=dim(cancer.mat)[1]),
                           'CANCER'=rep(row.names(cancer.mat),length(treatments)),
                           'Preclinical.Est'=preclinical.m1.df$pred,'Preclinical.LB'=preclinical.m1.df$ci.lb,
                           'Preclinical.UB'=preclinical.m1.df$ci.ub,
                           'N_CellLines'=rep(nCellLines,length(treatments)))

trans.m1.df<-join(preclinical.m1.df,clinical.mat[,c('STUDY_ID','CANCER','Treatments',clin.var,paste(clin.var,'LB',sep='_'),paste(clin.var,'UB',sep='_'),
                                                    by=c('Treatments','CANCER'))])
trans.m1.df$TreatmentLabels<-factor(trans.m1.df$Treatments,levels=treatments,
                               labels=c('anti-CTLA-4','anti-CTLA-4 + anti-PD-1',
                                        'anti-CTLA-4 + anti-PD-L1', 'anti-CTLA-4 + Chemotherapy',
                                        'anti-PD-1','anti-PD-1 + Chemotherapy','anti-PD-L1',
                                        'anti-PD-L1 + Chemotherapy'))
# Log-transform the OS_HR estimates 
if(clin.var=='OS_HR')trans.m1.df$y<-log(trans.m1.df[,clin.var]) else{
  trans.m1.df$y<-(trans.m1.df[,clin.var])
}
# Remove entries without confidence intervals
trans.m1.df<-trans.m1.df[!is.na(trans.m1.df[,paste(clin.var,'UB',sep='_')]),]
# Calculate squared errors between preclinical prediction and cv prediction
trans.m1.df$SqE<-with(trans.m1.df,((Preclinical.Est)-exp(y))^2)
# Sort by squared errors
transSortIndx<-order((trans.m1.df$SqE))
trans.m1.df<-trans.m1.df[transSortIndx,]
# Save predictions
write.xlsx(trans.m1.df,file=paste(outputFilePath,'Clinical_',clin.var,'_Preclinical_',pre.var,'_Predictions_Model1.xlsx',sep=''))

# Plot the predictions vs clinical dta
pred.m1.plot<-ggplot(trans.m1.df,aes(x=Preclinical.Est,y=exp(y),color=CANCER))+
  geom_point()+
  labs(title='Predictions of clinical efficacy', 
       subtitle=paste('Preclinical model predictions (',pre.var,') vs clinical efficacy (', clin.var,')',sep=''))+
  xlab('Preclinical predictions') +
  ylab('Clinical efficacy')+
  geom_abline(slope=1,intercept = 0,linetype=2,color='black')+
  scale_x_continuous(limits=c(0.3,1.2))+
  scale_y_continuous(limits=c(0.3,1.2))+
  theme_minimal()+
  scale_color_discrete()+
  theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
        axis.title=element_text(size=12,face="bold",family="Helvetica"))
pred.m1.plot
ggsave(pred.m1.plot, file=paste(outputFilePath,'Model1_Predictions_vs_Data_Preclinical_', pre.var, '_vs_Clinical_', clin.var,'.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')

# Plot by cancer type and treatment ( omit large errors to show distribution)
(trans.m1.plot<-ggplot(trans.m1.df[1:dim(trans.m1.df)[1],],aes(x=TreatmentLabels, y=exp(y)-(Preclinical.Est),color=CANCER))+
  geom_point()+
  labs(color='Cancer',  size='Number of cell lines',title='Prediction errors for clinical efficacy estimates based on preclinical model',
       subtitle=paste('Preclinical', x.label, 'vs', 'Clinical',y.label,sep=' '),
                       y='Data - prediction')+
    ylab(TeX('y - $\\hat{y}$ (Preclinical model prediction errors)'))+
    xlab('Treatments')+
  theme_minimal()+
  scale_color_discrete()+
  theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
        axis.title=element_text(size=12,face="bold",family="Helvetica"),
        axis.text.x=element_text(angle=45,vjust=.5)))
ggsave(trans.m1.plot, file=paste(outputFilePath,'Model1_Error_Preclinical_', pre.var, '_vs_Clinical_', clin.var,'.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')

# Calculate average SqE by cancer
sort(tapply(trans.m1.df$SqE,trans.m1.df$CANCER,mean))
f#Calculate average SqE by treatment
sort(tapply(trans.m1.df$SqE,trans.m1.df$Treatments,median))

# Number of discordant pairs based on upper interval
predTab<-table((trans.m1.df$Preclinical.UB<1)*1,trans.m1.df[,paste(clin.var,'UB',sep='_')]<1)
attributes(predTab)$dimnames[[1]]<-c('Prediction: Not significant', 'Prediction: Significant')
attributes(predTab)$dimnames[[2]]<-c('Data: Not significant', 'Data: Significant')
predTab
write.xlsx(predTab,file=paste(outputFilePath,'Clinical_',clin.var,'_Preclinical_',pre.var,'ConfusionMat.xlsx',sep=''))

# Compare predictions from preclinical model to loo predicted values from clinical model
clinical.m1.df<-pred.mat
clinical.m1.df$Study_Treatments<-with(clinical.m1.df,paste(Study,Treatments,sep='_'))
trans.m1.df$Study_Treatments<-paste(trans.m1.df$STUDY_ID,trans.m1.df$Treatments,sep='_')
trans.m1.df<-join(trans.m1.df,clinical.m1.df[,c('Study_Treatments','Prediction')],by='Study_Treatments')
trans.m1.df<-trans.m1.df[!is.na(trans.m1.df$Prediction),]

# Save unique cancer-type specific treatment estimates
trans.m1.df$Treat_Cancer<-paste(trans.m1.df$Treatments,trans.m1.df$CANCER,sep='_')
trans.m1.df.sub<-trans.m1.df[match(unique(trans.m1.df$Treat_Cancer),trans.m1.df$Treat_Cancer),]
write.xlsx(trans.m1.df.sub,file=paste(outputFilePath,'Clinical_',clin.var,'_Preclinical_',pre.var,'CancerTypeUniquePred.xlsx',sep=''))

# Plotting fitting errors to prediction errors
(errors.m1.plot<-ggplot(trans.m1.df[1:dim(trans.m1.df)[1],],aes(x=exp(y)-exp(Prediction),y=exp(y)-Preclinical.Est,color=TreatmentLabels))+
  geom_point()+
  labs(color='Treatment', title='LOO errors vs prediction errors',
       subtitle=paste('Preclinical', x.label, 'vs', 'Clinical',y.label,sep=' '))+
    xlab(TeX('y - $\\hat{y}$ (Clinical fitting errors)'))+
    ylab(TeX('y - $\\hat{y}$ (Preclinical to clinical prediction errors)'))+
    theme_minimal()+
  scale_color_discrete()+
    geom_abline(slope=1,intercept=0,size=.3)+
  theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
        axis.title=element_text(size=12,face="bold",family="Helvetica"),
        axis.text.x=element_text(angle=45,vjust=.5)))
ggsave(errors.m1.plot, file=paste(outputFilePath,'Model1_Prediction_Error_Preclinical_', pre.var, '_vs_LOO_Errors_Clinical_', clin.var,'.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')
mse.clin<-with(trans.m1.df,crossprod(exp(y)-exp(Prediction))/dim(trans.m1.df)[1])
mse.preclin<-with(trans.m1.df,(crossprod((Preclinical.Est)-exp(y)))/dim(trans.m1.df)[1])

## Result 2: Preclinical model estimates compared to clinical model estimates for treatment effects only
# Generate predictions for overall treatment effect with preclin model to compare to the overall clinical effect
preclinical.treatments<-sub('Treatments','',row.names(preclinical.m0$beta))
clinical.treatments<-sub('Treatments','',row.names(clinical.m0$beta))
# Identify which preclin treatments are present in the clinical treats
transIndx<-preclinical.treatments%in%clinical.treatments
treatments<-preclinical.treatments[transIndx]
newmods.mat<- diag(1,nrow=length(preclinical.m0$beta))*as.numeric(transIndx)
rowIndx<-apply(newmods.mat,1,function(x)all(x==0))
newmods.mat<-newmods.mat[!rowIndx,]
preclinical.m0.df<-(predict(preclinical.m0,newmods =newmods.mat,vcov=TRUE))
preclinical.m0.df<-data.frame('Treatments'=treatments,
                            'Preclinical.Est'=preclinical.m0.df$pred$pred,'Preclinical.LB'=preclinical.m0.df$pred$ci.lb,
                            'Preclinical.UB'=preclinical.m0.df$pred$ci.ub,'Preclinical.SE'=sqrt(diag(preclinical.m0.df$vcov)))
clinical.m0.df<-(predict(clinical.m0,newmods=diag(1,length(clinical.treatments))))
clinical.m0.df<-data.frame('Treatments'=clinical.treatments,'Clinical.Est'=clinical.m0.df$pred,'Clinical.SE'=clinical.m0.df$se,'Clinical.LB'=clinical.m0.df$ci.lb,
                        'Clinical.UB'=clinical.m0.df$ci.ub)
# Join prediction dfs according to treatments 
trans.m0.df<-join(clinical.m0.df,preclinical.m0.df,by='Treatments',type='inner')
trans.m0.df[,c(2,4:8)]<-exp(trans.m0.df[,c(2,4:8)])

## Fitting model
trans.m1<-glm(Clinical.Est~log(Preclinical.Est),gaussian(link='log'),start=c(0,1),weights=1/Preclinical.SE^2,data=trans.m0.df)
summary(trans.m1)

# Create prediction dataset
x.hat<-seq(0.1,3.5,.05)
y.hat<-predict(trans.m1,newdata=data.frame('Preclinical.Est'=c(x.hat)),type='response',se.fit=T)
predtrans.df<-data.frame('Treatments'=NA,'Preclinical.Est'=x.hat,'Preclinical.SE'=NA,'Preclinical.LB'=NA,'Preclinical.UB'=NA,
                              'Model'=NA,'Clinical.Est'=NA,'Clinical.SE'=NA,'Clinical.LB'=NA,'Clinical.UB'=NA,'x.hat'=x.hat,'y.hat'=NA)
predtrans.df$y.hat<-y.hat$fit
predtrans.df$y.hat.LB<-y.hat$fit-y.hat$se.fit
predtrans.df$y.hat.UB<-y.hat$fit+y.hat$se.fit

# Plotting results
# Plot preclinical MR vs clinical OS HR
(trans.plot<-ggplot(trans.m0.df,aes(as.numeric(Preclinical.Est),as.numeric(Clinical.Est),color=Treatments))+ 
  geom_point()+
  geom_errorbar(aes(ymin=(Clinical.LB), ymax=(Clinical.UB))) +
  geom_vline(xintercept = 1)+
  geom_hline(yintercept=1)+
  geom_abline(slope=1,intercept = 0,linetype=2,color='black')+
  geom_line(data=predtrans.df,aes(x=(x.hat),y=y.hat))+
  geom_ribbon(data=predtrans.df,aes(ymin=y.hat.LB,ymax=y.hat.UB),fill='blue',alpha=.1)+
  scale_x_continuous(limits=c(0,1.5))+
  scale_y_continuous(limits=c(0,2.0))+
  labs(color='Treatment',y='Clinical treatment effects (95% CI)',x="Preclinical treatment effects",
       title='Translatability potential of preclinical efficacy into clinical improvements in survival',
       subtitle=paste('Preclinical (', pre.var, ') vs clinical (',clin.var,') treatment effects',sep=''))+
  theme_minimal()+
  scale_color_discrete()+
  theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
        axis.title=element_text(size=12,face="bold",family="Helvetica")))

ggsave(trans.plot, file=paste(outputFilePath,'Preclinical_', pre.var, '_vs_Clinical_', clin.var,'.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')

# Average overestimation
trans.m0.df$Overestimation<-with(trans.m0.df,(Preclinical.Est-Clinical.Est)/Clinical.Est*100)
write.xlsx(trans.m0.df,file=paste(outputFilePath,'Clinical_',clin.var,'_Preclinical_',pre.var,'_TransEfficacyEst.xlsx',sep=''))

with(trans.m0.df,mean((Preclinical.Est-Clinical.Est)/Clinical.Est))
# For MR to PFS_HR: 14.85%
# For HR to PFS_HR: 66.58%
# For MR to OS_HR: 15.16%
# For HR to OS_HR: 67.14%

# Overall HRs were really bad at predicting individual studies HR estimates
# MRs were better at approximating the overall treatment-specific PFS_HRs, but no significant association was found
# MRs were better at predicting the individual studies of OS_HR than all other models,
# however the overall treatment effect in the clinic was overestimated
with(clinical.mat,plot(exp(HR),exp(PFS_HR)))
abline(0,1)

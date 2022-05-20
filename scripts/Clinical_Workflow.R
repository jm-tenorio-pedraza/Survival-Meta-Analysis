(rm(list=ls()))
# Clinical workflow
library(ggplot2)
library(metafor)
library(plyr)
library(openxlsx)
meta.df<-read.csv('~/Documents/GitHub/Survival-Meta-Analysis/output/clinical.red.df.csv',sep=',')
meta.df2<-read.csv("~/Documents/GitHub/Survival-Meta-Analysis/output/clinical.df.csv")
source('~/Documents/GitHub/Survival-Meta-Analysis/scripts/DataProcessingFunctions.R')
outputFilePath<-'~/Documents/Thesis/Results/'
# Unique treatment
treatments<-unique(meta.df$Treatments)
# Change the treatments to a factor
treatment.labels<-gsub('_',' + ',treatments)
treatment.labels<-gsub('anti','anti-',treatment.labels)
# Determine which dependent vars to study
var.y<-'OS_HR'
if(var.y=='OS_HR'){
  # For OS HRs:
  meta.df$y<-meta.df$OS_HR
  meta.df$sigma<-meta.df$SE_coef
  meta.df2$y<-meta.df2$OS_HR
  meta.df2$sigma<-meta.df2$SE_coef
  y.legend<-'OS hazard ratios'
  x.lab.legend<-'log(OS HRs)'
  y.lab.legend<-'Precision (1/SE)'
  } else {
    # For PFS HRss:
    meta.df$y<-exp(meta.df$PFS_HR)
    meta.df$sigma<-meta.df$PFS_SE
    meta.df<-meta.df[!is.na(meta.df$sigma),]
    
    meta.df2$y<-exp(meta.df2$PFS_HR)
    meta.df2$sigma<-meta.df2$PFS_SE
    meta.df2<-meta.df2[!is.na(meta.df2$sigma),]
    y.legend<-'PFS hazard ratios'
    x.lab.legend<-'log(PFS HRs)'
    y.lab.legend<-'Precision (1/SE)'
  }

## Result 1: Publication bias
# 1.1: Trim and fill analyses
## M0: fit mle models to treatment-based data subsets
m0.sub<-list()
m0.form<-list()
for(i in 1:length(treatments)){
  sub.df<-meta.df[is.element(meta.df$Treatments,treatments[i]),]
  m0.i<-rma(log(y), sei=sigma, slab=STUDY_ID,
            data=sub.df)
  m0.sub[[i]]<-m0.i
  m0.form[[i]]<-paste(treatments[i],'1',sep=' ~ ')
}
# Trim and fill analysis of M0 models,treatmented-based submodels
# Run trim-and-fill models for each subset of data
m0.sub.tf<-list()
for(i in 1:length(treatments)){
  try(tf.mi<-trimfill(m0.sub[[i]],'right'))
  m0.sub.tf[[i]]<-tf.mi
}
# Extract the omitted points, adjusted and original estimates for each treatment
m0.df<-getTFdf(m0.sub[[1]],m0.sub.tf[[1]])
m0.df$Treatment<-treatments[1]
for(i in 2:length(m0.sub.tf)){
  df<-getTFdf(m0.sub[[i]],m0.sub.tf[[i]])
  df$Treatment<-treatments[i]
  m0.df<-rbind(m0.df,df)
}
# Funnel plots for each treatment
m0.funnel<-plotTF(m0.df)
m0.funnel<-m0.funnel+labs(title=paste('Trim-and-fill analysis for ',y.legend,sep=''),
                          subtitle=('ICB as monotherapy and combination therapy'))+
  xlab(x.lab.legend) + ylab(y.lab.legend)
m0.funnel
# Save results
ggsave(m0.funnel, file=paste(outputFilePath,'Clinical_TrimFill_m0_', var.y,'.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')
# Number of omitted studies
k0<-sapply(m0.sub.tf,function(x) x$k0) # % 5 omitted studies in OSHR, 3 from antiCTLA4+antiPD1 and 2 from antiPD1; 1 omitted study in PFSHR from antiPDL1+Chemo

# 1.2: Adjusted estimates: effect of publication bias on estimates of treatment efficacy
estimate.unadj<-sapply(m0.sub,function(x)x$beta)
estimate.adj<-sapply(m0.sub.tf,function(x)x$beta)
ci.lb.adjusted<-sapply(m0.sub.tf,function(x)x$ci.lb)
ci.ub.adjusted<-sapply(m0.sub.tf,function(x)x$ci.ub)
N.studies<-sapply(m0.sub.tf,function(x)x$k)
est.adj.df<-data.frame('Original Estimate'=estimate.unadj,'Adjusted Estimate'=estimate.adj,'LB'=ci.lb.adjusted,'UB'=ci.ub.adjusted,
                       'Treatment' = treatments,'Model'='Treatment subset ~ 1', 'N'=N.studies,'k0'=k0)
est.adj.df
# Average overestimation of treatment effect due to publication bias
est.adj.df$Diff.Perc<-(exp(est.adj.df$Original.Estimate)-exp(est.adj.df$Adjusted.Estimate))/exp(est.adj.df$Original.Estimate)
dfIndx<-sort(est.adj.df$Diff.Perc,decreasing=F,index.return=T)$ix
est.adj.df[dfIndx,]
mean(est.adj.df$Diff.Perc)
write.xlsx(est.adj.df,file=paste(outputFilePath,'Clinical_',var.y,'_PublicationBiasAdjusted_EfficacyEst.xlsx',sep=''))

# Result 2 Heterogeneity
# 2.1: identification of heterogeneity-inducing experimental variables
# Proposed experimental vars:
model.terms<- c('CANCER','mAb','MASK')
h0<-list()
h0.m.form<-list()
for(i in 1:length(model.terms)){
  formula.mi<-as.formula(paste('~ Treatments +',model.terms[i],sep=''))
  mi<-rma(log(y),sei=sigma,mods=formula.mi,data=meta.df)
  h0[[i]]<-mi
  h0.m.form[[i]]<-paste("~Treatments +",model.terms[i],sep='')
}
h0[[i+1]]<-rma(log(y),sei=sigma,mods=~Treatments,data=meta.df)
h0.m.form[[i+1]]<-'~Treatments'
h0.Models.summ<-modelComparison(h0,h0.m.form)
h0.Models.summ

write.xlsx(h0.Models.summ,file=paste(outputFilePath,'Clinical_', var.y,'.HetSum.xlsx',sep=''))

# Including CANCER and masking:
h0.1<-rma(log(y),sei=sigma,mods=~Treatments+CANCER+MASK,data=meta.df)
h0.1

h0.unadjusted<-rma(log(y),sei=sigma,mods=~Treatments-1,data=meta.df)
h0.adjusted<-rma(log(y),sei=sigma,mods=~Treatments+CANCER+MASK-1,data=meta.df)
# 2.2: Adjusted estimates: effect of heterogeneity on estimates of treatment efficacy
estimate.unadj<-h0.unadjusted$beta
estimate.adj<-h0.adjusted$beta[1:length(treatments)]
ci.lb.adjusted<-h0.adjusted$ci.lb[1:length(treatments)]
ci.ub.adjusted<-h0.adjusted$ci.ub[1:length(treatments)]
N.studies<-h0.adjusted$k.all
est.adj.df<-data.frame('Original Estimate'=estimate.unadj,'Adjusted Estimate'=estimate.adj,'LB'=ci.lb.adjusted,'UB'=ci.ub.adjusted,
                       'Treatment' = treatments,'Model'='~Treatments +CANCER +MASK -1', 'N'=N.studies)
est.adj.df
# Average overestimation of treatment effect due to publication bias
est.adj.df$Diff.Perc<-(exp(est.adj.df$Original.Estimate)-exp(est.adj.df$Adjusted.Estimate))/exp(est.adj.df$Original.Estimate)
dfIndx<-sort(est.adj.df$Diff.Perc,decreasing=F,index.return=T)$ix
est.adj.df[dfIndx,]
mean(est.adj.df$Diff.Perc)
write.xlsx(est.adj.df,file=paste(outputFilePath,'Clinical_',var.y,'_HeterogeneityAdjusted_EfficacyEst.xlsx',sep=''))

## Result 3.3: Running trim-and-fill analyses again accounting for heterogeneity
h0.sub<-list()
listIndx<-1
factor.var<-'CANCER' # Make sure it is a factor variable in the data frame
meta.df[,factor.var]<-factor(meta.df[,factor.var])
h0.form<-data.frame(t(c(NA,NA)))
names(h0.form)<-c(factor.var,'Treatment')
# Fit model for each subset of data partitioned according to treatment
for(i in 1:length(treatments)){
  sub.df<-subset(meta.df,(meta.df$Treatments %in% treatments[i]))
  cellLines<-summary(sub.df[,factor.var])>2
  cellLines<-names(cellLines)[cellLines]
  sub.df<-subset(sub.df,sub.df[,factor.var] %in% cellLines)
  sub.df[,factor.var]<-factor(sub.df[,factor.var])
  sub.df[,factor.var]<-droplevels(sub.df[,factor.var])
  if(dim(sub.df)[1]>0)
    for(j in 1:length(cellLines)){
      sub.df.j<-subset(sub.df,sub.df[,factor.var] %in% cellLines[j])
      m0.i<-rma(log(y), sei=sigma, slab=STUDY_ID,
                data=sub.df.j)
      h0.sub[[listIndx]]<-m0.i
      h0.form[listIndx,factor.var]<-cellLines[j]
      h0.form[listIndx,'Treatment']<-treatments[i]
      listIndx<-listIndx+1
    }
  
}
h0.sub.tf<-list()
# Run trim-fill in each submodel
for(i in 1:length(h0.sub)){
  try(tf.mi<-trimfill(h0.sub[[i]],'right'))
  h0.sub.tf[[i]]<-tf.mi
}
h0.df<-getTFdf(h0.sub[[1]],h0.sub.tf[[1]])
h0.df$Treatment<-h0.form[1,'Treatment']
h0.df[,factor.var]<-h0.form[1,factor.var]
# Obtain df with omitted points and adjusted estiamtes for each treatment and concatanate them into 1 df
for(i in 2:length(h0.sub.tf)){
  df<-getTFdf(h0.sub[[i]],h0.sub.tf[[i]])
  df$Treatment<-h0.form[i,'Treatment']
  df[,factor.var]<-h0.form[i,factor.var]
  h0.df<-rbind(h0.df,df)
}
h0.df$Treatment<-factor(h0.df$Treatment,levels=treatments,labels=treatment.labels,ordered = T)
h0.df[,factor.var]<-factor(h0.df[,factor.var])
h0.df$Treatment<-droplevels(h0.df$Treatment)
# Plot and save funnel plots for each treatment
h0.funnel.list<-list()
for(i in (unique(h0.df$Treatment))){
  h0.subset<-subset(h0.df,h0.df$Treatment %in% i)
  plot.new()
  (h0.funnel.list[[i]]<-plotTFbyCell(h0.subset,factor.var))
  h0.funnel.list[[i]]<-h0.funnel.list[[i]]+title(main=paste('Trim-and-fill analysis for',y.legend,sep=' ')) +
    xlab(x.lab.legend) + ylab(y.lab.legend)
  ggsave(h0.funnel.list[[i]], file=paste(outputFilePath,'Clinical_TrimFill_', factor.var,'_', var.y, '_', i, '.png',sep=''),
         width = 10, height=8, dpi=300,bg='white')
}

# Result 4 - Model selection
## Model with treatments only
clinical.m1<-rma(log(y),sei=sigma,mods=~Treatments -1 ,data=meta.df)
clinical.p1<-data.frame(predict(clinical.m1,newmods=diag(x=1,nrow=length(clinical.m1$beta))))
clinical.p1$Treatments<-sub('Treatments','',(row.names(clinical.m1$beta)))
clinical.p1
sum(clinical.p1$ci.ub<0) # 5 treatments considered significant in OS HR; 3 significant treats in PFS HR
# Visualizing estimates
ggplot(clinical.p1,aes(y=Treatments,x=exp(pred),xmin=exp(ci.lb),xmax=exp(ci.ub)))+
  geom_point()+
  geom_errorbarh(height=.1)+
  labs(title=("Treatment estimates"),subtitle = "Model with cancer as the modifying variable")+
  xlab(x.lab.legend) +
  ylab('')+
  theme_minimal()+
  theme(text=element_text(family="Helvetica",size=8, color="black"),
        title=element_text(family='Helvetica',size=12,color='black',face='bold'))+
  theme(panel.spacing = unit(1, "lines"))+
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5)

# Model 2: Treatments + Cancer
clinical.m2<-rma(log(y),sei=sigma,mods=~Treatments+CANCER-1,data=meta.df,slab=paste(meta.df$STUDY_ID, meta.df$Treatments,sep='_'))
# identify the treatments and cell effects indexes
cancerIndx<-grep('CANCER',row.names(clinical.m2$beta))
treatIndx<-grep('Treatments',row.names(clinical.m2$beta))
# Generate prediciton matrix
newmods.mat<-(diag(1,nrow=length(treatments)))
cancerMods.mat<-matrix(rep(0,length(cancerIndx)*length(treatments)),nrow=length(treatments))
newmods.mat<-cbind(newmods.mat,cancerMods.mat)
# Generate predictions
clinical.p2<-data.frame(predict(clinical.m2,newmods=newmods.mat))
clinical.p2$Treatments<-sub('Treatments','',(row.names(clinical.m2$beta)[treatIndx]))
clinical.p2
sum(clinical.p2$ci.ub<0) # Only 4 treatments declared significant when accounting for heterogeneity and cancer-effects in OS_HR; 0 in PFS_HR

# Model 3: Treatments + Mask
clinical.m3<-rma(log(y),sei=sigma,mods=~Treatments+MASK-1,data=meta.df)
newmods.mat<-diag(1,nrow=length(treatments))
maskIndx<-grep('MASK',row.names(clinical.m3$beta))
maskMods.mat<-matrix(rep(0,length(maskIndx)*length(treatments)),nrow=length(treatments))
newmods.mat<-cbind(newmods.mat,maskMods.mat)
clinical.p3<-data.frame(predict(clinical.m3,newmods=newmods.mat))
clinical.p3$Treatments<-sub('Treatments','',(row.names(clinical.m3$beta)[treatIndx]))
clinical.p3
sum(clinical.p3$ci.ub<0) # 7 are significant in OS_HRs; 7 are significant in PFS_HRs

# Model 4: Treatments + CANCER + Mask
clinical.m4<-rma(log(y),sei=sigma,mods=~Treatments+CANCER+MASK-1,data=meta.df)
newmods.mat<-diag(1,nrow=length(treatments))
maskIndx<-grep('MASK',row.names(clinical.m4$beta))
cancerIndx<-grep('CANCER',row.names(clinical.m4$beta))
maskMods.mat<-matrix(rep(0,length(maskIndx)*length(treatments)),nrow=length(treatments))
cancerMods.mat<-matrix(rep(0,length(cancerIndx)*length(treatments)),nrow=length(treatments))

newmods.mat<-cbind(newmods.mat,cancerMods.mat,maskMods.mat)
clinical.p4<-data.frame(predict(clinical.m4,newmods=newmods.mat))
clinical.p4$Treatments<-sub('Treatments','',(row.names(clinical.m4$beta)[treatIndx]))
clinical.p4
sum(clinical.p4$ci.ub<0) #6 are significant in OS_HRs ; 4 are significant in PFS_HRs
# 4.5 Compare and choose model
AICs<-c(AIC(clinical.m1),AIC(clinical.m2),AIC(clinical.m3),AIC(clinical.m4))
I2<-c(clinical.m1$I2,clinical.m2$I2,clinical.m3$I2,clinical.m4$I2)
cbind(AICs,I2)
clinical.m0<-clinical.m4 # For OS_HR, m1 minimizes AIC, while m4 minimizes I2; For PFS_HR, m4 minimizes AIC and I^2
# However for both variables, the null model is chosen to compare against preclinical estiamtes
# Obtain estimates for chosen model
estimate.adj<-clinical.m0$beta[1:length(treatments)]
ci.lb.adjusted<-clinical.m0$ci.lb[1:length(treatments)]
ci.ub.adjusted<-clinical.m0$ci.ub[1:length(treatments)]
N.studies<-clinical.m0$k
est.adj.df<-data.frame('Original Estimate'=exp(clinical.m1$beta[1:length(treatments)]),
                       'Original.LB'=exp(clinical.m1$ci.lb[1:length(treatments)]),
                       'Original.UB'=exp(clinical.m1$ci.ub[1:length(treatments)]),
                       'Adjusted Estimate'=exp(estimate.adj),
                       'LB'=exp(ci.lb.adjusted),
                       'UB'=exp(ci.ub.adjusted),
                       'Treatment' = treatments,
                       'Model'=paste(as.character(clinical.m0$formula.mods),collapse=''),
                       'N'=N.studies)
est.adj.df
write.xlsx(est.adj.df,file=paste(outputFilePath,'Clinical_',var.y,'_MASK_adjusted_Estimates.xlsx'))

## Result 5: Fit model to all data for comparison to preclinical estimates
clinical.m0<-clinical.m1
clinical.m1<-clinical.m2

# Leave-one-out CV
cancerType<-summary(factor(meta.df$CANCER))
cancerType<-names(cancerType[cancerType>1])
meta.df<-meta.df[meta.df$CANCER%in%cancerType,]
meta.df<-meta.df[!is.na(meta.df$sigma),]
k<-dim(meta.df)[1]
pred.df<-data.frame('Study'=NA,'Treatments'=NA,'Cancer'=NA,y=NA,'Prediction'=NA)
for(i in 1:k){
  sub.df<-meta.df[-i,]
  model.i<-rma(log(y),sei=sigma,mods=~Treatments+CANCER-1,data=sub.df)
  newmods.i<-(t(sub('Treatments','',row.names(model.i$beta))%in%meta.df[i,'Treatments'])+t(sub('CANCER','',row.names(model.i$beta))%in%meta.df[i,'CANCER']))
  pred.i<-predict(model.i,newmods=newmods.i)
  pred.df[i,]<-c('Study'=meta.df[i,"STUDY_ID"],'Treatments'=meta.df[i,'Treatments'],'Cancer'=meta.df[i,'CANCER'],y=meta.df[i,'y'],'Prediction'=pred.i$pred)
}

## Save chosen model
save(clinical.m0,file=paste(outputFilePath,var.y,'_','clinical.m0.RData',sep=''))
save(clinical.m1,file=paste(outputFilePath,var.y,'_','clinical.m1.RData',sep=''))
write.table(pred.df,file=paste(outputFilePath,var.y,'_','Clinicical_LOO_Predictions.m1.csv',sep=''))


rm(list=ls())
load('~/Documents/GitHub/Survival-Meta-Analysis/output/PreclinicalSubset.RData')
source('~/Documents/GitHub/Survival-Meta-Analysis/functions/DataProcessingFunctions.R')
library(ggplot2)
library(openxlsx)
library(metafor)
outputFilePath<-'~/Documents/Thesis/Results/'
plotTreatments<-c('antiCTLA4','antiPD1','antiPDL1','antiCTLA4_antiPD1','antiCTLA4_antiPDL1',
                  'antiCTLA4_Chemotherapy','antiPD1_Chemotherapy','antiPDL1_Chemotherapy')
## Workflow for analysis of HRs and MRs
meta.df<-HRsubset.mAb
# Unique treatments
treatments<-unique(meta.df$Treatments)
# Change the treatments to a factor
treatment.labels<-gsub('_',' + ',treatments)
treatment.labels<-gsub('anti','anti-',treatment.labels)
# Determine which dependent vars to study
var.y<-'HR'
if(var.y=='HR'){
  # For HRs:
  meta.df$y<-meta.df$HR
  meta.df$sigma<-meta.df$SE_coef} else {
  # For MRs:
  meta.df$y<-meta.df$MR
  meta.df$sigma<-sqrt(1/meta.df$N)
}
# Determine title, x and y axes labels based on studied dependent variable
if(var.y=='HR'){
  y.legend<-'log-Hazard ratios'
  x.lab.legend<-'log(HR)'
  y.lab.legend<-'Precision (1/SE)'
}else
{ y.legend<-'log-Median survival ratios'
x.lab.legend<-'log(MR)'
y.lab.legend<-'Precision (N)'}
## Result 1: Effect of randomisation, blinded assessment of outcome and model replication
# Effect of randomization on estimates: Not significant
m0.rand<-rma(log(y),sei=sigma,mods=~RANDOM_ALLOCATION,data=meta.df)
m0.rand
# Effect of blinded assessment: Not significant
m0.blind<-rma(log(y),sei=sigma,mods=~BLINDED_ASSESMENT,data=meta.df)
m0.blind
# Effect of model-replication: Not significant
m0.modelrep<-rma(log(y),sei=sigma,mods=~MODEL_REPLICATION,data=meta.df)
m0.modelrep

## Result 2: Publication bias
# 2.1: Trim and fill analyses
## M0: fit mle models to treatment-based data subsets
m0.sub<-list()
m0.form<-list()
for(i in 1:length(treatments)){
  sub.df<-meta.df[is.element(meta.df$Treatments,treatments[i]),]
  m0.i<-rma(log(y), sei=sigma, slab=Study,
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
m0.df$Fill.in<-factor(m0.df$Fill.in)
# Funnel plots for each treatment
m0.funnel<-plotTF(m0.df)
m0.funnel<-m0.funnel+labs(title=paste('Trim-and-fill analysis for ',y.legend,sep=''),
                          subtitle=('ICB as monotherapy and combination therapy'))+
  xlab(x.lab.legend) + ylab(y.lab.legend)
m0.funnel
# Save results
ggsave(m0.funnel, file=paste(outputFilePath,'Preclinical_TrimFill_Full_m0_', var.y,'.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')
# Plot Subset of treatments
m0.df.sub<-subset(m0.df,m0.df$Treatment%in%plotTreatments)
# Add 1 obs that is fill-in to correct color coding
m0.df.sub<-rbind(m0.df.sub[1,],m0.df.sub)
m0.df.sub[1,'Fill.in']<-'Fill-in data'
m0.df.sub[1,c('points','se')]<-c(0,0)
m0.funnel.sub<-plotTF(m0.df.sub)
m0.funnel.sub<-m0.funnel.sub + labs(title=paste('Trim-and-fill analysis for ',y.legend,sep=''),
                     subtitle=('ICB as monotherapy and combination therapy'))+
  xlab(x.lab.legend) + ylab(
    y.lab.legend)
m0.funnel.sub 
ggsave(m0.funnel.sub, file=paste(outputFilePath,'Preclinical_TrimFill_Subset_m0_', var.y,'.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')
# Number of omitted studies
k0<-sapply(m0.sub.tf,function(x) x$k0)

# 2.2: Adjusted estimates: effect of publication bias on estimates of treatment efficacy
estimate.unadj<-sapply(m0.sub,function(x)x$beta)
estimate.adj<-sapply(m0.sub.tf,function(x)x$beta)
ci.lb.adjusted<-sapply(m0.sub.tf,function(x)x$ci.lb)
ci.ub.adjusted<-sapply(m0.sub.tf,function(x)x$ci.ub)
N.studies<-sapply(m0.sub.tf,function(x)x$k)
est.adj.df<-data.frame('Original Estimate'=estimate.unadj,'Adjusted Estimate'=estimate.adj,'LB'=ci.lb.adjusted,'UB'=ci.ub.adjusted,
                       'Treatment' = treatments,'Model'='Treatment subset ~ 1', 'N'=N.studies,'k0'=k0)
est.adj.df
write.xlsx(est.adj.df,file=paste(outputFilePath,'Preclinical_',var.y,'_PublicationBiasAdjusted_EfficacyEst.xlsx',sep=''))
# Average overestimation of treatment effect due to publication bias
est.adj.df$Diff.Perc<-(exp(est.adj.df$Original.Estimate)-exp(est.adj.df$Adjusted.Estimate))/exp(est.adj.df$Original.Estimate)
dfIndx<-sort(est.adj.df$Diff.Perc,decreasing=F,index.return=T)$ix
est.adj.df[dfIndx,]
mean(est.adj.df$Diff.Perc)

# Relation between number of ommited studies and magnitude of difference between adj. and original estimates
with(est.adj.df,plot(k0,-Diff.Perc))
abline(a=.2,b=0,h=.2)

# Relation between number of studies and magnitude of difference: Larger N non-linear association with lower difference for above-average deviations
with(est.adj.df,plot(N-k0,(-Diff.Perc)))
abline(a=.2,b=0,h=.2)

# Result 3: Heterogeneity
# 3.1: identification of heterogeneity-inducing experimental variables
# Proposed experimental vars:
model.terms<- c('CELL_FAMILY','SEX_1_female','T0_cells','BASELINE_absdiff','DOSE',
                'SITE','ROUTE','RANDOM_ALLOCATION','BLINDED_ASSESMENT','INSTITUTE_ID',
                'LAB_ID','Study','CANCER_DISEASE','CANCER_TYPE','STRAIN')
h0<-list()
h0.m.form<-list()
for(i in 1:length(model.terms)){
  formula.mi<-as.formula(paste('~Treatments',model.terms[i],sep='+'))
  mi<-rma(log(y),sei=sigma,mods=formula.mi,data=meta.df)
  h0[[i]]<-mi
  h0.m.form[[i]]<-paste('~Treatments',model.terms[i],sep=' + ')
}
h0.Models.summ<-modelComparison(h0,h0.m.form)
r2.Indx<-sort(h0.Models.summ$R.2,decreasing=T,index.return=T)$ix
h0.Models.summ<-h0.Models.summ[r2.Indx,]
h0.Models.summ
write.xlsx(h0.Models.summ,file=paste(outputFilePath,'Preclinical_', var.y,'.HetSum.xlsx',sep=''))

# Result 3.2: Model with CELL_FAMILY + Study
h0.1<-rma(log(y),sei=sigma,mods=~Treatments+CELL_FAMILY+LAB_ID,data=meta.df) # HS: Marginally better at explaining heterogeneity than the single var model

h0.unadjusted<-rma(log(y),sei=sigma,mods=~Treatments-1,data=meta.df)
h0.adjusted<-rma(log(y),sei=sigma,mods=~Treatments+CELL_FAMILY+Study-1,data=meta.df)
# 2.2: Adjusted estimates: effect of heterogeneity on estimates of treatment efficacy
estimate.unadj<-h0.unadjusted$beta
estimate.adj<-h0.adjusted$beta[1:length(treatments)]
ci.lb.adjusted<-h0.adjusted$ci.lb[1:length(treatments)]
ci.ub.adjusted<-h0.adjusted$ci.ub[1:length(treatments)]
N.studies<-h0.adjusted$k.all
est.adj.df<-data.frame('Original Estimate'=estimate.unadj,'Adjusted Estimate'=estimate.adj,'LB'=ci.lb.adjusted,'UB'=ci.ub.adjusted,
                       'Treatment' = treatments,'Model'='~Treatments +CELL_FAMILY +Study -1', 'N'=N.studies)
est.adj.df
# Average overestimation of treatment effect due to publication bias
est.adj.df$Diff.Perc<-(exp(est.adj.df$Original.Estimate)-exp(est.adj.df$Adjusted.Estimate))/exp(est.adj.df$Original.Estimate)
dfIndx<-sort(est.adj.df$Diff.Perc,decreasing=F,index.return=T)$ix
est.adj.df[dfIndx,]
mean(est.adj.df$Diff.Perc)
write.xlsx(est.adj.df,file=paste(outputFilePath,'Preclinical_',var.y,'_HeterogeneityAdjusted_EfficacyEst.xlsx',sep=''))


## Result 3.3: Running trim-and-fill analyses again accounting for heterogeneity
h0.sub<-list()
listIndx<-1
factor.var<-'Study'
h0.form<-data.frame(t(c(NA,NA)))
names(h0.form)<-c(factor.var,'Treatment')
# Fit model for each subset of data partitioned according to treatment
for(i in 1:length(treatments)){
  sub.df<-subset(meta.df,(meta.df$Treatments %in% treatments[i]))
  cellLines<-summary(sub.df[,factor.var])>2
  cellLines<-names(cellLines)[cellLines]
  sub.df<-subset(sub.df,sub.df[,factor.var] %in% cellLines)
  sub.df[,factor.var]<-droplevels(sub.df[,factor.var])
  if(dim(sub.df)[1]>0)
    for(j in 1:length(cellLines)){
      sub.df.j<-subset(sub.df,sub.df[,factor.var] %in% cellLines[j])
      m0.i<-rma(log(y), sei=sigma, slab=Study,
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
  ggsave(h0.funnel.list[[i]], file=paste(outputFilePath,'Preclinicical_TrimFill_', factor.var,'_', var.y, '_', i, '.png',sep=''),
         width = 10, height=8, dpi=300,bg='white')
}

# Result 4 - Model selection
## Model 1: with treatments only
preclinical.m1<-rma(log(y),sei=sigma,mods=~Treatments-1,data=meta.df)
# Generate prediction matrix
newmods<-diag(length(treatments))
# Compute predictions
preclinical.p1<-data.frame(predict(preclinical.m1,newmods=newmods))
preclinical.p1$Treatments<-sub('Treatments','',row.names(preclinical.m1$beta))
preclinical.p1
sum(preclinical.p1$pi.ub<0) # 8 treatment predictions considered significant in HRs, 9 in MRs

# Model 2: with treatments and cell lines
preclinical.m2<-rma(log(y),sei=sigma,mods=~Treatments+CELL_FAMILY - 1,data=meta.df)
# identify the treatments and cell effects indexes
cellIndx<-grep('CELL_FAMILY',row.names(preclinical.m2$beta))
treatIndx<-grep('Treatments',row.names(preclinical.m2$beta))

# Predictions for treatments with the cell line with the most positive effect
medianmods.mat<-matrix(rep(0,length(cellIndx)*length(treatments)),nrow=length(treatments))
newmods.mat<-cbind(newmods,medianmods.mat)
preclinical.p2<-data.frame(predict(preclinical.m2,newmods=newmods.mat))
preclinical.p2$Treatments<-sub('Treatments','',row.names(preclinical.m2$beta)[1:length(treatments)])
preclinical.p2
sum(preclinical.p2$pi.ub<0) # 33 treatments declared significant when accounting for heterogeneity and cell line-effects in HRs,6 in MRs

# Visualizing predictions
ggplot(preclinical.p2,aes(y=Treatments,x=exp(pred),xmin=exp(pi.lb),xmax=exp(pi.ub)))+
  geom_point()+
  geom_errorbarh(height=.1)+
  labs(title=("Treatment predictions"),subtitle = "Model with cell line as the modifying variable")+
  xlab(x.lab.legend) +
  ylab('')+
  theme_minimal()+
  theme(text=element_text(family="Helvetica",size=8, color="black"),
        title=element_text(family='Helvetica',size=12,color='black',face='bold'))+
  theme(panel.spacing = unit(1, "lines"))+
  geom_vline(xintercept=1, color="black", linetype="dashed", alpha=.5)

# Model 3: CELL + INSTITUTE_ID
preclinical.m3<-rma(log(y),sei=sigma,mods=~Treatments + CELL_FAMILY + INSTITUTE_ID- 1,data=meta.df)
anova(preclinical.m3,btt='INSTITUTE_ID') # Wald-type test to determine whether factor is significant

treatIndx<-grep('Treatments',row.names(preclinical.m3$beta))
cellIndx<-grep('CELL_FAMILY',row.names(preclinical.m3$beta))
instIndx<-grep('INSTITUTE_ID',row.names(preclinical.m3$beta))

# Predictions for treatments wiht the cell line with the most positive effect
cellmedianmods.mat<-matrix(rep(0,length(treatments)*length(cellIndx)),nrow=length(treatments))
instmedianmods.mat<-matrix(rep(0,length(treatments)*length(instIndx)),nrow=length(treatments))
newmods.mat<-cbind(newmods,cellmedianmods.mat,instmedianmods.mat)
preclinical.p3<-data.frame(predict(preclinical.m3,newmods=newmods.mat))
preclinical.p3$Treatments<-sub('Treatments','',row.names(preclinical.m3$beta)[1:length(treatments)])
preclinical.p3
sum(preclinical.p3$pi.ub<0) # None are considered significant, higher uncertainty and heterogeneity in HRs. None are considered significant in MRs

# Model 4: CELL_FAMILY + BASELINE_absdiff
preclinical.m4<-rma(log(y),sei=sigma,mods=~Treatments + CELL_FAMILY + BASELINE_absdiff-1,data=meta.df)
anova(preclinical.m4,btt='BASELINE_absdiff') # not significant

# Model 5: CANCER_TYPE 
preclinical.m5<-rma(log(y),sei=sigma,mods=~Treatments + CANCER_TYPE -1,data=meta.df)

# Model comparison
AICs<-c(AIC(preclinical.m1),AIC(preclinical.m2),AIC(preclinical.m3),AIC(preclinical.m4),AIC(preclinical.m5))
I2<-c(preclinical.m1$I2,preclinical.m2$I2,preclinical.m3$I2,preclinical.m4$I2,preclinical.m5$I2)
cbind(AICs,I2) #For HRs m3 minimizes AIC and I2;  For MRs, m2 minimizes AIC and has relatively low I2
preclinical.m0<-preclinical.m2
# Obtain estimates for chosen model
estimate.adj<-preclinical.m0$beta[1:length(treatments)]
ci.lb.adjusted<-preclinical.m0$ci.lb[1:length(treatments)]
ci.ub.adjusted<-preclinical.m0$ci.ub[1:length(treatments)]
N.studies<-preclinical.m0$k
est.adj.df<-data.frame('Original Estimate'=exp(preclinical.m1$beta[1:length(treatments)]),
                       'Original.LB'=exp(preclinical.m1$ci.lb[1:length(treatments)]),
                       'Original.UB'=exp(preclinical.m1$ci.ub[1:length(treatments)]),
                       'Adjusted Estimate'=exp(estimate.adj),
                       'LB'=exp(ci.lb.adjusted),
                       'UB'=exp(ci.ub.adjusted),
                       'Treatment' = treatments,
                       'Model'=paste(as.character(preclinical.m0$formula.mods),collapse=''),
                       'N'=N.studies)
est.adj.df
# Average overestimation of treatment effect due to publication bias
est.adj.df$Diff.Perc<-(exp(est.adj.df$Original.Estimate)-exp(est.adj.df$Adjusted.Estimate))/exp(est.adj.df$Original.Estimate)
dfIndx<-sort(est.adj.df$Diff.Perc,decreasing=F,index.return=T)$ix
est.adj.df[dfIndx,]
mean(est.adj.df$Diff.Perc)
write.xlsx(est.adj.df,file=paste(outputFilePath,'Preclinical_',var.y,'_CELL_INSTITUTE_adjusted_Estimates.xlsx'))
# Set m0 to the best model
preclinical.m0<-preclinical.m2 # m2 for MR and m3 for HR
## Match cell line to cancer_type
cell<-data.frame(meta.df$CELL_FAMILY,meta.df$CANCER_TYPE)
cancer.type<-sort(unique(cell$meta.df.CANCER_TYPE))
cancer.cellLines<-lapply(cancer.type,function(x)as.character(unique(droplevels(cell[cell$meta.df.CANCER_TYPE==x,1]))))
cell.family<-sub('CELL_FAMILY','',row.names(preclinical.m0$beta)[grep('CELL_FAMILY',row.names(preclinical.m0$beta))])
# Create prediction matrix with 1 col for each cell line and 1 row for each cancer type
cancer.mat<-data.frame(t(sapply(cancer.cellLines,function(x)cell.family%in%x))*1,row.names = cancer.type)
names(cancer.mat)<-cell.family
cancer.mat<-cancer.mat/rowSums(cancer.mat)
# Assign output preclinical models
preclinical.m1<-preclinical.m2
preclinical.m0<-preclinical.m1
# Save outputs
save.image(file=paste(outputFilePath,'Preclinical_',var.y,'_workflowOutput.RData',sep=''))
save(preclinical.m0,file=paste(outputFilePath,var.y,'_','preclinical.m0.RData',sep=''))
save(preclinical.m1,file=paste(outputFilePath,var.y,'_','preclinical.m1.RData',sep=''))
write.table(cancer.mat,file=paste(outputFilePath,'Preclinical_',var.y,'CELL_predMat.csv'),sep='')
save(h0.1,file=paste(outputFilePath,var.y,'_','simModel.RData',sep=''))


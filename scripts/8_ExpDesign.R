rm(list=ls())
set.seed(123)
## Simulation for experiment optimization
library(reshape2)
library(parallel)
library(metafor)
library(ggplot2)
library(openxlsx)
setwd('/') # Set working directory
outputFilePath<-paste(getwd(),'/output',sep='')
# Choose variable of interest {MR for median survival ratios or HR for hazard ratios}
var.y<-'MR'
if(var.y=='HR'){
  y.legend<-'log-Hazard ratios'
  x.lab.legend<-'log(HR)'
  y.lab.legend<-'Precision (1/SE)'
}else
{ y.legend<-'log-Median survival ratios'
x.lab.legend<-'log(MR)'
y.lab.legend<-'Precision (N)'}

# Load model to simulate
sim.model<-load(file=paste(outputFilePath,var.y,'_simModel.RData',sep=''))
# Extract heterogeneity
tau<-sim.model$tau2
# Extract variability
xi<-sim.model$vi
# Extract estimated coefficients
beta.est<-sim.model$beta
# Identify cell-related coefs
cell_indx<-grep('CELL_FAMILY.*',row.names(beta.est))
cell.est<-beta.est[cell_indx,1]
cell.psi<-sd(cell.est)
names(cell.est)<-sub('CELL_FAMILY','',names(cell.est))
n_cell<-length(cell.est)

# Identify lab-related coefs
lab_indx<-grep('Study.*',row.names(beta.est))
lab.est<-beta.est[lab_indx,1]
lab.psi<-sd(lab.est)
names(lab.est)<-sub('Study','',names(lab.est))
n_lab<-length(lab.est)

# identify treatment-related coefs
treat_indx<-grep('Treatments.*',row.names(beta.est))
treat.est<-beta.est[treat_indx,1]
names(treat.est)<-sub('Treatments','',names(treat.est))
# Select treatments of interest
treatments<-c('antiCTLA4','antiPD1','antiPDL1','antiCTLA4_antiPD1','antiCTLA4_antiPDL1','antiCTLA4_Chemotherapy',
              'antiPD1_Chemotherapy','antiPDL1_Chemotherapy')
treat.est<-treat.est[names(treat.est)%in%treatments]
n_treat<-length(treat.est)

# Number of artificial datasets
N_samples<-2000
N_cell<-6
N_lab<-10
N_replicates<-3
# Number of cores
numCores<-detectCores() #useful only if parallelization is available
# Generate samples for each exp design
MSE.exp.design<-list()
HET.exp.design<-list()
exp.indx<-1
fx<-function(x,model.formula){ # Define function that extracts relevant estimates from model
  model.i<-rma(MR,sei=SD,mods=model.formula,data=x)
  beta.ij<-model.i$beta[1:length(treat.est)]
  se.ij<-model.i$se[1:length(treat.est)]
  pval.ij<-model.i$pval[1:length(treat.est)]
  tau.ij<-model.i$tau2
  I2.ij<-model.i$I2
  output<-data.frame('Treatments'=names(treat.est),'Treatment.Est'=beta.ij,
                     'SE.Est'=se.ij,'MSE'=(beta.ij-treat.est)^2,'Power'=1-pval.ij,'Tau'=tau.ij,'I2'=I2.ij)
  return(output)
}
default.fx<-function(){ # Define function to output a df when models cannot be fitted due to low number of factor levels 
  output<-data.frame('Treatments'=names(treat.est),'Treatment.Est'=NaN,
                     'SE.Est'=NaN,'MSE'=NaN,'Power'=NaN,'Tau'=NaN,'I2'=NaN)
  return(output)
}

for(i in 1:N_cell){ # Loop to simulate data and estimate models
  for(j in 1:N_lab){
    # Sample cell-related effects
    sample_cell<-matrix(rep(matrix(sample(cell.est,N_samples*i,replace=T),nrow=N_samples),j),nrow=N_samples)
    # sample_cell<-matrix(rep(matrix(rnorm(N_samples*i,sd=cell.psi),nrow=N_samples),j),nrow=N_samples)
    # Sample lab-related effects
    sample_lab<-matrix(sample(lab.est,N_samples*j,replace=T),nrow=N_samples)
    # sample_lab<-matrix(rnorm(N_samples*j,sd=lab.psi),nrow=N_samples)
    sample_lab <- matrix(data = apply(sample_lab, 2, function(x) rep(x, i)), ncol = ncol(sample_lab)*i)
    # Sample precision estimates
    sample_var<-sample(sqrt(xi),N_samples*i*j*n_treat,replace=T)
    # Sample residual heterogeneity-effect
    het_error<-(rnorm(N_samples*n_treat*i*j*N_replicates,sd=sqrt(tau)))
    # Calculate random effects
    rand.effs<-sample_cell+sample_lab
    fixed.effs<-matrix(rep(rep(treat.est,N_samples),dim(rand.effs)[2]),nrow=N_samples*n_treat)
    # Caculate the obs
    rand.effs<-matrix(data = apply(rand.effs, 1, function(x) rep(x, each=n_treat)), ncol = dim(rand.effs)[2])
    yi<-fixed.effs+rand.effs
    # Create labels for the cell and lab combinations
    cell_labels<-rep(paste("CELL",1:i,sep='.'),j)
    lab_labels<-rep(paste('LAB',1:j,sep='.'),each=i)
    # Create labels for each study dataset
    study_labels<-rep(paste('Study',1:N_samples,sep='_'),each=n_treat)
    # Consolidate observations into dataframe
    hr.df<-data.frame(yi)
    names(hr.df)<-paste(cell_labels,lab_labels,sep='_')
    # Add study and treatment
    hr.df$Study_id<-study_labels
    hr.df$Treatment<-rep(names(treat.est),N_samples)
    # Reshape into long format with 1 column for the observations
    hr.df.long<-melt(hr.df,id=c('Study_id','Treatment'))
    # Split variable vector into cell- and lab-related effects
    cell_lab_var<-unlist(strsplit(as.character(hr.df.long$variable),'_'))
    cell_var<-cell_lab_var[seq(1,length(cell_lab_var),2)]
    lab_var<-cell_lab_var[seq(2,length(cell_lab_var),2)]
    # Add character vectors of cells and labs
    hr.df.long$CELL<-cell_var
    hr.df.long$LAB<-lab_var
    # Add the sampling SD
    hr.df.long$SD<-sample_var
    # Replicate each experiment
    hr.df.long<-hr.df.long[rep(seq_len(nrow(hr.df.long)), each = N_replicates), ]
    # Add the sampling SD to the yi's
    hr.df.long$MR<-hr.df.long$value+rnorm(N_samples*n_treat*i*j*N_replicates,sd=hr.df.long$SD)+het_error
    hr.list<-split(hr.df.long,as.factor(hr.df.long$Study_id))
    rm(list=c('rand.effs','fixed.effs','yi','hr.df','hr.df.long',
              'lab_var','cell_var','sample_var','sample_cell','sample_lab','cell_lab_var'))
    # Estimate model
    if(i==1 & j==1){model.form<-as.formula('~Treatment -1')}
    if(i==1 & j!=1){model.form<-as.formula('~Treatment +LAB -1')}
    if(i!=1 & j==1){model.form<-as.formula('~Treatment + CELL -1')}
    if(i!=1 & j!=1){model.form<-as.formula('~Treatment +LAB + CELL-1')}
    fitModel<-function(x,treat.names){
      out<-tryCatch(fx(x,model.form), 
                    error=function(cond){
                      
                      out<-data.frame('Treatments'=treat.names,'Treatment.Est'=NaN,
                                      'SE.Est'=NaN,'MSE'=NaN,'Power'=NaN)
                      return(out)
                    },
                    warning=function(cond){
                      out<-data.frame('Treatments'=treat.names,'Treatment.Est'=NaN,
                                      'SE.Est'=NaN,'MSE'=NaN,'Power'=NaN)
                      return(out)
                    },
                    finally={
                      message('Processed')
                    })
      return(out)
    }
    # rma.i.j<-mclapply(hr.list,fitModel,mc.cores = numCores) Only useful in mac
    cl<-makePSOCKcluster(8)
    setDefaultCluster(cl)
    adder<-function(x) fitModel(x,names(treat.est))
    clusterExport(NULL,c('adder','fitModel','fx','treat.est','model.form'))
    clusterEvalQ(NULL,library(metafor))
    clusterEvalQ(NULL,library(base))
    rma.i.j<-parLapply(NULL,hr.list,adder) # Only useful in windows
    stopCluster(cl)
    # Extract treatment effects
    beta.ij<-t(sapply(rma.i.j,function(x)x$Treatment.Est))
    #Extract SE of treatment effects
    se.ij<-t(sapply(rma.i.j,function(x)x$SE.Est))
    # Extract the p-values
    power.ij<-t(sapply(rma.i.j,function(x)x$Power))
    # Extract tau
    tau.ij<-t(sapply(rma.i.j,function(x)x$Tau))
    # Extract I2
    I2.ij<-t(sapply(rma.i.j,function(x)x$I2))
    # Dataframe 
    MSE.exp.design[[exp.indx]]<-data.frame('Treatments'=names(treat.est),
                                           'Treatment.Est'=colMeans(beta.ij,na.rm=T),
                                           'SE.est'=colMeans(se.ij,na.rm = T),
                                           'MSE'=colMeans((beta.ij-((treat.est)))^2,na.rm=T),
                                           'Power'=colMeans(power.ij,na.rm=T))
    HET.exp.design[[exp.indx]]<-data.frame('Tau'=tau.ij[,1],I2=I2.ij[,1])
    names(MSE.exp.design)[exp.indx]<-paste(paste('Design',exp.indx,sep='.'), '- Cell:',i,'Lab:',j,sep=' ')
    names(HET.exp.design)[exp.indx]<-paste(paste('Design',exp.indx,sep='.'), '- Cell:',i,'Lab:',j,sep=' ')
    
    exp.indx<-exp.indx+1 

  }
}

# Power
power.est<-sapply(MSE.exp.design,function(x)x$Power)
power.3d<-power.est
# Tall-version
power.est<-melt(power.est)
## MSE
MSE.est<-sapply(MSE.exp.design,function(x)x$MSE)
MSE.3d<-MSE.est
# Tall-version 
MSE.est<-melt(MSE.est)
MSE.est$Treatment<-rep(names(treat.est),dim(MSE.est)[1]/n_treat)
names(MSE.est)<-c('ID','Design','MSE','Treatment')
head(MSE.est)
cell_lab_var<-sub('.*- ','',as.character(MSE.est$Design))
cell_lab_var
cell<-as.numeric(sub(' Lab: .*','',sub('Cell: ','',cell_lab_var)))
lab<-as.numeric(sub('Lab: ','',sub('Cell: [0 1 2 3 4 5 6 7 8 9 10]*','',cell_lab_var)))
MSE.est$CELL<-(cell)
MSE.est$LAB<-lab
MSE.est$Power<-power.est$value
MSE.est$Treatment<-factor(MSE.est$Treatment,levels=treatments,labels=c('anti-CTLA-4','anti-PD-1','anti-PD-L1','anti-CTLA-4 + anti-PD-1',
                                                                   'anti-CTLA-4 + anti-PD-L1','anti-CTLA-4 + Chemotherapy',
                                                                   'anti-PD-1 + Chemotherapy', 'anti-PD-L1 + Chemotherapy'),ordered = T)

write.xlsx(MSE.est,file=paste(outputFilePath,var.y,'_MSE.est.summary.xlsx',sep = ''))
#MSE.est.xlsx<-read.xlsx(paste(outputFilePath,var.y, '_MSE.est.summary.xlsx',sep=''))
#MSE.3d<-MSE.est.xlsx

# 3D version of MSE
MSE.3d<-colMeans(MSE.3d,na.rm=T)
cell_lab_var<-sub('.*- ','',names(MSE.3d))
cell_lab_var
cell<-as.numeric(sub(' Lab: .*','',sub('Cell: ','',cell_lab_var)))
lab<-as.numeric(sub('Lab: ','',sub('Cell: [0 1 2 3 4 5 6 7 8 9 10]*','',cell_lab_var)))

MSE.3d<-data.frame(MSE.3d)
MSE.3d$CELL<-cell
MSE.3d$LAB<-lab
MSE.3d$variable <-'MSE'
names(MSE.3d)<-c('value','CELL','LAB','variable')
MSE.3d<-acast(MSE.3d,CELL~LAB)
write.xlsx(MSE.3d,file=paste(outputFilePath,var.y,'_MSE.3d.xlsx',sep=''))
#MSE.3d.xlsx<-read.xlsx('/Users/migueltenorio/Documents/Thesis/Results/MSE.3d.xlsx')

# 3D version of Power
power.3d<-colMeans(power.3d,na.rm=T)
power.3d<-data.frame(power.3d)
power.3d$CELL<-cell
power.3d$LAB<-lab
power.3d$variable<-'Power'
names(power.3d)<-c('value','CELL','LAB','variable')

power.3d<-acast(power.3d,CELL~LAB)
write.xlsx(power.3d,file=paste(outputFilePath,var.y,'_Power.3d.xlsx'))

# Power.3d.xlsx<-read.xlsx('/Users/migueltenorio/Documents/Thesis/Results/Power.3d.xlsx')

## 2-D plotof MSE
MSE.ggplot<-ggplot(MSE.est,aes(CELL,MSE,color=LAB))+
  geom_point()
MSE.ggplot
# SUbplot for each cell number
MSE.ggplot<-ggplot(MSE.est,aes(LAB,MSE,color=Treatment))+
  geom_point(size=1)+
  geom_line()+
  facet_wrap(.~CELL,nrow=1,ncol=6,labeller = label_both)

MSE.ggplot<-MSE.ggplot+labs(color='Treatment',  y='MSE',x="Number of studies",
                            title=paste('Mean squared error of treatment effects ', '(',var.y,')', sep=''),
                            subtitle='Experimental design: 1-6 cell lines + 1-10 studies')+
  theme_minimal()+
  scale_color_discrete()+
  theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
        axis.title=element_text(size=12,face="bold",family="Helvetica"))+
  scale_x_continuous(limits=c(0,11))
MSE.ggplot
ggsave(MSE.ggplot, file=paste(outputFilePath, var.y,'_MSE_CELL_LAB.png',sep=''),
       width = 10, height=8, dpi=300,bg = 'white')
## 2-D plot of Power
Power.ggplot<-ggplot(MSE.est,aes(LAB,Power,color=Treatment))+
  geom_point()+
  geom_line()+
  facet_wrap(~CELL)
Power.ggplot<-Power.ggplot+labs(color='Treatment',  y='Power',x="Number of studies",
                                title=paste('Power of the experimental designs for ICB treatments ','(',var.y,')',sep = ''),
                                subtitle='Experimental design: 1-6 cell lines + 1-10 studies')+
  theme_minimal()+
  scale_color_discrete()+
  theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
        axis.title=element_text(size=12,face="bold",family="Helvetica"))+
  scale_x_continuous(limits=c(0,11))
Power.ggplot
ggsave(Power.ggplot, file=paste(outputFilePath, var.y, '_Power_CELL_LAB.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')
MSE.est$CELL_factor<-as.factor(MSE.est$CELL)
# Plot power and MSE
MSE.Power.ggplot<-ggplot(MSE.est,aes(x=LAB,color=CELL_factor))+
  geom_line(aes(y=MSE*2),linetype=1)+
  geom_point(aes(y=MSE*2))+
  geom_line(aes(y=Power),linetype=2)+
  scale_color_discrete()+
  scale_y_continuous(
    name='Power', sec.axis = sec_axis(~.*1/2,name='MSE')) +
  facet_wrap(~Treatment) + theme_minimal()+
  theme(plot.title=element_text(size=16,face="bold",family="Helvetica"),
        axis.title=element_text(size=12,face="bold",family="Helvetica"))

MSE.Power.ggplot<-MSE.Power.ggplot+labs(
  title='Power and MSE of the experimental designs for ICB treatments',
  subtitle='Experimental design: 1-6 cell lines + 1-10 labs')

MSE.Power.ggplot
ggsave(MSE.Power.ggplot, file=paste(outputFilePath,var.y,'_MSE_Power_CELL_LAB.png',sep=''),
       width = 10, height=8, dpi=300,bg='white')


# 3-D plot of MSE
library(plotly)
yax <- list(title='Number of cell lines',
            ticketmode='array',
            tickvals=0:(dim(MSE.3d)[1]-1),
            ticktext=row.names(MSE.3d),
            titlefont=list(family='Helvetica',size=16))
xax <- list(title='Number of labs',
            ticketmode='array',
            tickvals=0:(dim(MSE.3d)[2]-1),
            ticktext=colnames(MSE.3d),
            titlefont=list(family='Helvetica',size=16))

MSE.plotly<-plot_ly(z=~MSE.3d)%>% 
  add_surface()
MSE.plotly<-MSE.plotly %>% plotly::layout(title='Average MSE over all treatments',
                                          scene=list(yaxis=yax,
                                                     xaxis=xax,
                                                     zaxis=list(title='MSE')))
MSE.plotly
# 3-D plot of Power
Power.plotly<-plot_ly(z=~power.3d)%>% 
  add_surface()
Power.plotly<-Power.plotly %>% plotly::layout(title='Average Power over all treatments',
                                              scene=list(yaxis=yax,
                                                         xaxis=xax,
                                                         zaxis=list(title='Power')))
Power.plotly
save.image(file=paste(outputFilePath,var.y,'_ExpDesign.RData',sep = ''))
load(paste(outputFilePath,var.y,'_ExpDesign.RData',sep = ''))

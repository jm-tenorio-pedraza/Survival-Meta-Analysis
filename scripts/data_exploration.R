## Exploration of data set for meta-analysis
# Load libraries
library(ggplot2)
library(plyr)
# Load dataset
survival<-read.csv('~/Documents/Data/CSV/Survival/Survival_Preclinical.csv',
                   header = TRUE, sep = ';', quote = "\"", dec = ',')
survival_clinical<-read.csv('~/Documents/Data/CSV/Survival/Survival_Clinical.csv',
                   header = TRUE, sep = ';', quote = "\"", dec = ',')

head(survival,10)

## Subset relevant treatments
therapy_indx <- survival$MAB=='antiPDL1'|survival$MAB=='antiPD1'|survival$MAB=='antiCTLA4'|
  survival$MAB=='antiCTLA4+antiPD1'|survival$MAB=='antiCTLA4+antiPDL1'|survival$MAB=='Control'|
  survival$MAB=='iIDO'|survival$MAB=='antiCTLA4+iIDO'|survival$MAB=='antiPDL1+iIDO'|survival$MAB=='antiCTLA4+antiPDL1+iIDO'

survival_subset <-subset(survival, is.finite(survival$MOS_COR)&therapy_indx, 
                         select = c(STUDY_ID,MOS_COR, CELL, STRAIN, MAB, LAB,N, EUTHANIZED,
                                    DOSE_mg.kg, DOSE_2, DOSE_3,T_0_cells))
head(survival_subset)
summary(survival_subset)
## Drop irrelevant levels from each category
survival_subset$MAB<-droplevels(survival_subset$MAB)
survival_subset$CELL<-droplevels(survival_subset$CELL)
survival_subset$LAB<-droplevels(survival_subset$LAB)
survival_subset$EUTHANIZED<-droplevels(survival_subset$EUTHANIZED)
survival_subset$DOSE_mg.kg<-droplevels(survival_subset$DOSE_mg.kg)

## Set default levels in each of the categories
survival<-within(survival, MAB<-relevel(MAB, 'Control'))
survival<-within(survival, STRAIN<-relevel(STRAIN,'C57BL/6'))
survival<-within(survival, CELL<-relevel(CELL,'B16F10'))

## Check levels
levels(survival_subset$MAB)
levels(survival_subset$LAB)
levels(survival_subset$CELL)

## Plotting MOS by cell line
survival_MAB<-aggregate(log(survival_subset[,2]),list(survival_subset$MAB),mean)
survival_MAB$N<-aggregate(survival_subset[,7],list(survival_subset$MAB),sum)$x
colnames(survival_MAB)<-c("MAB", "MOS_COR","N")
survival_CELL<-aggregate(log(survival_subset[,2]),list(survival_subset$CELL),mean)
survival_CELL$N<-aggregate(survival_subset[,7],list(survival_subset$CELL),sum)$x
colnames(survival_CELL)<-c("CELL", "MOS_COR","N")
survival_CELL
## Plotting by treatment
cell<-ggplot(survival_subset,aes(y=log(MOS_COR),x=MAB,colour=CELL,size=N)) +
  geom_point() +
  geom_point(data=survival_MAB,aes(y=MOS_COR),colour='red', size=5)+ 
  theme(axis.text.x = element_text(angle = 45))
cell+labs(color='Cell line', size='Number of mice', y='Log of Median Survival',
          title='Survival in preclinical experiments',
          subtitle='Effect of cell line')
## Plotting by cell
cell<-ggplot(survival_subset,aes(y=log(MOS_COR),x=CELL,colour=MAB,size=N)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45))

cell+labs(color='Treatment', size='Number of mice', y='Log of Median Survival',
          title='Survival in preclinical experiments',
          subtitle='Effect of cell line')
  
## Plotting by lab
lab<-ggplot(survival_subset,aes(y=log(MOS_COR),x=LAB,colour=MAB,size=N)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 45))
lab+labs(color='Treatment', size='Number of mice', y='Log of Median Survival',
           title='Survival in preclinical experiments',
           subtitle='Effect of lab')

## Plotting by dose
doses<-sort(unique(survival_subset$DOSE_mg.kg))
dose<-ggplot(survival_subset,aes(y=log(MOS_COR),x=MAB,colour=log(DOSE_mg.kg),size=N)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45))+
  scale_colour_gradient2(breaks=sort(log(doses[seq(1,length(doses),2)])), 
                         labels=as.character(doses[seq(1,length(doses),2)]))+
  geom_point(data=survival_MAB,aes(y=MOS_COR),colour='red',size=5)
dose+labs(color='Dose (mg/kg)', size='Number of mice',
          title='Survival in preclinical experiments',
          subtitle='Effect of dose',y='Log of Median Survival')

# Counting factors
count(survival_subset, c('LAB', 'STUDY_ID'))
count(survival_subset, c('STUDY_ID', 'MAB'))
count(survival_subset, c('LAB', 'CELL'))
count(survival_subset,c('EUTHANIZED','STUDY_ID'))
count(survival,c('CELL'))
# Plotting survival
time<-seq(0,200, 10)
s<-function(x) exp(-time*1/x)
survival_curves<-apply(survival_subset[2],1, s)
f<-function(x)
  survival_curves<-by(t(survival_curves),survival_subset$MAB,colMeans)
survival_curves<-t(matrix(unlist(survival_curves),ncol=length(time),byrow=FALSE))
survival_curves.df<-data.frame(time, survival_curves)


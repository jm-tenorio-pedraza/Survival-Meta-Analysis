# Define extensions of .csv files containing the survival curves from preclinical expts.

setwd("~/Documents/Data/CSV")

fileExt<-c('Allard_2013/Allard_2013_Fig3.csv',
           'Duraiswamy_2013/Duraiswamy_2013_Fig3i.csv',
           'Duraiswamy_2013/Duraiswamy_2013_Fig5i.csv',
           'Duraiswamy_2013/Duraiswamy_2013_Fig5ii.csv',
           'Duraiswamy_2013/Duraiswamy_2013_Fig5iii.csv',
           'Duraiswamy_2013/Duraiswamy_2013_Fig5iv.csv',
           'Lu_2014/Lu_2014_Fig1c.csv',
           'Lu_2014/Lu_2014_Fig5b.csv',
           'Lu_2014/Lu_2014_Fig5d.csv',
           'Yu_2010/Yu_2010_Fig6.csv')
fileExt <- (paste(getwd(), fileExt, sep='/'))
study <- c('ALLARD_2013', 'DURAISWAMY_2013','DURAISWAMY_2013', 'DURAISWAMY_2013',
           'DURAISWAMY_2013', 'DURAISWAMY_2013', 'LU_2014','LU_2014','LU_2014', 'YU_2010')
datasets <- lapply(fileExt, read.csv)
N_mice <- list(c(15, 15, 15, 15, 15,15),
          c(12, 12, 12, 12),
          c(12, 12, 12, 12, 12, 12, 12, 12),
          c(12, 12, 12, 12),
          c(12, 12, 12, 12),
          c(12, 12),
          c(10, 10, 10, 10),
          c(10, 10, 10, 10),
          c(10, 10, 10, 10, 10, 10),
          c(6, 6, 6, 6, 6, 6, 6))
finalT <- list(80,
               120,
               120,
               120,
               120,
               120,
               90,
               90,
               90,
               100)
dataSurvival <- list()
for(i in 1:length(datasets)){
  dataSurvival[[i]] = processSurvivalData(datasets[[i]], N_mice[[i]], finalT[[i]])
}
cell<-list()
cell<-lapply(dataSurvival, plotSurvival)

processSurvival.df <-lapply(dataSurvival,expandSurvivalFunction)
survival.df <- processSurvival.df[[1]]
survival.df$Study <- rep(study[1], dim(survival.df)[1])
for(i in 2:length(study)){
  subset_i <- processSurvival.df[[i]]
  subset_i$Study <- rep(study[i], dim(subset_i)[1])
  survival.df <- rbind(survival.df, subset_i)
}

survivalMedian.df <-lapply(processSurvival.df,function(x)by(x, x$Experiment, SurvMedians))

for(i in 1:length(study)){
  survivalMedian.df[[i]] <- createDF(survivalMedian.df[[i]], processSurvival.df[[i]], study[i])
}

survMedian.df <- survivalMedian.df[[1]]
for(i in 2:length(study)){
  survMedian.df <- rbind(survMedian.df, survivalMedian.df[[i]])
}

OS_PFS<-ggplot(survMedian.df,aes(y=Median,x=Cell,colour=paste(survMedian.df$Treatment_1,survMedian.df$Treatment_2,survMedian.df$Treatment_3,sep = '+')))+
  geom_point()+
  geom_errorbar(aes(ymin=LB, ymax=UB))

OS_PFS+labs(color='Treatment',  y='Time [days]',x="Treatment",
            title='Median survival in preclinical experiments',
            subtitle='Progression-free vs overall survival')
  
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


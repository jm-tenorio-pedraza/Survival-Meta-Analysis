# Define extensions of .csv files containing the survival curves from preclinical expts.

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

EXPERIMENT <- paste(Survival_Data$TREATMENT_1, 
                                  Survival_Data$TREATMENT_2, 
                                  Survival_Data$TREATMENT_3,
                                  Survival_Data$CELL,sep='_')
EXPERIMENT <- sub('__', '_', EXPERIMENT) # apply twice
EXPERIMENT <- sub('__', '_', EXPERIMENT) # apply twice

Survival_Data$EXPERIMENT <- EXPERIMENT
head(Survival_Data)
survival.subset <-match_df(Survival_Data, experiments.factor, on = c('STUDY_ID', 'FIG', 'EXPERIMENT'))
differences <- data.frame(setdiff(experiments.factor$EXPERIMENT, survival.subset$EXPERIMENT))
setdiff(experiments.factor$STUDY_ID, survival.subset$STUDY_ID)
setdiff(survival.subset$EXPERIMENT, experiments.factor$EXPERIMENT)
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


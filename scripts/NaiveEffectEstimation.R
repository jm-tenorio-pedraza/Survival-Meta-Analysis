## Naive estimation of HR of effect
beta_survival<-aggregate(log(survival_subset$MOS_COR),list(survival_subset$MAB), mean)
# a-CTLA4
exp(beta_survival[1,2])/(exp(beta_survival[2,2]))
# a-CTLA4+antiPD1
exp(beta_survival[1,2])/(exp(beta_survival[3,2]))
# a-CTLA4+antiPDL1
exp(beta_survival[1,2])/(exp(beta_survival[4,2]))
# a-CTLA4+antiPDL1+iIDO
exp(beta_survival[1,2])/(exp(beta_survival[5,2]))
# a-PD1
exp(beta_survival[1,2])/(exp(beta_survival[6,2]))
# a-PDL1
exp(beta_survival[1,2])/(exp(beta_survival[7,2]))
# i-IDO
exp(beta_survival[1,2])/(exp(beta_survival[8,2]))

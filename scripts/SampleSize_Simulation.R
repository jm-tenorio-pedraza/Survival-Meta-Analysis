# Simulation study to calculate sample sizes

library(MCMCglmm)
# Previous results

MAB <- as.matrix(model3_glm_red$Sol[,2:10])
colnames(MAB)<- colnames(model3_glm_red$Sol[,2:10])
b0 = as.matrix(model3_glm_red$Sol[,1])
colnames(b0)<- 'Control'
V_LAB <- as.matrix(model3_glm_red$VCV[,1])
V_CELL <- as.matrix(model3_glm_red$VCV[,2])
Verror <- as.matrix(model3_glm_red$VCV[,3])

 # Asssuming a complete cross-design
sim1 <- function(MAB=0, b0=1, V_LAB=1, V_CELL=1, Verror=1,n_mab, n_labs, n_cell) { # Inputs to this function are the posterior samples
  LAB <- rep( 1:n_labs ,each=n_cell*(n_mab+1))
  CELL <- rep( 1:n_cell, n_labs*(n_mab+1))
  treatment <- rep(rep(0:n_mab, each=n_cell),n_labs)
  
  # Randomly choose an effect size and an intercept
  b_MAB = MAB[sample(1:nrow(MAB), 1, replace = TRUE),]
  b_0 <- b0[sample(1:nrow(b0),1),]
  # Subset V_LAB and V_CELL to be the size of the sample
  V_LAB_hat = V_LAB[sample(1:nrow(V_LAB), n_labs, replace = TRUE),]
  V_CELL_hat = V_CELL[sample(1:nrow(V_CELL),n_cell, replace = TRUE),]
  
  # Randomly choose an error variance term
  V_error_hat = Verror[sample(1:nrow(Verror), n_labs*n_cell*n_mab, replace=TRUE),]
  
  # random effects per lab
  S.re <- rnorm(n_labs, 0, sqrt(V_LAB_hat)) # Check to see if this produces 1 observation per lab sampled
  
  # random effects per cell
  W.re <- rnorm(n_cell, 0, sqrt(V_CELL))
  
  # epsilons
  eps <- rnorm(n_labs*n_cell*(n_mab+1),0, sqrt(Verror))
  b = rep(0,n_labs*n_cell*(n_mab+1));
  b[treatment!=0] = b_MAB[treatment]
  
  # put it all together
  MOS <- b_0 + b + S.re[LAB] + W.re[CELL] + eps
  
  MAB_char <- rep(0, n_labs*n_cell*(n_mab+1))
  MAB_char[treatment==0] <- colnames(b0)
  MAB_char[treatment!=0] <- colnames(MAB)[treatment]
  # put into a data frame
  mydata <- data.frame( LAB  = paste('LAB',LAB, sep='_'), 
                        CELL = paste('CELL', CELL, sep='_'), 
                        MAB  = factor(MAB_char),
                        MOS = MOS)
  mydata<-within(mydata, MAB<-relevel(MAB, 'Control'))
  
  # analyze looking at interaction term with LR test
  prior3 = list(R = list(V = 1, nu = 2), G = list(G1 = list(V = 1,
                                                                    nu =2), G2 = list(V = 1, nu = 2)))
  model3_glm_red <- MCMCglmm((MOS)~MAB, random=~CELL+LAB, data=mydata,
                         prior = prior3,pr=FALSE,nitt=1e4,burnin=1e3, thin=10)
  summary(model3_glm_red)[5]$solutions[,5] # p-values for the two-tailed tests
}
out <- sim1(MAB, b0, V_LAB, V_CELL, Verror, n_mab = 9, n_labs = 3, n_cell=2)
pb <- txtProgressBar(max=100) # or tkProgressBar or txtProgressbar

setTxtProgressBar(pb, 0)
PowerSampleSizes <- matrix(rep(0,10*6),nrow=10)
PowerSampleSD <- matrix(rep(0,10*6),nrow=10)

out1 <- replicate( 100, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
  sim1( MAB, b0, V_LAB, V_CELL, Verror, n_mab= 9, n_labs = 1, n_cell = 2)})
PowerSampleSizes[,1]<-rowMeans(out1)
PowerSampleSD[,1]<-apply(1-out1, 1, 'sd')

out2 <- replicate( 100, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
  sim1( MAB, b0, V_LAB, V_CELL, Verror, n_mab= 9, n_labs = 2, n_cell = 2)})
PowerSampleSizes[,2]<-rowMeans(out2)
PowerSampleSD[,2]<-apply(1-out2, 1, 'sd')

out4 <- replicate( 100, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
  sim1( MAB, b0, V_LAB, V_CELL, Verror, n_mab= 9, n_labs = 4, n_cell = 2)})
PowerSampleSizes[,3]<-rowMeans(out4)
PowerSampleSD[,3]<-apply(1-out4, 1, 'sd')

out6 <- replicate( 100, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
  sim1( MAB, b0, V_LAB, V_CELL, Verror, n_mab= 9, n_labs = 6, n_cell = 2)})
PowerSampleSizes[,4]<-rowMeans(out6)
PowerSampleSD[,4]<-apply(1-out6, 1, 'sd')

out8 <- replicate( 100, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
  sim1( MAB, b0, V_LAB, V_CELL, Verror, n_mab= 9, n_labs = 8, n_cell = 2)})
PowerSampleSizes[,5]<-rowMeans(out8)
PowerSampleSD[,5]<-apply(1-out8, 1, 'sd')

out16 <- replicate( 100, {setTxtProgressBar(pb, getTxtProgressBar(pb)+1);
  sim1( MAB, b0, V_LAB, V_CELL, Verror, n_mab= 9, n_labs = 16, n_cell = 2)})
PowerSampleSizes[,6]<-rowMeans(out16)
PowerSampleSD[,6]<-apply(1-out16, 1, 'sd')

# Plotting results
rownames(PowerSampleSizes) <- levels(survival_subset$MAB)
nlab_sampleSizes<- c(1,2,4,6,8,16)
power<- data.frame(t(rbind(nlab_sampleSizes, 1-PowerSampleSizes)))
power_sd <- data.frame(t(rbind(nlab_sampleSizes, PowerSampleSD)))
library(ggplot2)
library(reshape2)
power <- melt(power, id.vars="nlab_sampleSizes")
power_sd <- melt(power_sd, id.vars = 'nlab_sampleSizes')
power$sd <- power_sd$value
# Everything on the same plot
nlab_power <- ggplot(power, aes(nlab_sampleSizes,value, col=variable)) + 
  geom_point() +
  geom_line()
nlab_power+labs(color='Treatment',  y='Power',x="Number of sampled labs",
            title='Lab random effect on power to detect treatment effect',
            subtitle='Mean power of n=100 artificial datasets')

nlab_power +  geom_errorbar(aes(ymin= value-sd, ymax= value+sd))+
  labs(color='Treatment',  y='Power',x="Number of sampled labs",
                                                                      title='Lab random effect on power to detect treatment effect',
                                                                      subtitle='Mean power of n=100 artificial datasets')

# Separate plots
ggplot(d, aes(Xax,value)) + 
  geom_point() + 
  stat_smooth() +
  facet_wrap(~variable)


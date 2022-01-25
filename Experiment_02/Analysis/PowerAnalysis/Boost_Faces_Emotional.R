
library("ggplot2")
library("lme4")
library("lm.beta")
library(Rcmdr)

rm(list=ls()) 

my_path = "/Users/tomsalomon/Drive/Experiment_Israel/Codes/Boost_faces_emotional/Analysis/PowerAnalysis/"
my_datafile = "SummaryTable.txt"
PHQ_data=read.table((paste(my_path,my_datafile,sep = "")),header=T)

# Plotting workspace
quartz(width=5,height=4,pointsize=12,dpi=100)

# PHQ ~ PreferenceBias + CAT_effect
model=lm(PHQ_rank ~ PreferenceBias_rank + CAT_effect_rank,data=PHQ_data)
PHQ_data$Prediction=model$fitted.values
summary(lm.beta(lm(PHQ_rank ~ PreferenceBias_rank + CAT_effect_rank,data=PHQ_data)))

# Prediction versus actual values plot
ggplot(data=PHQ_data,aes(x=Prediction,y=PHQ,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_smooth(method = "lm") + 
  xlab("Model Predicion - PHQ Rank") +
  ylab("PHQ (Rank)") +
  theme_bw()

# PHQ ~ PreferenceBias
summary(lm.beta(lm(PHQ_rank ~ PreferenceBias_rank,data=PHQ_data)))

# Raw data plot
ggplot(data=PHQ_data,aes(x=PreferenceBias,y=PHQ,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_smooth(method = "lm") + 
  xlab("Preferences Bias") +
  geom_vline(xintercept=0,linetype=2) +
  theme_bw()

# Ranks plot
ggplot(data=PHQ_data,aes(x=PreferenceBias_rank,y=PHQ_rank,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_smooth(method = "lm") + 
  xlab("Preferences Bias (rank)") +
  ylab("PHQ (rank)") +
  theme_bw()

# PHQ ~ CAT_effect
summary(lm.beta(lm(PHQ_rank ~ CAT_effect_rank,data=PHQ_data)))

# Raw data plot
ggplot(data=PHQ_data,aes(x=CAT_effect,y=PHQ,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_smooth(method = "lm") + 
  xlab("Preferences Bias") +
  geom_vline(xintercept=0.5,linetype=2) +
  theme_bw()

# Ranks plot
ggplot(data=PHQ_data,aes(x=CAT_effect_rank,y=PHQ_rank,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_smooth(method = "lm") + 
  xlab("CAT Effect (rank)") +
  ylab("PHQ (rank)") +
  theme_bw()

# 3D Plot
scatter3d(PHQ_data$PreferenceBias_rank, PHQ_data$PHQ_rank, PHQ_data$CAT_effect_rank, size=10)
attach(PHQ_data)
scatter3d(PHQ_rank ~ PreferenceBias_rank + CAT_effect_rank, data = PHQ_data, size=10)
scatter3d(PHQ_rank ~ SubjectID + CAT_effect_rank, data = PHQ_data, size=10)

PHQ_data2=rbind(PHQ_data,PHQ_data[1:10,])
PHQ_data2$nullmodel= 0
PHQ_data2$nullmodel[1:length(PHQ_data[,1])]=1
PHQ_data2$nullmodel=as.factor(PHQ_data2$nullmodel)
PHQ_data2$PHQ_rank[PHQ_data2$nullmodel==0]=mean(PHQ_data$PHQ_rank)
PHQ_data2$PreferenceBias_rank[PHQ_data2$nullmodel==0]=mean(PHQ_data$CAT_effect_rank)
PHQ_data2$CAT_effect_rank[PHQ_data2$nullmodel==0]=mean(PHQ_data$CAT_effect_rank)

PHQ_data2$z_PHQ_rank[PHQ_data2$nullmodel==0]=0
PHQ_data2$z_PreferenceBias_rank[PHQ_data2$nullmodel==0]=0
PHQ_data2$z_CAT_effect_rank[PHQ_data2$nullmodel==0]=0

scatter3d(z_PHQ_rank ~ z_PreferenceBias_rank + z_CAT_effect_rank | nullmodel, data = PHQ_data2, 
          parallel=FALSE, fit="linear",
          surface.col=list("black","blue") ,surface.alpha = list(0.3),axis.scales=FALSE,
          xlab = c("Preference Bias (rank)"), 
          ylab = c("PHQ (rank)"), 
          zlab = c("CAT effect (rank)"), 
          fov=90) 


# 3D animation 
for (i in seq(1,360,2)) {
  rgl.viewpoint(i, 15)
  # Save image
  snapshot3d(paste(my_path,1000+i,".png",sep=""),fmt="png")
}

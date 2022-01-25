rm(list=ls())

library("reshape")
library("lm.beta")
library("ggplot2")

## Original Sample
#For iMac
path="~/Drive/Experiment_Israel/Codes/Boost_faces_emotional/Output/"
# For PC
#path="  "
subjects=c(113,132,301:312,314:326,328:331,333:346,348:359) # Technical only.
subjects=c(113,301,303:304,306:307,309:312,314:316,318:326,328:331,333:346,348:359) # Define here your subjects' codes.

#exclude (group):
# 308 (1), 317 (2) - Training: false alarm
# 132/332 (2), 302 (1), 305 (1) - Training: minimal ladder; 132 - also Training: misses.

# 356 (PHQ=7), 359 (PHQ=6) - can't assign group

# not really excluded
# 313 (2) - was coded as 113 - not excluded
# 332 (2) - was coded as 132. Training: minimal ladder (stopped hearing cue).
# 327 (2) - did not run subject with that code
# 347 (1) - did not complete the experiment


filelist=c()
for (s in subjects){
  # These files were created with the MATAB code: 'CreateItemRanking.m'
  filelist=c(filelist,Sys.glob(paste(path, "*",s,"*ItemRankingResults_combined*.txt",sep="")))
}
Data=c()
for (s in subjects){
  tmp_file=Sys.glob(paste(path, "*",s,"*ItemRankingResults_combined*.txt",sep=""))
  tmp_data=read.table(tmp_file,header=T,na.strings=c(999,999000))
  tmp_data$SubjectID=s
  tmp_data$Colley_score_ranked=rank(tmp_data$Rank)
  Data=rbind(Data,tmp_data)
}

Data$SubjectID[Data$SubjectID<300]=Data$SubjectID[Data$SubjectID<300]+200 # Conver subjects who were misscoded as 1XX to 3XX
Data$Face_Affect[Data$StimNum<=40]=1
Data$Face_Affect[Data$StimNum>40]=2
Data$Face_Affect=factor(Data$Face_Affect,labels=c("Happy","Neutral"))

PHQ_data=read.table(Sys.glob(paste(path, "/PHQ_Scores.txt",sep="")),header=T)
PHQ_data$SubjectID_name=as.character(PHQ_data$SubjectID)
PHQ_data$SubjectID_num=as.numeric(substr(PHQ_data$SubjectID_name,nchar(PHQ_data$SubjectID_name)-2,nchar(PHQ_data$SubjectID_name)))
PHQ_data$SubjectID=PHQ_data$SubjectID_num

Data=merge(Data,PHQ_data[c("SubjectID","PHQ")],by = "SubjectID")
model_data=melt(tapply(Data$Rank, list(SubjectID=Data$SubjectID,PHQ=Data$PHQ,Face_Affect=Data$Face_Affect), mean,na.rm=T),id="Subject")
model_data=model_data[!is.na(model_data$value),]
model_data=model_data[model_data$Face_Affect=="Happy",]
model_data$PHQ_ordinal=rank(model_data$PHQ)
model_data$value_ordinal=rank(model_data$value)

model_linear=lm(PHQ ~ value,data=model_data)
summary(lm.beta(model_linear))

dev.new()
ggplot(data=model_data,aes(value,PHQ,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_vline(xintercept=0.5,linetype=2) +
  geom_smooth(method='lm') +
  xlab("Initial Preferences - Bias towards Happy Faces") +
  theme_bw()

model_spearman=lm(PHQ_ordinal ~ value_ordinal,data=model_data)
summary(lm.beta(model_spearman))

dev.new()
ggplot(data=model_data,aes(value_ordinal,PHQ_ordinal,color=PHQ)) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  geom_smooth(method='lm') +
  xlab("Initial Preferences - Bias towards Happy Faces (Rank)") +
  ylab("PHQ (Rank)") +
  theme_bw()





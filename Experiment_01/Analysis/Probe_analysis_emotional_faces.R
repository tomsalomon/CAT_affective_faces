
library(lme4)
library("ggplot2")
library(rockchalk)
library("lm.beta")
library("e1071")

rm(list=ls())

# Standard error function
se = function(x) { out=sqrt((var(x, na.rm = TRUE))/(length(which(!is.na(x))))) }

## Original Sample
#For iMac
path="~/Drive/Experiment_Israel/Codes/Boost_faces_emotional_II/Output/"
# For PC
#path="  "

#subjects_group1=c(301,303,307,310,314:316,323,326,328,334:336,338,340,342,345,346,351:355) # Define here your subjects' codes.
#subjects_group2=c(113,304,306,309,311:312,318:322,324:325,329:331,333,337,339,341,343,344,348:350,356:359)
subjects=c(101:102,104:114,116:118,120:122,124:125,127:139,142:152,155:157,159:166,168:178) # Technical only.



subjects_group1 = c(101, 103, 104, 106, 110, 112, 115, 120, 121, 123, 124, 126, 129, 132, 135, 137, 142, 146, 147, 149, 151, 155, 156, 158, 162, 163, 165, 167, 169, 170, 173)
subjects_group2 = subjects[!is.element(subjects,subjects_group1)]
#exclude (group):
# 103, 123, 126, 140, 154, 167 - Training: false alarm
# 158 - Training: misses
# 115, 158 - Training: minimal ladder
# 119, 141, 153 - quit in the middle

# Load in the data - Group I
filelist=c()
for (s in subjects_group1){
	filelist=c(filelist,Sys.glob(paste(path, "/*",s,"*probe_block*.txt",sep="")))
}
Data_group1=c()
for (f in filelist){
  Data_group1=rbind(Data_group1,read.table(f,header=T,na.strings=c(999,999000)))
}

# Load in the data - Group II
filelist=c()
for (s in subjects_group2){
  filelist=c(filelist,Sys.glob(paste(path, "/*",s,"*probe_block*.txt",sep="")))
}
Data_group2=c()
for (f in filelist){
  Data_group2=rbind(Data_group2,read.table(f,header=T,na.strings=c(999,999000)))
}

# Merge both groups
Data_group1$group=1
Data_group2$group=2
Results=rbind(Data_group1,Data_group2)
Results$SubNumeric=as.numeric(gsub("[^0-9]", "",Results$subjectID)) # Convert SubjectID into numbers only
Results$Order=as.factor(2-Results$SubNumeric%%2)
Results$group=as.factor(Results$group)

Results$Outcome[Results$Outcome>1]=NA # transform no choice (999) and non Go-NoGo (99) outcome into NA
Results$PairType2[Results$PairType==1]="Happy_Faces"
Results$PairType2[Results$PairType==2]="Neutral_Faces"
Results$PairType2[Results$PairType==3]="Sanity_Happy"
Results$PairType2[Results$PairType==4]="Sanity_Neutral"

tapply(Results$Outcome,Results$PairType2,mean,na.rm=T)
tapply(Results$Outcome, list(group=Results$group,PairType2=Results$PairType2), mean,na.rm=T)
tapply(Results$Outcome, list(SubjectID=Results$subjectID,PairType2=Results$PairType2,group=Results$group), mean,na.rm=T)

# Main analysis - logistic regression
# Group 1: Happy Faces
summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==1 & group==1 )),na.action=na.omit,family=binomial)) 
# Group 1: Neutral Faces
summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==2 & group==1 )),na.action=na.omit,family=binomial)) 
# Group 2: Happy Faces
summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==1 & group==2 )),na.action=na.omit,family=binomial)) 
# Group 2: Neutral Faces
summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==2 & group==2 )),na.action=na.omit,family=binomial)) 

# Group 1: Happy-Neutral difference
summary(glmer(Outcome ~ 1 + PairType + (1|subjectID),data=subset(Results,(Results$PairType<=2 & group==1 )),na.action=na.omit,family=binomial)) 
# Group 2: Happy-Neutral difference
summary(glmer(Outcome ~ 1 + PairType + (1|subjectID),data=subset(Results,(Results$PairType<=2 & group==2 )),na.action=na.omit,family=binomial)) 

# Full model - Interaction: group and stimulus
summary(glmer(Outcome ~ 1 + PairType * group + (1|subjectID),data=subset(Results,(Results$PairType<=2)),na.action=na.omit,family=binomial)) 
# Interaction: group and stimulus in sanity
summary(glmer(Outcome ~ 1 + PairType + group + (1|subjectID),data=subset(Results,(Results$PairType==3|Results$PairType==4)),na.action=na.omit,family=binomial)) 

# Happy faces: group differences
summary(glmer(Outcome ~ 1 + group + (1|subjectID),data=subset(Results,(Results$PairType==1)),na.action=na.omit,family=binomial)) 
# Neutral faces: group differences
summary(glmer(Outcome ~ 1 + group + (1|subjectID),data=subset(Results,(Results$PairType==2)),na.action=na.omit,family=binomial)) 
# Sanity : group differences
summary(glmer(Outcome ~ 1 + group + (1|subjectID),data=subset(Results,(Results$PairType==3|Results$PairType==4)),na.action=na.omit,family=binomial)) 

# PHQ prediction
PHQ_data=read.table(Sys.glob(paste(path, "/PHQ_Scores.txt",sep="")),header=T)
PHQ_data$SubjectID_name=as.character(PHQ_data$SubjectID)
PHQ_data$SubjectID_num=as.numeric(substr(PHQ_data$SubjectID_name,nchar(PHQ_data$SubjectID_name)-2,nchar(PHQ_data$SubjectID_name)))

data_by_subject_group=as.data.frame(tapply(Results$Outcome, list(SubjectID=Results$subjectID,PairType2=Results$PairType2,group=Results$group), mean,na.rm=T))
data_by_subject_group$group[data_by_subject_group$Happy_Faces.1>=0]=1
data_by_subject_group$group[data_by_subject_group$Happy_Faces.2>=0]=2
data_by_subject_group$Happy_Faces=0
data_by_subject_group$Happy_Faces[data_by_subject_group$group==1]=data_by_subject_group$Happy_Faces.1[data_by_subject_group$group==1]
data_by_subject_group$Happy_Faces[data_by_subject_group$group==2]=data_by_subject_group$Happy_Faces.2[data_by_subject_group$group==2]
data_by_subject_group$Neutral_Faces[data_by_subject_group$group==1]=data_by_subject_group$Neutral_Faces.1[data_by_subject_group$group==1]
data_by_subject_group$Neutral_Faces[data_by_subject_group$group==2]=data_by_subject_group$Neutral_Faces.2[data_by_subject_group$group==2]
data_by_subject_group$group=as.factor(data_by_subject_group$group)
data_by_subject_group$SubjectID_name=as.character(row.names(data_by_subject_group))
data_by_subject_group$SubjectID_num=as.numeric(substr(data_by_subject_group$SubjectID_name,nchar(data_by_subject_group$SubjectID_name)-2,nchar(data_by_subject_group$SubjectID_name)))
# convert miscalssified id's of 100's as 300'
#data_by_subject_group$SubjectID_num[data_by_subject_group$SubjectID_num<200]=data_by_subject_group$SubjectID_num[data_by_subject_group$SubjectID_num<200]+200
data_by_subject_group=merge(data_by_subject_group, PHQ_data, by="SubjectID_num")
data_by_subject_group$diff=data_by_subject_group$Happy_Faces - data_by_subject_group$Neutral_Faces
data_by_subject_group$PHQ_binary=1
data_by_subject_group$PHQ_binary[data_by_subject_group$PHQ<9]=-1
data_by_subject_group$PHQ_binary=as.factor(data_by_subject_group$PHQ_binary)
levels(data_by_subject_group$PHQ_binary)=c("Low PHQ","High PHQ")
data_by_subject_group$diff_binary[data_by_subject_group$diff>0]="Delta > 0"
data_by_subject_group$diff_binary[data_by_subject_group$diff<=0]="Delta <= 0"

# plot PHQ distribution
dev.new()
ggplot(data_by_subject_group,aes(PHQ_binary,PHQ,color=PHQ_binary)) +
  geom_point(shape = 16, size = 4, alpha= 0.6,position = position_jitter(w = 0.2, h = 0)) +
  geom_hline(yintercept = 9, linetype=2) +
  geom_hline(yintercept = 5, linetype=2) +
  xlab("Experimental group") +
  theme_bw()


# plot PHQ distribution
dev.new()
ggplot(data_by_subject_group,aes(PHQ,fill=PHQ_binary,color=PHQ_binary)) +
  geom_histogram(breaks=seq(-1,27,1),alpha=0.6, position="identity", lwd=0.2)+
  xlab("PHQ") +
  theme_bw()

descriptives_all=setNames(data.frame(tapply(Results$Outcome,Results$PairType2,mean,na.rm=T)),'means')
descriptives_all$pairtype=as.factor(rownames(descriptives_all))
descriptives_all$SEM=c(se(data_by_subject_group$Happy_Faces),se(data_by_subject_group$Neutral_Faces),
                       se(c(data_by_subject_group$Sanity_Happy.1,data_by_subject_group$Sanity_Happy.2)),
                       se(c(data_by_subject_group$Sanity_Neutral.1,data_by_subject_group$Sanity_Neutral.2)))
descriptives_all$asteriks="***"
# Bar Plot - Proportions all
dev.new()
ggplot(descriptives_all, aes(pairtype,means)) +
  #theme_bw() + # white background
  geom_bar(width=.7,position=position_dodge(0.7), stat="identity") + # Bar plot
  theme(legend.position="top",legend.title=element_blank()) + # position legend
  geom_errorbar(position=position_dodge(.7), width=.7/4, aes(ymin=means-SEM, ymax=means+SEM))  + # add error bar of SEM
  scale_y_continuous("Proportion of trials Go items were chosen",limit=c(0,1),breaks=seq(0, 1, 0.1),expand = c(0,0)) + # define y axis properties
  geom_abline(intercept = (0.5),slope=0,linetype =2, size = 1,show.legend=TRUE,aes()) + # chace level 50% reference line
  geom_text(position=position_dodge(.7),aes(y=means+SEM+0.01,label=(asteriks)),size=8) # significance asteriks

descriptives_by_group=data.frame(melt(tapply(Results$Outcome, list(group=Results$group,PairType2=Results$PairType2), mean,na.rm=T),id="group"))
descriptives_by_group$group=as.factor(descriptives_by_group$group)
levels(descriptives_by_group$group)=c("Low PHQ","High PHQ")
colnames(descriptives_by_group)[2]="pairtype"
colnames(descriptives_by_group)[3]="means"
descriptives_by_group$SEM=c(se(data_by_subject_group$Happy_Faces.1),
                            se(data_by_subject_group$Happy_Faces.2),
                            se(data_by_subject_group$Neutral_Faces.1),
                            se(data_by_subject_group$Neutral_Faces.2),
                            se(data_by_subject_group$Sanity_Happy.1),
                            se(data_by_subject_group$Sanity_Happy.2),
                            se(data_by_subject_group$Sanity_Neutral.1),
                            se(data_by_subject_group$Sanity_Neutral.2))
descriptives_by_group$p=0

m1=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==1 & group==1 )),na.action=na.omit,family=binomial)) # Group 1: Happy Faces
m2=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==1 & group==2 )),na.action=na.omit,family=binomial)) # Group 2: Happy Faces
m3=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==2 & group==1 )),na.action=na.omit,family=binomial)) # Group 1: Neutral Faces
m4=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(Results,(Results$PairType==2 & group==2 )),na.action=na.omit,family=binomial)) # Group 2: Neutral Faces

descriptives_by_group$p[1:4]=c(m1$coefficients[4],m2$coefficients[4],m3$coefficients[4],m4$coefficients[4])
descriptives_by_group$asteriks=""
descriptives_by_group$asteriks[descriptives_by_group$p<0.05]="*"
descriptives_by_group$asteriks[descriptives_by_group$p<0.01]="**"
descriptives_by_group$asteriks[descriptives_by_group$p<0.001]="***"

# Bar Plot - Proportions by group
dev.new()
ggplot(descriptives_by_group[1:4,], aes(x=pairtype,y=means,fill=group)) +
  theme_bw() + # white background
  geom_bar(width=.7,position=position_dodge(0.7), stat="identity") + # Bar plot
  # scale_fill_brewer(palette="Set2") +
  theme(legend.position="top",legend.title=element_blank()) + # position legend
  geom_errorbar(position=position_dodge(.7), width=.7/4, aes(ymin=means-SEM, ymax=means+SEM))  + # add error bar of SEM
  scale_y_continuous("Proportion of trials Go items were chosen",limit=c(0,1),breaks=seq(0, 1, 0.1),expand = c(0,0)) + # define y axis properties
  geom_abline(intercept = (0.5),slope=0,linetype =2, size = 1,show.legend=TRUE,aes()) + # chace level 50% reference line
 geom_text(position=position_dodge(.7),aes(y=means+SEM+0.01,label=(asteriks)),size=5) # significance asteriks

# Bar Plot - Proportions by group
dev.new()
ggplot(descriptives_by_group[5:8,], aes(x=pairtype,y=means,fill=group)) +
  theme_bw() + # white background
  geom_bar(width=.7,position=position_dodge(0.7), stat="identity") + # Bar plot
  theme(legend.position="top",legend.title=element_blank()) + # position legend
  geom_errorbar(position=position_dodge(.7), width=.7/4, aes(ymin=means-SEM, ymax=means+SEM))  + # add error bar of SEM
  scale_y_continuous("Proportion of trials Go items were chosen",limit=c(0,1),breaks=seq(0, 1, 0.1),expand = c(0,0)) + # define y axis properties
  geom_abline(intercept = (0.5),slope=0,linetype =2, size = 1,show.legend=TRUE,aes()) + # chace level 50% reference line
  geom_text(position=position_dodge(.7),aes(y=means+SEM+0.01,label=(asteriks)),size=5) # significance asteriks

 summarySE(tg, measurevar="len", groupvars=c("supp","dose"))


# Plot CAT effect for happy and neutral faces for each group
dev.new()
ggplot(data_by_subject_group,aes(Happy_Faces, Neutral_Faces,color=PHQ_binary))+
  xlim(0,1) +
  ylim(0,1) +
  xlab("CAT effect - Happy Faces") +
  ylab("CAT effect - Neutral Faces") +
  geom_hline(yintercept = .5) +
  geom_vline(xintercept = .5) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_point(shape = 16, size = 5, alpha= 0.6) +
  theme_bw()

# Plot PHQ versus delta effect
dev.new()
ggplot(data_by_subject_group,aes(diff, PHQ,color=PHQ_binary))+
  #  xlim(-1,1) +
  xlab(expression(Delta*" CAT effect (Happy Faces - Neutral Faces)")) +
  geom_hline(yintercept = 9, linetype=2) +
  geom_hline(yintercept = 5, linetype=2) +
  geom_vline(xintercept = 0) +
  geom_point(shape = 16, size = 5, alpha= 0.6) +
  theme_bw()

PHQ_CAT_table=table(data_by_subject_group$PHQ_binary,data_by_subject_group$diff_binary)
chisq.test(PHQ_CAT_table)

PHQ_pred_model <- lm(PHQ ~ Happy_Faces + Neutral_Faces, data = data_by_subject_group)
dev.new()
plotPlane(PHQ_pred_model, plotx1 = "Happy_Faces", plotx2 = "Neutral_Faces", drawArrows = TRUE, pcol="PHQ")


cor(data_by_subject_group$Happy_Faces.1,data_by_subject_group$Neutral_Faces.1,use = "pairwise.complete.obs")
cor(data_by_subject_group$Happy_Faces.2,data_by_subject_group$Neutral_Faces.2,use = "pairwise.complete.obs")

cor(data_by_subject_group$Happy_Faces,data_by_subject_group$PHQ,use = "pairwise.complete.obs")
cor(data_by_subject_group$Neutral_Faces,data_by_subject_group$PHQ,use = "pairwise.complete.obs")
cor(data_by_subject_group$diff,data_by_subject_group$PHQ,use = "pairwise.complete.obs")

data_by_subject_group$PHQ_rank=rank(data_by_subject_group$PHQ)
data_by_subject_group$Happy_Faces_rank=rank(data_by_subject_group$Happy_Faces)
data_by_subject_group$Neutral_Faces_rank=rank(data_by_subject_group$Neutral_Faces)

model = lm(PHQ ~ 1 + Happy_Faces + Neutral_Faces,data=data_by_subject_group,na.action=na.omit)
model_logistic = glm(PHQ_binary ~ Happy_Faces * Neutral_Faces,data=data_by_subject_group,na.action=na.omit,family = "binomial")
data_by_subject_group$prediction=predict(model_logistic)
data_by_subject_group$prediction_prop=exp(data_by_subject_group$prediction)/(1+exp(data_by_subject_group$prediction))
anova(model, model_logistic, test="Chisq")
data_by_subject_group$predicted_PHQ=predict(model)
model_rank = lm(PHQ_rank ~ 1 + Happy_Faces_rank ,data=data_by_subject_group,na.action=na.omit)
model_rank_full = lm(PHQ_rank ~ 1 + Happy_Faces_rank*Neutral_Faces_rank ,data=data_by_subject_group,na.action=na.omit)

anova(model_rank,model_rank_full)

summary (model)
summary(lm.beta(model_rank))

# Plot PHQ versus predicted PHQ (linear model)
dev.new()
ggplot(data_by_subject_group,aes(Happy_Faces_rank, PHQ_rank,color=PHQ))+
  xlab("CAT Effect: Happy Faces (Rank)") +
  ylab("PHQ (Rank)") +
  geom_smooth(method='lm') +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  theme_bw()

dev.new()
ggplot(data_by_subject_group,aes(Happy_Faces, PHQ,color=PHQ))+
  xlab("CAT Effect: Happy Faces") +
  ylab("PHQ") +
  geom_smooth(method='lm') +
  geom_vline(xintercept = 0.5, linetype=2) +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  theme_bw()

# Plot PHQ versus predicted PHQ (binary)
dev.new()
ggplot(data_by_subject_group,aes(prediction_prop, PHQ,color=PHQ_binary))+
  xlim(0,1) +
  xlab(expression("Model prediction - probability of high PHQ")) +
  geom_hline(yintercept = 9, linetype=2) +
  geom_hline(yintercept = 5, linetype=2) +
  geom_vline(xintercept = 0.5, linetype=3) +
  geom_point(shape = 16, size = 5, alpha= 0.6) +
  theme_bw()

PHQ_CAT_table2=table(data_by_subject_group$PHQ_binary,data_by_subject_group$prediction_prop_bin)
chisq.test(PHQ_CAT_table)

subjects=c(subjects_group1,subjects_group2)
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

#Data$SubjectID[Data$SubjectID<300]=Data$SubjectID[Data$SubjectID<300]+200 # Conver subjects who were misscoded as 1XX to 3XX
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

model_data$SubjectID_num=model_data$SubjectID
model_data=merge(model_data,data_by_subject_group[,c("SubjectID_num","Happy_Faces","Neutral_Faces","Happy_Faces_rank","Neutral_Faces_rank")],by="SubjectID_num")
names_tmp=colnames(model_data)
names_tmp[names_tmp=="value"]="Happy_Faces_initial_prefernces"
names_tmp[names_tmp=="value_ordinal"]="Happy_Faces_initial_prefernces_ordinal"
names_tmp[names_tmp=="Happy_Faces"]="Happy_Faces_CAT_effect"
names_tmp[names_tmp=="Happy_Faces_rank"]="Happy_Faces_CAT_effect_ordinal"

colnames(model_data)=names_tmp

model_PHQ_by_CAT_and_prefernces=lm(PHQ_ordinal ~ Happy_Faces_CAT_effect_ordinal + Happy_Faces_initial_prefernces_ordinal,data=model_data)
model_PHQ_by_CAT=lm(PHQ_ordinal ~ Happy_Faces_CAT_effect_ordinal ,data=model_data)
model_PHQ_by_prefernces=lm(PHQ_ordinal ~ Happy_Faces_initial_prefernces_ordinal,data=model_data)

summary(lm.beta(model_PHQ_by_CAT_and_prefernces))
anova(model_PHQ_by_CAT_and_prefernces,model_PHQ_by_prefernces)

# Plot model prediction versus actual values
model_data$predicted=predict(model_PHQ_by_CAT_and_prefernces, model_data[,c("Happy_Faces_CAT_effect_ordinal","Happy_Faces_initial_prefernces_ordinal")])
dev.new()
ggplot(model_data,aes(predicted,PHQ_ordinal,color=PHQ))+
  scale_colour_gradient(low = "red", high = "yellow") +
  xlab(expression("Model prediction - PHQ Rank")) +
  ylab(expression("Actual PHQ Rank")) +
  geom_smooth(method='lm',color="red") +
  geom_point(shape = 16, size = 5, alpha= 0.8) +
  theme_bw()

# 3D plot
colors=model_data$PHQ_ordinal
colors=rank(model_data$Happy_Faces_initial_prefernces_ordinal*model_data$Happy_Faces_CAT_effect_ordinal)
colfunc<-colorRampPalette(c("red","yellow"))
colormap=colfunc(length(unique(colors)))
colors=colormap[c(colors)]
dev.new()
s3d <- scatterplot3d(model_data$Happy_Faces_initial_prefernces_ordinal, model_data$Happy_Faces_CAT_effect_ordinal, model_data$PHQ_ordinal, type = "h", lab = c(3, 3),
                     pch=19,color=colors,
                     xlab="Initial Preference Bias towards Happy Faces (Rank)",
                     ylab="CAT Effect - Happy Faces (Rank)",
                     zlab="PHQ (Rank)") # lab: number of tickmarks on x-/y-axes
s3d$plane3d(model_PHQ_by_CAT_and_prefernces)

# # 3D animation 
#  for (i in seq(1,360,2)) {
# 
#    #saves the plot as a .png file in the working directory
#  png(paste(1000+i,".png",sep=""))
#    s3d <- scatterplot3d(model_data$Happy_Faces_initial_prefernces_ordinal, model_data$Happy_Faces_CAT_effect_ordinal, model_data$PHQ_ordinal, type = "h", lab = c(3, 3),
#                         pch=19,color=colors,
#                         xlab="Initial Preference Bias towards Happy Faces (Rank)",
#                         ylab="CAT Effect - Happy Faces (Rank)",
#                         zlab="PHQ (Rank)",
#                         angle=i) # lab: number of tickmarks on x-/y-axes
#    s3d$plane3d(model_PHQ_by_CAT_and_prefernces)
#   dev.off()
#  }
# my_command <- 'convert *.png -delay 1 -loop 0 3d.gif'
# system(my_command)

 # ANIMATE: https://davetang.org/muse/2015/02/12/animated-plots-using-r/

#svm_dat=data_by_subject_group[ , c("Happy_Faces","Neutral_Faces","PHQ_binary")]
#svm_fit=svm(PHQ_binary ~ .,data=svm_dat,kernel="linear",cost= 1e-13,scale = FALSE)
#tuned=tune(svm, PHQ_binary ~ .,data = svm_dat, kernel="linear",ranges=list(cost=c(1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1e-00, 1e01, 1e02)),scale=FALSE)
#summary(tuned)
#print(svm_fit)
#plot(svm_fit,svm_dat[,c("Happy_Faces","Neutral_Faces","PHQ_binary")])
#plot(svm_fit,svm_dat)


